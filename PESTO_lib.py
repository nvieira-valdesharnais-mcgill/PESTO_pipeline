import sys
import os
from subprocess import run
from astropy.io import fits
import numpy as np
import time

class PESTO_data:
     """
     Input: locations_list: A list of strings, one for each directory with PESTO data
            image_type: A list whose elements are the strings 'calibration' or 'object', 
                           describing the content of the files in each location in locations_list
            name: A name to use in the creation of new files
            Example input: (['path/to/object','path/to/flat'],['object','calibration'], foo)
                           (['homes/huanca/aramirez/SkyFlat','homes/huanca/aramirez/GRB180618A'],['calibration','object'], bar)
     Return: PESTO_data object
     """
     def __init__(self, locations_list,image_type_list,name):
          self.loc = locations_list
          self.imtype = dict([(locations_list[i],image_type_list[i]) for i in range(0,len(locations_list))])
          self.name = name
          self.r_obj = []
          self.r_cal = []
          self.g_obj = []
          self.g_cal = []
          self.i_obj = []
          self.i_cal = []
          self.z_obj = []
          self.z_cal = []
          self.NA_cal = []
          self.bias = []
          self.hdr_ind = {} # ADDED BY NICK: dict with location and which header index to use in said location
 
     def untar(self):
          """
          Input: None
          Return: None      
          Decompresses all the (compressed) files associated with
          the PESTO_data object.
          """
          for l in self.loc:
               files = os.listdir(l)
               for f in files:
                    if '.gz' in f:
                         run(['gzip', '-d','-f',l+'/'+f])

     def tar(self):
         """
         Input: None
         Return: None
         Compresses all the (decompressed) files associated with
         the PESTO_data object.
         """
         for l in self.loc:
              files = os.listdir(l)
              for f in files:
                   if '.gz' not in f:
                        run(['gzip','-f',l+'/'+f])

     def flush(self): 
         """
         Input:None
         Output:None
         Flushes the working directory (but doesn't delete it) of unneeded files
         and optionally deletes the segmented image and csv. 
         Added by Nick and Val.
         """
         temp = os.getcwd()
         yes_no = 'x'
         valid_inputs = ["YES","NO","Y","N","y","n","yes","no","Yes","No"]
         while not(yes_no in valid_inputs) and self.reduced:           
             yes_no = input('Delete the segmented image and csv, too? [y/n]')
         if self.reduced:
             os.chdir(self.loc[0]+'/'+self.name+'/..')
             run('rm -rf '+self.name, shell=True)
         else:
             os.chdir(self.tgt+'/'+self.name)
         #run('rm -rf *', shell=True)    
         run('rm -rf *.py', shell=True)
         run('rm -rf *.fits', shell=True)
         run('rm -rf *.txt', shell=True)
         run('rm -rf *.par', shell=True)
         run('rm -rf logfile', shell=True)
         run('rm -rf minired', shell=True)    # HARD FIX
         os.chdir(temp)
         if yes_no in ["YES","Y","y","yes","Yes"]:
             run('rm -rf *.png', shell=True)
             run('rm -rf *.csv', shell=True)
         run('rm -rf minired', shell=True)
         print("\nFlush!")

##########################

class raw_PESTO_data(PESTO_data):
     
    def __init__(self,locations_list,image_type_list,name, target_directory):
        super(raw_PESTO_data, self).__init__(locations_list,image_type_list,name)
        self.tgt = target_directory
        self.reduced = False
        self.list_made = False
        self.workdir_present = False

        # Allows new data to be read regardless of header indexing.
        # Assumption: only need to look at one .fits per location, and if extended,
        # only need to look at the first extension.
        for l in self.loc:
            files = os.listdir(l)
            for f in files:
                if '.fits' in f: # Look until a .fits is found
                    tempfile = f
                hdr_temp = fits.open(l+'/'+tempfile)
                
                # First: case of an unextended file
                if ('FILTRE' in hdr_temp[0].header): # If 0th header contains 'FILTRE' keyword
                    self.hdr_ind[l] = 0 # Use index 0 for all header reading in this location
                # Second: case of an extended file
                elif ('FILTRE' in hdr_temp[1].header):  
                    self.hdr_ind[l] = 1  
           
    def extended_header_cleanup(self):
        """
        Input:None
        Output:None
        In certain extended .fits files, the SIMPLE, EXTEND, and NEXTEND headers 
        are misplaced. This corrects that. For a given data object, only needs to be run once.
        If the files in a given location have already been modified (i.e., you initialized 
        an identical object and used this function on the object at any point in time during 
        any python session), this does not need to run. 
        Added by Nick and Val.
        WARNING: This PERMANENTLY modifies data in a given directory.
        """
        for l in self.loc:
            files = os.listdir(l)
            for f in files: 
                if '.fits' in f:
                    if self.hdr_ind[l] != 0: # i.e., looking at extended .fits files 
                       n_extend = fits.open(l+'/'+f)[0].header['NEXTEND'] # Get the no. of extensions
                       hdr_temp = fits.open(l+'/'+f,mode='update')
                       for n in range(1, n_extend+1):
                           hdu = hdr_temp[n]
                           if 'SIMPLE' in hdu.header:
                               temp_simple = hdu.header['SIMPLE']
                               hdu.header.remove('SIMPLE')  # Shifts all other cards up by 1
                               hdu.header.append(('SIMPLE',temp_simple)) # Append header to end (just in case)
                           if 'EXTEND' in hdu.header:
                               hdu.header.remove('EXTEND')   
                               hdu.header.append(('EXTEND','T')) # Shouldn't be needed, but to be safe
                           if 'NEXTEND' in hdu.header:
                               hdu.header.remove('NEXTEND') # Shifts
                               hdu.header.append(('NEXTEND',n_extend)) # Be safe
                hdr_temp.close()

    def resize(self, sf, df,  mode, sizes):
        """
        Input: source filepath, destination filepath, trim or extend mode, list of pixel dimensions [left, right, top, bottom]
        Output: None
        Resizes the image to the appropriate dimensions by either trimming or extending by NaNs the pixel dimensions passed.
        Written by Val and Nick.
        WARNING: This method may PERMANENTLY ALTER/OVERWRITE DATA in both source and destination folders. 
        """
        if '.fits' in sf:
            hdr = fits.open(sf,mode='update')
            n_extend = range(1)
            if not 'FILTRE' in hdr[0].header:
                print("Extended file detected.")
                fits.append(df,hdr[0].data,hdr[0].header)
                n_extend = range(1,hdr[0].header['NEXTEND'])
            for n in n_extend:                
                image_header = hdr[n].header
                image_data = hdr[n].data
                if 'NAXIS1' in image_header and 'NAXIS2' in image_header:
                    x = image_header['NAXIS1']
                    y = image_header['NAXIS2']
                    if mode == "trim":
                        if (not all(s>=0 for s in sizes)) or (sizes[0]+sizes[1])>x or (sizes[2]+sizes[3])>y:
                            print("Dimensions to trim are out of range for: "+sf)
                            return
                        new_image_data = image_data[(sizes[2]):(y-sizes[3]), (sizes[0]):(x-sizes[1])]
                        image_header["NAXIS1"]=x-(sizes[0]+sizes[1])
                        image_header["NAXIS2"]=y-(sizes[2]+sizes[3])
                        fits.append(df,new_image_data,image_header)
                    elif mode == "extend":
                        if not all(s>=0 for s in sizes):
                            print("Dimensions to extend by are out of range for: "+sf)
                            return
                        new_image_data = np.hstack((np.full((y,sizes[0]),np.nan), image_data, np.full((y,sizes[1]),np.nan)))
                        new_image_data = np.vstack((np.full((sizes[2],(x+sizes[0]+sizes[1])),np.nan), new_image_data, np.full((sizes[3],(x+sizes[0]+sizes[1])),np.nan)))
                        image_header["NAXIS1"]=x+(sizes[0]+sizes[1])
                        image_header["NAXIS2"]=y+(sizes[2]+sizes[3])
                        fits.append(df,new_image_data,image_header)
                    else:
                        print("Incorrect mode argument passed. Please use 'trim' or 'extend'.")                               
                else:
                    print("Missing image dimensions (NAXIS1,NAXIS2) in header["+str(n)+"] for: "+sf)
                    return
            hdr.close()

    def resize_folder(self, source, destination_name, mode, sizes):
        """
        Input: source folder path, destination folder name, trim or extend mode, list of pixel dimensions [left, right, top, bottom]
        Output: None
        Resizes all fit files in a folder to the appropriate dimensions by either trimming or extending by NaNs the pixel dimensions passed and deposits in subfolder.
        Written by Val and Nick.
        WARNING: This method may PERMANENTLY ALTER/OVERWRITE DATA in both source and destination folders. 
        """
        files = os.listdir(source)
        if not os.path.exists(source+'/'+destination_name):
            os.mkdir(source+'/'+destination_name)
        for f in files:
            self.resize(source+'/'+f,source+'/'+destination_name+'/'+mode+'_'+f,mode,sizes)    
 
################################################################################

    def produce_lists(self,NA=False):
        """
        Input: None
        Output: None
        Produces lists of .fits files associated to the raw_PESTO_data object. 
        e.g all images of the object taken in the r band are added to self.r_obj,
        all calibration images in the r band are added to self.r_cal, etc.
        """
        if not self.list_made:
            for l in self.loc:			
                files = os.listdir(l)
                run_type = self.imtype[l]   
                for f in files:
                    if '.fits' in f:
                        hdr_temp = fits.open(l+'/'+f, mode = 'update')
                        hdu = hdr_temp[self.hdr_ind[l]] # Modified: was 0
                        hdu.header['NAXIS'] = 3

                        #Adds the type of frame (flat, bias, object, etc...) to the FITS header
                        #to help IRAF's data reduction.
                        if run_type == 'calibration':
                            # Check 2 conditions to see if it's a flat
                            isflat_condition1 = ('OBJECT' in hdu.header) and ('flat' in hdu.header['object'])
                            isflat_condition2 = ('IMAGETYP' in hdu.header) and ('flat' in hdu.header['imagetyp'])
                            isflat = isflat_condition1 or isflat_condition2
                            #To produce a list of biases
                            if ('IMAGETYP' in hdu.header) and ('zero' in hdu.header['imagetyp']):
                                self.bias.append(f.replace('.gz',''))
                            #To produce lists of which flats are in which band
                            elif isflat: 
                                if 'r' in hdu.header['filtre']:
                                    self.r_cal.append(f.replace('.gz',''))
                                elif 'g' in hdu.header['filtre']:
                                    self.g_cal.append(f.replace('.gz',''))
                                elif 'i' in hdu.header['filtre']:
                                    self.i_cal.append(f.replace('.gz',''))
                                elif 'z' in hdu.header['filtre']:
                                    self.z_cal.append(f.replace('.gz',''))
                                elif 'N/A' in hdu.header['filtre']:
                                    if NA:
                                        self.r_cal.append(f.replace('.gz',''))
                                        self.g_cal.append(f.replace('.gz',''))
                                        self.i_cal.append(f.replace('.gz',''))
                                        self.z_cal.append(f.replace('.gz',''))
                                    self.NA_cal.append(f.replace('.gz',''))

                        elif run_type == 'object':
                            hdu.header['imagetyp'] = 'object'
                            if 'r' in hdu.header['filtre']:
                                self.r_obj.append(f.replace('.gz',''))
                            elif 'g' in hdu.header['filtre']:
                                self.g_obj.append(f.replace('.gz',''))
                            elif 'i' in hdu.header['filtre']:
                                self.i_obj.append(f.replace('.gz',''))
                            elif 'z' in hdu.header['filtre']:
                                self.z_obj.append(f.replace('.gz',''))


                        hdr_temp.close()

                self.list_made=True
    
    def make_working_directory(self):
        """
        Input: None
        Output: None
        Creates a directory containing all .fits files associated to the object and one .txt file
        for each list bound to the raw_PESTO_data object. The created directory contains everything PyRAF needs for the data reduction.   
        """
        if not self.workdir_present:
             l = self.tgt+'/'+self.name

             run(['mkdir','-p', l])
             run('cp -f ~/login.cl '+l, shell=True)
             run('mkdir -p '+l+'/uparm', shell=True)
             run('cp -f ~/uparm/* '+l+'/uparm', shell=True) 
             
             #Prints the list for each filter to the working directory
             np.savetxt(l+'/r_list.txt',self.r_cal, fmt = "%s")
             np.savetxt(l+'/g_list.txt', self.g_cal, fmt = "%s")
             np.savetxt(l+'/i_list.txt', self.i_cal, fmt = "%s")
             np.savetxt(l+'/z_list.txt', self.z_cal, fmt = "%s")
             np.savetxt(l+'/NA_list.txt', self.NA_cal, fmt = "%s")
             np.savetxt(l+'/bias_list.txt', self.bias, fmt = "%s")
             np.savetxt(l+'/object_list_r.txt', self.r_obj, fmt = "%s")
             np.savetxt(l+'/object_list_g.txt', self.g_obj, fmt = "%s")
             np.savetxt(l+'/object_list_i.txt', self.i_obj, fmt = "%s")
             np.savetxt(l+'/object_list_z.txt', self.z_obj, fmt = "%s")
             self.workdir_present=True

             for l in self.loc:
                 run(['rsync','-r','-p',l+'/',self.tgt+'/'+self.name]) # p option added by Nick: should sync without changing permissions
                                                                       # which would mean all copied files have all permissions
             run('chmod 777 '+self.tgt+'/'+self.name+'/*', shell=True) # Should give everything permissions
 
    def delete_working_directory(self):
         """
         Input: None
         Output: None
         Removes the working directory if it is present
         """
         if self.workdir_present:
              run(['rm', '-r', self.tgt+'/'+self.name])
              self.workdir_present=False

    def pyraf_reduction(self):
        """
        Input: None
        Output: None
        If th
        """
        if self.list_made==False:
            return 'Please make sure lists for each band were produced'
        if self.workdir_present==False:
            return 'Please make sure a working directory exists'
        temp = os.getcwd()
        run(['rsync','datared.py',self.tgt+'/'+self.name])
        os.chdir(self.tgt+'/'+self.name)
        run("bash -c 'source activate iraf27 && python2 datared.py && source deactivate'", shell=True)
        os.chdir(temp)

    def extract_reduced_images(self,reduced_data_folder_loc, reduced_data_folder_name):
        run(['mkdir', reduced_data_folder_loc+'/'+reduced_data_folder_name])
        run('cp -r '+self.tgt+'/'+self.name+'/*reduced.fits '+reduced_data_folder_loc+'/'+reduced_data_folder_name, shell=True) 
        return reduced_PESTO_data([reduced_data_folder_loc], ['object'], reduced_data_folder_name)

##########################

import string

class Del:
  def __init__(self, keep=string.digits):
    self.comp = dict((ord(c),c) for c in keep)
  def __getitem__(self, k):
    return self.comp.get(k)

keep_digits = Del()


#########################

class reduced_PESTO_data(PESTO_data):
     def __init__(self, locations_list, image_type_list, name):
        super(reduced_PESTO_data, self).__init__(locations_list, image_type_list, name)
        self.reduced = True


     def WCS_preparation(self, angle=0):
        """
        Input: The angle of the frame relative to a frame where X=RA, Y=DEC, counterclockwise, in degrees.
               Defaults to 0.
        Output: None. 
        Modifies the reduced fits file to provide it all the information the WCS conversion needs.
        """
        files = os.listdir(self.loc[0]+'/'+self.name)
        for f in files: 
             hdu_temp = fits.open(self.loc[0]+'/'+self.name+'/'+f,mode='update')[0]

             # Picking out the reference pixel (located around 80, 300 for MAXI in July/Sept data)
             #image_data = hdu_temp.data[0:95,300:350] # Clipping the top to avoid the horizontal dead pixel
             #brightest_pix = np.unravel_index(np.argmax(image_data), image_data.shape)
             #brightest_pix = [brightest_pix[0],brightest_pix[1]+300] # To adjust for trim to 300:350
             #ref_RA = 275.1101511 
             #ref_DEC = 7.169812515

             # For 221000-271000 in July, need a different reference pixel
             image_data = hdu_temp.data[0:95,20:80] # Clipping the top to avoid the horizontal dead pixel
             brightest_pix = np.unravel_index(np.argmax(image_data), image_data.shape)
             brightest_pix = [brightest_pix[0],brightest_pix[1]+20] # To adjust for trim to 20:80
             ref_RA = 275.13258 
             ref_DEC = 7.140235475

             print("\nThe pixel which will be used as reference is at "+str((brightest_pix[1],brightest_pix[0])))
             print("Which corresponds to (RA,DEC) = ", str((ref_RA,ref_DEC)),"\n")
             
             # Indicates that X=RA, Y=DEC, and TAN (gnomonic) is the projection to use 
             hdu_temp.header.append(('CTYPE1','RA---TAN'))
             hdu_temp.header.append(('CTYPE2','DEC--TAN'))

             # The change in RA, DEC in X and Y when projecting
             # 0.466''/px == 0.000129 degrees/px
             hdu_temp.header.append(('CDELT1', -0.000129))
             hdu_temp.header.append(('CDELT2', 0.000129))

             # The parameters of the rotation matrix
             # The angle accounts for the rotation of the frame 
             angle_rad = angle * np.pi/180.0
             hdu_temp.header.append(('PC1_1',np.cos(angle_rad)))
             hdu_temp.header.append(('PC1_2',-np.sin(angle_rad)))
             hdu_temp.header.append(('PC2_1',np.sin(angle_rad)))
             hdu_temp.header.append(('PC2_2',np.cos(angle_rad)))

             #hdu_temp.header.append(('PC1_1',np.cos(angle_rad)+np.sin(angle_rad)))
             #hdu_temp.header.append(('PC1_2',np.cos(angle_rad)+np.sin(angle_rad)))
             #hdu_temp.header.append(('PC2_1',-np.sin(angle_rad)+np.cos(angle_rad)))
             #hdu_temp.header.append(('PC2_2',-np.sin(angle_rad)+np.cos(angle_rad)))

             # OMM's values for 180709 data
             #hdu_temp.header.append(('PC1_1',-0.0000943745257614))
             #hdu_temp.header.append(('PC1_2',0.0000851308715725))
             #hdu_temp.header.append(('PC2_1',0.0000855725732586))
             #hdu_temp.header.append(('PC2_2',0.0000945508646851))

             # Pixel values of the reference coord
             hdu_temp.header.append(('CRPIX1', brightest_pix[1])) 
             hdu_temp.header.append(('CRPIX2', brightest_pix[0])) 
             
             # WCS values of the reference coord 
             hdu_temp.header.append(('CRVAL1', ref_RA)) 
             hdu_temp.header.append(('CRVAL2', ref_DEC))

             hdu_temp.writeto(self.loc[0]+'/'+self.name+'/'+f,'warn',overwrite=True) 

     def photometry(self,thresh_factor,RA_bounds,DEC_bounds):     # Added the threshold factor
        """
        Input: a threshold factor to be used in selecting what amount of background to ignore during 
        image segmentation and 2 arrays denoting the RA and DEC boundaries (in degrees) of the source 
        we wish to detect. 
        Example input: reduced_dataset.photometry(1.5, [275.1, 276.2], [7.10, 7.18])
        Produces a segmented image and a .csv with source properties (see apeturephotometry.py for more)
        """
        import aperturephotometry
        files = os.listdir(self.loc[0]+'/'+self.name)
        for f in files: 
             hdu_temp = fits.open(self.loc[0]+'/'+self.name+'/'+f,mode='update')[0]
             hdu_temp.header['EPOCH'] = float(hdu_temp.header['EPOCH']) # EPOCH needs to be a float   
             hdu_temp.writeto(self.loc[0]+'/'+self.name+'/'+f,'warn',overwrite=True) 
             hdu = fits.open(self.loc[0]+'/'+self.name+'/'+f)    
             aperturephotometry.photometry(hdu[0].header,hdu[0].data,f.replace('.fits', ''),thresh_factor,RA_bounds, DEC_bounds)

########################

