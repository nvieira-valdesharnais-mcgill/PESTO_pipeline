"""
@authors: Valérie Desharnais & Nicholas Vieira
@PESTO_lib.py

A library of functions to reduce, stack, and calibrate images from the Planètes 
Extra-Solaires en Transit et Occultations (PESTO) optical telescope located at 
l’Observatoire du Mont-Mégantic (OMM)
"""
import os
from subprocess import run
from astropy.io import fits
import numpy as np


class PESTO_data:
     """
     Input: 
     locations_list: a list of strings, one for each directory with PESTO data
     image_type: a list whose elements are the strings 'calibration' or 
                 'object', describing the content of the files in each location 
                 in locations_list
     name: a name to use in the creation of new files
     
     Example input (1):
     (['path/to/object','path/to/flat'],['object','calibration'], foo)
     Example input (2): 
         (['homes/huanca/aramirez/SkyFlat','homes/huanca/aramirez/GRB180618A'],
         ['calibration','object'], bar)
         
     Output: PESTO_data object
     """
     def __init__(self, locations_list, image_type_list, name):
          self.loc = locations_list
          self.imtype = dict([(self.loc[i],image_type_list[i]) for 
                              i in range(0,len(self.loc))])
          self.name = name
          self.r_obj = []
          self.r_cal = []
          self.g_obj = []
          self.g_cal = []
          self.i_obj = []
          self.i_cal = []
          self.z_obj = []
          self.z_cal = []
          self.bias = []
          self.hdr_ind = {} # dictionary of locations and indices to read when
                            # looking at headers in these locations 
 
     def untar(self):
          """
          Input: None
          Output: None      
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
         Output: None
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
         Flushes the working directory of unneeded files and directories.
         """
         temp = os.getcwd()
    
         os.chdir(self.tgt+'/'+self.name)
         run('rm -rf ./*/', shell=True) # directories  
         run('rm -rf *.py', shell=True) # scripts
         run('rm -rf *.fits', shell=True) # fits files 
         run('rm -rf *.txt', shell=True) # text files (lists) 
         run('rm -rf *.par', shell=True) # related to iraf
         run('rm -rf logfile', shell=True) # related to iraf

         os.chdir(temp)
         print("\nFlush!")

###############################################################################

class raw_PESTO_data(PESTO_data):
     
    def __init__(self,locations_list,image_type_list,name, target_directory):
        super(raw_PESTO_data, self).__init__(locations_list,image_type_list,
             name)
        self.tgt = target_directory
        self.reduced = False # updated later 
        self.list_made = False
        self.workdir_present = False

        # allow new data to be read regardless of header indexing.
        # assumption: only need to look at one .fits per location, and if 
        # extended, only need to look at the first extension
        for l in self.loc: # for all files in data directories 
            files = os.listdir(l)
            for f in files:
                if '.fits' in f: # look until a .fits is found
                    tempfile = f
                hdr_temp = fits.open(l+'/'+tempfile)
                
                # unextended file
                # if 0th header contains FILTRE keyword, use index 0 for all
                # header reading in this directory
                if ('FILTRE' in hdr_temp[0].header): 
                    self.hdr_ind[l] = 0
                # extended file
                elif ('FILTRE' in hdr_temp[1].header):  
                    self.hdr_ind[l] = 1  

###############################################################################
           
    def extended_header_cleanup(self):
        """
        Input: None
        Output: None
        In certain extended .fits files, the SIMPLE, EXTEND, and/or NEXTEND 
        headers are misplaced. This function repairs these files, PERMANENTLY 
        modifying them. For a given data object, only needs to be run once.
        WARNING: This PERMANENTLY modifies data in a given directory.
        """
        for l in self.loc: # for all files in data directories 
            files = os.listdir(l)
            for f in files: 
                if '.fits' in f:
                    if self.hdr_ind[l] != 0: # if looking at extended files 
                        # get no. of headers and begin updating headers 
                       n_extend = fits.open(l+'/'+f)[0].header['NEXTEND']
                       hdr_temp = fits.open(l+'/'+f,mode='update')
                       for n in range(1, n_extend+1):
                           hdu = hdr_temp[n]
                           if 'SIMPLE' in hdu.header:
                               temp_simple = hdu.header['SIMPLE']
                               # shift headers up by 1, append header to end
                               hdu.header.remove('SIMPLE')  
                               hdu.header.append(('SIMPLE',temp_simple)) 
                           if 'EXTEND' in hdu.header:
                               # shift headers up by 1, append header to end
                               hdu.header.remove('EXTEND')   
                               hdu.header.append(('EXTEND','T')) 
                           if 'NEXTEND' in hdu.header:
                               # shift headers up by 1, append header to end
                               hdu.header.remove('NEXTEND') 
                               hdu.header.append(('NEXTEND',n_extend)) 
                hdr_temp.close()

    def resize(self, sf, df, mode, sizes):
        """
        Input: source filepath, destination filepath, 'trim' or 'extend' mode, 
        list of pixel dimensions [left, right, top, bottom]
        Output: None
        Resizes the image to the appropriate dimensions by either trimming or 
        extending by NaNs the pixel dimensions passed. 
        WARNING: This method may PERMANENTLY ALTER/OVERWRITE DATA in both 
        source and destination folders. 
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
                        if (not all(s>=0 for s in sizes)) or (
                                sizes[0]+sizes[1])>x or (sizes[2]+sizes[3])>y:
                            print("Dimensions to trim are out of range for: "+sf)
                            return
                        new_image_data = image_data[(sizes[2]):(y-sizes[3]), (
                                sizes[0]):(x-sizes[1])]
                        image_header["NAXIS1"]=x-(sizes[0]+sizes[1])
                        image_header["NAXIS2"]=y-(sizes[2]+sizes[3])
                        fits.append(df,new_image_data,image_header)
                    elif mode == "extend":
                        if not all(s>=0 for s in sizes):
                            print("Dimensions to extend by are out of range for: "+sf)
                            return
                        new_image_data = np.hstack((np.full((y,sizes[0]),np.nan), 
                                                    image_data, 
                                                    np.full((y,sizes[1]),np.nan)))
                        new_image_data = np.vstack(
                                (np.full((sizes[2],(x+sizes[0]+sizes[1])),np.nan), 
                                 new_image_data, 
                                 np.full((sizes[3],(x+sizes[0]+sizes[1])),np.nan)))
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
        Input: source folder path, destination folder name, 'trim' or 'extend' 
        mode, list of pixel dimensions [left, right, top, bottom]
        Output: None
        Resizes all fit files in a folder to the appropriate dimensions by 
        either trimming or extending by NaNs the pixel dimensions passed and 
        deposits in subfolder.
        WARNING: This method may PERMANENTLY ALTER/OVERWRITE DATA in both 
        source and destination folders. 
        """
        files = os.listdir(source)
        if not os.path.exists(source+'/'+destination_name):
            os.mkdir(source+'/'+destination_name)
        for f in files:
            self.resize(source+'/'+f,source+'/'+destination_name+
                        '/'+mode+'_'+f,mode,sizes)    
 
###############################################################################

    def produce_lists(self):
        """
        Input: None
        Output: None
        Produces lists of .fits files associated with the raw_PESTO_data 
        object. e.g all images of the object taken in the r band are added to 
        self.r_obj, all calibration images in the r band are added to 
        self.r_cal, etc.
        """
        if not self.list_made:
            for l in self.loc:			
                files = os.listdir(l)
                run_type = self.imtype[l]   
                for f in files:
                    if '.fits' in f:
                        hdr_temp = fits.open(l+'/'+f, mode = 'update')
                        hdu = hdr_temp[self.hdr_ind[l]] 
                        hdu.header['NAXIS'] = 3

                        # adds the type of frame (flat, bias, object, etc...) 
                        # to the fits header to help iraf's data reduction
                        ## CALIBRATION files (biases, flats)
                        if run_type == 'calibration':
                            # check 2 conditions to see if it's a flat
                            condition1 = ('OBJECT' in hdu.header) and (
                                    'flat' in hdu.header['object'])
                            condition2 = ('IMAGETYP' in hdu.header) and (
                                    'flat' in hdu.header['imagetyp'])
                            isflat = condition1 or condition2
                            # produce a list of biases
                            if ('IMAGETYP' in hdu.header) and (
                                    'zero' in hdu.header['imagetyp']):
                                self.bias.append(f.replace('.gz',''))
                            # roduce lists of flats in each band 
                            elif isflat:
                                if 'r' in hdu.header['filtre']:
                                    self.r_cal.append(f.replace('.gz',''))
                                elif 'g' in hdu.header['filtre']:
                                    self.g_cal.append(f.replace('.gz',''))
                                elif 'i' in hdu.header['filtre']:
                                    self.i_cal.append(f.replace('.gz',''))
                                elif 'z' in hdu.header['filtre']:
                                    self.z_cal.append(f.replace('.gz',''))

                        ## OBJECT files 
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

                self.list_made=True # update this bool
    
    def make_working_directory(self):
        """
        Input: None
        Output: None
        Creates a directory containing all .fits files associated to the 
        object and one .txt file for each list bound to the raw_PESTO_data 
        object. The created directory contains everything PyRAF needs for the 
        data reduction.   
        """
        if not self.workdir_present:
             l = self.tgt+'/'+self.name

             ## IMPORTANT: copy over everything needed for iraf 
             # the files must be in a directory in the user's home, and the 
             # directory must be named iraf 
             # most important: login.cl, which sets the conditions for iraf 
             # see the github for an example file 
             run(['mkdir','-p', l])
             run('cp -f ~/iraf/login.cl '+l, shell=True)
             run('mkdir -p '+l+'/uparm', shell=True)
             run('cp -f ~/iraf/uparm/* '+l+'/uparm', shell=True) 
             
             # writes the list for each filter to the working directory
             np.savetxt(l+'/r_list.txt',self.r_cal, fmt = "%s")
             np.savetxt(l+'/g_list.txt', self.g_cal, fmt = "%s")
             np.savetxt(l+'/i_list.txt', self.i_cal, fmt = "%s")
             np.savetxt(l+'/z_list.txt', self.z_cal, fmt = "%s")
             np.savetxt(l+'/bias_list.txt', self.bias, fmt = "%s")
             np.savetxt(l+'/object_list_r.txt', self.r_obj, fmt = "%s")
             np.savetxt(l+'/object_list_g.txt', self.g_obj, fmt = "%s")
             np.savetxt(l+'/object_list_i.txt', self.i_obj, fmt = "%s")
             np.savetxt(l+'/object_list_z.txt', self.z_obj, fmt = "%s")
            
             self.workdir_present = True

             # copy all files in loc directories to the working directory
             # give full permissions to everything 
             for l in self.loc:
                 run(['rsync','-r','-p',l+'/',self.tgt+'/'+self.name]) 
             run('chmod 777 '+self.tgt+'/'+self.name+'/*', shell=True)
 
    def delete_working_directory(self):
         """
         Input: None
         Output: None
         Removes the working directory if it is present
         """
         if self.workdir_present:
              run(['rm', '-r', self.tgt+'/'+self.name])
              self.workdir_present=False

    def pyraf_reduction(self, results_file="results.txt"):
        """
        Input: The name of the results file to which we save the data 
        (optional; "results.txt" by default)
        Output: None
        """
        if self.list_made==False:
            return 'Please make sure lists for each band were produced'
        if self.workdir_present==False:
            return 'Please make sure a working directory exists'
        temp = os.getcwd()
        run(['rsync','datared.py',self.tgt+'/'+self.name])
        os.chdir(self.tgt+'/'+self.name)
        run("bash -c 'source activate iraf27 && python2 datared.py "+
            results_file+" && source deactivate'", shell=True)
        os.chdir(temp)

    def extract_reduced_images(self, reduced_data_folder_loc, 
                               reduced_data_folder_name):
        """
        Input: the parent directory of the reduced data directory, and the name 
        of the directory itself 
        Output: a reduced_PESTO_data object 
        Creates the necessary directories to begin WCS/photometric calibration
        of the stacked image.
        """
        # create necessary directory, copy over stack's .fits file 
        run(['mkdir', reduced_data_folder_loc+'/'+reduced_data_folder_name])
        run('cp -r '+self.tgt+'/'+self.name+'/*reduced.fits '+
            reduced_data_folder_loc+'/'+reduced_data_folder_name, shell=True) 
        return reduced_PESTO_data([reduced_data_folder_loc], ['object'], 
                                  reduced_data_folder_name)

###############################################################################

class reduced_PESTO_data(PESTO_data):
     def __init__(self, locations_list, image_type_list, name):
        super(reduced_PESTO_data, self).__init__(locations_list, 
             image_type_list, name)
        self.reduced = True

     def astrometry(self, roi_x, roi_y, ra_est=0, dec_est=0, rad_est=0, 
                    min_scale=0, max_scale=0, units_scale='x'):
       """
       Input: the bounds of the region of interest in x pixels and in y pixels 
       as arrays. **Optional: estimates for the RA, DEC, radius (in degrees) 
       and minimum/maximum scales of the image(s). 
       Output: None
       Uses the tool astrometry.net to solve the field of the reduced data 
       object(s) to prepare for source detection. Default units for scale is 
       arcseconds per pixel. RA, Dec, and radius guesses are ignored by 
       astrometry unless all 3 quantities are estimated. Same is true for 
       min_scale and max_scale.
       """
       temp = os.getcwd()
       os.chdir(self.loc[0]+'/'+self.name) # move to reduced data directory
       files = os.listdir()
      
       # reference pixel is at the center of the region of interest
       center_x = roi_x[0] + (roi_x[1] - roi_x[0])/2.0 
       center_y = roi_y[0] + (roi_y[1] - roi_y[0])/2.0 

       for f in files:
            # new filename includes _wcs to distinguish it
            newname = f.replace(".fits","_wcs.fits") 
            
            # give options to solve-field command: overwrite existing WCS 
            # solutions, produce no .png plots, tell solver that the input is a
            # fits image, and produce a new fits file with name newname
            options = " --overwrite --no-plots --fits-image"
            options += " --new-fits "+newname+" --cancel "+newname 
            options += " --crpix-x "+str(center_x)+" --crpix-y "+str(center_y)

            # if ra, dec, radius estimates are given
            if ((ra_est != 0) and (dec_est != 0) and (rad_est != 0)):
                 options += " --ra "+str(ra_est)+" --dec "+str(dec_est)
                 options += " --rad "+str(rad_est)

            # if scales estimates are given
            # does not check if units given are valid, but astrometry notices 
            if ((max_scale !=0) and (min_scale != 0) and (units_scale != 'x')):
                 options += " --scale-low "+str(min_scale)
                 options += " --scale-high "+str(max_scale)
                 options += " --scale-units "+units_scale

            # more options (experimenting)
            #options += " -v" # verbose
            #options += " --no-verify"  # speed up CPU time by not looking at WCS headers
            options += " --downsample 1" # decrease downsample -> increase star count (works, keep this)
            options += " --pixel-error 0.1" # decrease pixel error (default 1) -> increase star count 
            #options += " --nsigma 6" # decrease sigma required for source detection (default 8) -> increase source count
            #options += " --odds-to-solve 1e5" # decrease odds required to come to a soln (default 1e9) -> increase odds

            run("solve-field "+str(options)+" "+f, shell=True)
            run("find . -type f -not -name '*wcs*' -print0 | xargs -0 rm --",
                shell=True) # remove all files not in format *wcs*

       os.chdir(temp) # return to original directory

     def WCS_merge(self, wcs_location, delta_x=0, delta_y=0):
        """
        Input: The location of the WCS solution to be merged with the reduced 
        data object, and values for the difference in x=0 and y=0 pixel 
        positions between the WCS solution and the object images. 
        Output: None
        Takes an existing WCS solution and 
        """ 

        hdu_wcs = fits.open(wcs_location, mode='readonly')[0]

        temp = os.getcwd()
        os.chdir(self.loc[0]+'/'+self.name) # move to reduced data directory
        files = os.listdir()
        for f in files:
             hdu_temp = fits.open(
                     self.loc[0]+'/'+self.name+'/'+f,mode='update')[0]

             # in practice, delta_y encoded this way
             # delta_x is unfortunately not recorded anywhere and can only be 
             # observed by examining the fits files directly (via e.g. DS9)
             if delta_y == 0: # if no argument given to override this
                  delta_y = (-1.0)*hdu_temp.header['ROI_Y_2'] 
             
             hdu_temp.header.append(('CTYPE1','RA---TAN-SIP')) # projection type
             hdu_temp.header.append(('CTYPE2','DEC--TAN-SIP'))
             # x ref from WCS soln minus delta_x
             hdu_temp.header.append(('CRPIX1', 
                                     (float(hdu_wcs.header['CRPIX1'])+delta_x)))  
             hdu_temp.header.append(('CRPIX2', 
                                     (float(hdu_wcs.header['CRPIX2'])+delta_y)))
             # y ref from WCS soln minus delta_y
             hdu_temp.header.append(('CRVAL1', 
                                     float(hdu_wcs.header['CRVAL1']))) # in WCS 
             hdu_temp.header.append(('CRVAL2', 
                                     float(hdu_wcs.header['CRVAL2']))) # in WCS

             # the change in RA, Dec in X and Y when projecting
             # 0.46 arcsec per pixel == 0.000128 degrees/px
             hdu_temp.header.append(('CDELT1', -0.000128))
             hdu_temp.header.append(('CDELT2', 0.000128))

             # parameters of the rotation matrix
             # the angle accounts for the rotation of the frame 
             hdu_temp.header.append(('CD1_1', float(hdu_wcs.header['CD1_1'])))
             hdu_temp.header.append(('CD1_2', float(hdu_wcs.header['CD1_2'])))
             hdu_temp.header.append(('CD2_1', float(hdu_wcs.header['CD2_1'])))
             hdu_temp.header.append(('CD2_2', float(hdu_wcs.header['CD2_2'])))
 
             hdu_temp.writeto(self.loc[0]+'/'+self.name+'/'+f,'warn',
                              overwrite=True) # merge the WCS solution

        os.chdir(temp)

     def WCS_preparation(self, angle=0):
        """
        Input: The angle of the frame relative to a frame where x=RA, y=DEC, 
        counterclockwise, in degrees (optional; default is 0)
        Output: None
        Modifies the reduced fits file to provide it all the information the 
        WCS conversion needs.
        
        ** This method is necessary when fullframe data is not present, as was 
        sometimes the case for PESTO data. It requires manual modification of
        the image clipping arguments and those that immediately follow 
        (flagged with ### *MOD*)
        """
        files = os.listdir(self.loc[0]+'/'+self.name)
        for f in files: 
             hdu_temp = fits.open(
                     self.loc[0]+'/'+self.name+'/'+f,mode='update')[0]

             ### *MOD* # these lines need to be modified if function is used
             image_data = hdu_temp.data[0:95,20:80] # clip top of image 
             brightest_pix = np.unravel_index(np.argmax(image_data), 
                                              image_data.shape)
             # adjust for trim to 20:80
             brightest_pix = [brightest_pix[0],brightest_pix[1]+20] 
             ref_RA = 275.13258 
             ref_DEC = 7.140235475
             # *MOD* ###

             print("\nThe pixel which will be used as reference is at "+
                   str((brightest_pix[1],brightest_pix[0])))
             print("Which corresponds to (RA,DEC) = ", str((ref_RA,ref_DEC)),
                   "\n")
             
             # indicates that X=RA, Y=DEC, and TAN (gnomonic) is the projection
             hdu_temp.header.append(('CTYPE1','RA---TAN'))
             hdu_temp.header.append(('CTYPE2','DEC--TAN'))

             # the change in RA, DEC in X and Y when projecting
             # 0.466''/px == 0.000129 degrees/px
             hdu_temp.header.append(('CDELT1', -0.000129))
             hdu_temp.header.append(('CDELT2', 0.000129))

             # the parameters of the rotation matrix
             # the angle accounts for the rotation of the frame 
             angle_rad = angle * np.pi/180.0
             hdu_temp.header.append(('PC1_1',np.cos(angle_rad)))
             hdu_temp.header.append(('PC1_2',-np.sin(angle_rad)))
             hdu_temp.header.append(('PC2_1',np.sin(angle_rad)))
             hdu_temp.header.append(('PC2_2',np.cos(angle_rad)))

             # Pixel values of the reference coord
             hdu_temp.header.append(('CRPIX1', brightest_pix[1])) 
             hdu_temp.header.append(('CRPIX2', brightest_pix[0])) 
             
             # WCS values of the reference coord 
             hdu_temp.header.append(('CRVAL1', ref_RA)) 
             hdu_temp.header.append(('CRVAL2', ref_DEC))

             hdu_temp.writeto(self.loc[0]+'/'+self.name+'/'+f,'warn',
                              overwrite=True) 

     def photometry(self, RA_bounds, DEC_bounds, thresh_factor=3.0,
                    results_file="results.txt"):
        """
        Input: a threshold factor to be used in selecting what level of 
        background to ignore during image segmentation, 2 arrays denoting the 
        RA and Dec boundaries (in degrees) of the source to detect, and the 
        results file to append to (optional; default is 'results.txt' 
        ***If a non-default name is used in pyraf_reduction(), the same 
        filename must be used here.
        Output: None
        
        e.g. reduced_dataset.photometry([275.1,276.2], [7.10,7.18], 3.5)
        Produces a segmented image and a .csv with properties for all the 
        sources and appends properties of the specified source to the results 
        file (if it is found)
        """
        import aperturephotometry
        files = os.listdir(self.loc[0]+'/'+self.name)
        for f in files: 
             hdu_temp = fits.open(
                     self.loc[0]+'/'+self.name+'/'+f,mode='update')[0]
             # adjust EPOCH header to be a float 
             hdu_temp.header['EPOCH'] = float(hdu_temp.header['EPOCH'])
             hdu_temp.writeto(
                     self.loc[0]+'/'+self.name+'/'+f,'warn',overwrite=True) 
             hdu = fits.open(self.loc[0]+'/'+self.name+'/'+f)    
             aperturephotometry.photometry(hdu[0].header,
                                           hdu[0].data,f.replace('.fits', ''),
                                           RA_bounds, DEC_bounds, thresh_factor,
                                           results_file)

###############################################################################

