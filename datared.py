import sys
import os

#Run iraf in terminal-only mode (no graphics)
from stsci.tools import capable
capable.OF_GRAPHICS = False

from pyraf import iraf 
from iraf import noao, imred, ccdred 
from iraf import immatch
from iraf import zerocombine, ccdproc, imcombine

#Move to raw_data directory
#print("\n\nCurrent working directory for data reduction: "+os.getcwd())

#Performs the bias correction and the flat correction for each band

#bias corrections are performed automatically at OMM, so there are no biases to stack
#if os.stat("bias_list.txt").st_size!=0: 
#    folder = os.getcwd()
#    os.chmod(folder+"/bias_list.txt",0o777)
#    zerocombine('@bias_list.txt', output='bias.fits', combine='median')

for filter_type in ['r','g','i','z']:
    name = "object_list_"+str(filter_type)+".txt"
    subname = str(filter_type)+"_list.txt"
    if 'object_list_'+filter_type+'.txt' in os.listdir(os.getcwd()) and os.stat(name).st_size!=0 and os.stat(subname).st_size!=0:

        #bias corrections are performed automatically at OMM, so this correction is unnecessary
        #ccdproc('@object_list_'+filter_type+'.txt',
        #             zero='bias.fits',
        #             zerocor = iraf.yes,
        #             darkcor=iraf.no,
        #             flatcor=iraf.no,
        #             fixpix=iraf.no,
        #             overscan=iraf.no,     
        #             trim=iraf.no)

        #copying the stack number, and the average and stdevs of exposure (ms) and timestamp (s) to results.txt
        f = open('object_list_'+filter_type+'.txt','r')
        contents = f.readlines()
        f.close()
        stack = len(contents)
        import numpy as np
        from astropy.io import fits
        exp = []
        time = []
        for line in contents:
            hdu = fits.open(line[:-1])[0]
            exp.append(hdu.header['EXPOSURE'])
            date = hdu.header['DATE']
            time.append(float(date[11:-13])*3600.0+float(date[14:-10])*60.0+float(date[17:-7])+float(date[20:])/1000000.0)
        exp = np.asarray(exp)
        time = np.asarray(time)
        estd = np.std(exp)
        if estd < 1e-5:
            estd = 0.0
        tf = open('/data/irulan/omm_transients/results.txt','a')

        # TO WRITE STACK #, EXPOSURE, ERROR ON EXPSOSURE, TIMESTAMP, AND ERROR ON TIMESTAMP TO results.txt
        tf.write(str(stack) + "\t" + str(np.mean(exp)) + "\t" + str(estd) + "\t" + str(np.mean(time)) + "\t" + str(np.std(time)) + "\t")

        # TO WRITE ONLY THE TIMESTAMP AND ERROR ON TIMESTAMP TO results.txt:
        #tf.write(str(np.mean(time)) + "\t" + str(np.std(time)))

        #tf.write("NA"+"\t"+"NA"+"\t"+"NA"+"\t"+str(np.mean(time))+"\t"+str(np.std(time))+"\t")
        tf.close()

        #performing the stacked flat correction to each object image, then stacking the corrected objects
        ccdred.flatcombine('@'+filter_type+'_list.txt',
                           output='Flat_'+filter_type+'.fits',
                           combine='median',
                           reject='crreject')
        ccdproc('@object_list_'+filter_type+'.txt',
                     flat='Flat_'+filter_type+'.fits',
                     zerocor = iraf.no,
                     darkcor=iraf.no,
                     flatcor=iraf.yes,
                     fixpix=iraf.no,
                     overscan=iraf.no,
                     trim=iraf.no)
        #if bias correction, change to zero='bias.fits',zerocor=iraf.yes for above ccdproc()
        imcombine('@object_list_'+filter_type+'.txt',
                       output='object_'+filter_type+'_reduced.fits')

