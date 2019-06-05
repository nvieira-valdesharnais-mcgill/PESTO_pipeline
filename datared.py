"""
@authors: Andy Ramirez-Cote, Valerie Desharnais & Nicholas Vieira
@datared.py

Stack domeflats and apply this stack to correct individual object images. Then, 
stack the corrected object images.
"""


import sys
import os

#Run iraf in terminal-only mode (no graphics)
from stsci.tools import capable
capable.OF_GRAPHICS = False

from pyraf import iraf 
from iraf import noao, imred, ccdred 
from iraf import immatch
from iraf import zerocombine, ccdproc, imcombine

# if desired, can output a file which is named something other than 
# "results.txt" by passing it as a cmd line arg
results_file = str(sys.argv[1]) 

for filter_type in ['r','g','i','z']:
    name = "object_list_"+str(filter_type)+".txt"
    subname = str(filter_type)+"_list.txt"
    if ('object_list_'+filter_type+'.txt' in os.listdir(os.getcwd())) and (
            os.stat(name).st_size!=0) and (os.stat(subname).st_size!=0):

        # bias corrections performed automatically at OMM, so this correction 
        # is unnecessary, but it is left here in case it is useful 
        #ccdproc('@object_list_'+filter_type+'.txt',
        #             zero='bias.fits',
        #             zerocor = iraf.yes,
        #             darkcor=iraf.no,
        #             flatcor=iraf.no,
        #             fixpix=iraf.no,
        #             overscan=iraf.no,     
        #             trim=iraf.no)

        #copying the stack number, and the average and stdevs of exposure (ms) 
        # and timestamp (s) to results.txt
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
            time.append(float(date[11:-13])*3600.0+
                        float(date[14:-10])*60.0+
                        float(date[17:-7])+
                        float(date[20:])/1000000.0) # convert to seconds
        exp = np.asarray(exp)
        time = np.asarray(time)
        estd = np.std(exp)
        if estd < 1e-5:
            estd = 0.0
        tf = open('/data/irulan/omm_transients/'+results_file,'a')

        # write the stack size, exposure time, error on exposure, timestamp, 
        # and error on timestamp
        tf.write(str(stack) + "\t" + str(np.mean(exp)) + "\t" +
                 str(estd) + "\t" + str(np.mean(time)) + "\t" +
                 str(np.std(time)) + "\t")
        tf.close()

        # apply the stacked flat correction to each object image, 
        # then stack the corrected objects
        ccdred.flatcombine('@'+filter_type+'_list.txt',
                           output='Flat_'+filter_type+'.fits',
                           combine='median',
                           reject='crreject')
        # zerocor=iraf.yes if bias.fits is supplied
        # include zero="bias.fits" arg if bias.fits is supplied 
        ccdproc('@object_list_'+filter_type+'.txt',
                     flat='Flat_'+filter_type+'.fits',
                     zerocor = iraf.no, 
                     darkcor=iraf.no,   
                     flatcor=iraf.yes, 
                     fixpix=iraf.no,
                     overscan=iraf.no,
                     trim=iraf.no)
        imcombine('@object_list_'+filter_type+'.txt', 
                  output='object_'+filter_type+'_reduced.fits')