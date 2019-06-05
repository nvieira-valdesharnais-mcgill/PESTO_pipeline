"""
@authors: Nicholas Vieira & Valérie Desharnais 
@wcs_maker.py

A script for running the PESTO pipeline on full-frame data. It creates a stack 
of all dome flats, subtracts this stack individually from each fullframe image, 
and then stacks the full frames. The resultant reduced data file contains a WCS 
solution which can then be applied to all of the object images, provided that 
the instrument does not significantly drift during observations. 

This script was designed specifically to work with the source MAXIJ1820+070, 
but serves as an example for writing your own scripts to automate the process. 
"""

import PESTO_lib

whichdate = 'x'

while not(whichdate in ['190312','190317','190318', '190326', '190404']):
    whichdate = input("What date do you want? Enter 190312, 190317, 190318, 190326 and 190404 \n> ") 

if (whichdate == '190312'):
    location = ["/exports/scratch/MAXIJ1820/190312_fullframes",
                "/exports/scratch/MAXIJ1820/190312_calibs"]
    roi_x = [0.0,1024.0]
    roi_y = [385.0, 490.0]
elif (whichdate == '190317'):
    location = ["/exports/scratch/MAXIJ1820/190317_fullframes",
                "/exports/scratch/MAXIJ1820/190317_calibs"]
    roi_x = [0.0,1024.0]
    roi_y = [460.0, 566.0] 
elif (whichdate == '190318'):
    location = ["/exports/scratch/MAXIJ1820/190318_fullframes",
                "/exports/scratch/MAXIJ1820/190318_calibs"]
    roi_x = [0.0,1024.0]
    roi_y = [420.0, 526.0]
elif (whichdate == '190326'):
    location = ["/exports/scratch/MAXIJ1820/190326_fullframes",
                "/exports/scratch/MAXIJ1820/190326_calibs"]
    roi_x = [45.0,1069.0] # have to shift 45 pix because PESTO moves 45 pix rightward after taking fullframes
    roi_y = [423.0, 523.0]
elif (whichdate == '190404'):
    location = ["/exports/scratch/MAXIJ1820/190404_fullframes",
                "/exports/scratch/MAXIJ1820/190404_calibs"]
    roi_x = [0.0,1024.0]
    roi_y = [415.0, 520.0]

name = "mini_reduced_data" # reduced data directory
target = "/exports/scratch/MAXIJ1820" # location of reduced data directory
imtype = ["object","calibration"] # types of files being used 
# object files: the full frame images
# calibration files: dome flats 

data = PESTO_lib.raw_PESTO_data(location,imtype,name,target)
data.flush()
data.produce_lists()
data.make_working_directory()
data.pyraf_reduction()

reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, 
                                           "reduc")

ra_guess = 275.09
dec_guess = 7.18
rad_guess = 0.1
min_scale_guess = 0.1
max_scale_guess = 1.0
scale_units = "app" # arcsec per pixel

reduced_data.astrometry(roi_x, roi_y, ra_guess, dec_guess, rad_guess, 
                        min_scale_guess, max_scale_guess, scale_units)
