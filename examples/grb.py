"""
@authors: Nicholas Vieira & Valérie Desharnais
@grb.py

A script to correct and stack several object images for the source GRB180618A. 
"""

import PESTO_lib 
import subprocess


location = ["/exports/scratch/GRB180618A/grb_pesto_objs","/exports/scratch/GRB180618A/grb_pesto_domeflats"] 
target = "/exports/scratch/GRB180618A"
type = ["object","calibration"]

# Here, select the RA and DEC bounds, depending on whether we're looking at MAXI or an alt source
# Uncomment the arrays which you want to use
# THIS HAS TO BE UPDATED FOR THE GRB 

# FOR GRB
RA = [169.930, 169.940]
DEC = [73.842, 73.852]

data = PESTO_lib.raw_PESTO_data(location,type,name,target) 
data.flush()
data.produce_lists(NA=True) 
data.make_working_directory() 
data.pyraf_reduction()
    
reduced_data = data.extract_reduced_images("/exports/scratch/GRB180618A/"+name, "minired")

# as we do not have full frame images, the clunky manual technique for WCS calibration must be used
# this requirs that we know the angle of rotation (relative to cardinal points) of the telescope
# the following angle is not correct
reduced_data.WCS_preparation(42.0934) 
reduced_data.photometry(2.0, RA, DEC) 
