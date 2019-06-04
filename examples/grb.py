"""
"""

import PESTO_lib 
import subprocess

whomstdve = 'x'

while not(whomstdve in ['N','n','V','v']):
    whomstdve = input("\nAre you Nick [N/n] or Val [V/v]? \n> ")

if whomstdve in ["V","v"]:
    name = "reduced_val"
elif whomstdve in ["N","n"]:
    name = "mini_reduced_data" 

location = ["/exports/scratch/GRB180618A/grb_pesto_objs","/exports/scratch/GRB180618A/grb_pesto_domeflats"] 
target = "/exports/scratch/GRB180618A"
type = ["object","calibration"]

# Here, select the RA and DEC bounds, depending on whether we're looking at MAXI or an alt source
# Uncomment the arrays which you want to use
# THIS HAS TO BE UPDATED FOR THE GRB 

# FOR MAXI
RA = [275.085, 275.097]
DEC = [7.1849, 7.1860]

# Alternate source #1  (for MAXI, needs to be updated)
#RA = [275.128, 275.138]
#DEC = [7.1398, 7.1410]

data = PESTO_lib.raw_PESTO_data(location,type,name,target) 
data.flush()
data.produce_lists(NA=True) 
data.make_working_directory() 
data.pyraf_reduction()
    
reduced_data = data.extract_reduced_images("/exports/scratch/GRB180618A/"+name, "minired") 
reduced_data.WCS_preparation(42.0934) # THIS ANGLE NEEDS TO BE UPDATED
reduced_data.photometry(2.0, RA, DEC) 
