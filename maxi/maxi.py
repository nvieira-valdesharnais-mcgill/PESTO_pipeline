"""
This script was used specifically to pipe data on the rapid binary black hole system MAXI J1820+070 in 2018 and 2019. 
It was used on the irulan server of the physics department of McGill University. Itallowed the users (Nick, or Val) to 
obtain lightcurves of this object using either stacks of 100 or 1000 images for 2 epochs: 9 July 2018 or 28 September 2018. 
The script runs for either 100 stacks of 100 images or 10 stacks of 1000 images. 

A similar script which can be run indefinitely (until forced to stop), maxi_overnight.py, was also used. 
"""

import PESTO_lib 
import subprocess

whichdate = 'x'
whomstdve = 'x'
filestdve = 'x'

while not(whichdate in ['180709','180928']):
    whichdate = input("What date do you want? Enter 180709 or 180928 \n> ") 
while not(whomstdve in ['N','n','V','v']):
    whomstdve = input("\nAre you Nick [N/n] or Val [V/v]? \n> ")
while (filestdve == 'x'):
    filestdve = input("\nWhat is the integer starting point? \nFor stacks of 1000, 180709: 97 to 527; 180928: 31 to 462.\nFor stacks of 100 (currently in use) , 180709: 970 to 5284; 180928: 310 to 4634.\n> ")

filestdve_num = int(filestdve)

if whomstdve in ["V","v"]:
    name = "mini_reduced_val"
elif whomstdve in ["N","n"]:
    name = "mini_reduced_data" 

###### FORMAT
# JULY: 096931 to 528554
# SEPT: 030054 to 463776

if whichdate == '180709':
    location = ["/exports/scratch/MAXIJ1830/180709_works","/exports/scratch/MAXIJ1830/180709_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1830/180709/Target/MAXIJ1820+070/180709_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1830/180709_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1830/180709_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1830/180709_works/*"
elif whichdate == '180928':
    location = ["/exports/scratch/MAXIJ1830/180928_works","/exports/scratch/MAXIJ1830/180928_calibs"] 
    bash_copy_start = "cp /exports/scratch/MAXIJ1830/180928/MAXI1820+070-180928-OMM/Target/180928_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1830/180928_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1830/180928_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1830/180928_works/*"

og_filestdve_num = filestdve_num
target = "/exports/scratch/MAXIJ1830/"
type = ["object","calibration"]

# Here, select the RA and DEC bounds, depending on whether we're looking at MAXI or an alt source
# Uncomment the arrays which you want to use

# FOR MAXI
#RA = [275.085, 275.097]
#DEC = [7.1849, 7.1860]

# Source 1 (ref pix) in Nick's finderchart 
# Very bright, sometimes saturates the CCD
# Worked in July, not September
#RA = [275.128, 275.138]
#DEC = [7.1398, 7.1410]

# Source 2 in Nick's finderchart                                                                                       
# Worked in July, not in September 
#RA = [275.108, 275.115]
#DEC = [7.1655, 7.1775] 

# Source 3 in Nick's finderchart, used with threshfactor 3, for July and Sept
#RA = [275.09, 275.09]
#DEC = [7.180, 7.184]
                                                                                                                         
# Source 5 in Nick's finderchart                                                                                       
# Not tested on July, used for September
#RA = [275.1195, 275.1205]
#DEC = [7.1540, 7.1552]                                                                                                    
                                                                                                                           
# Source 6 in Nick's finderchart                                                                                       
# Not tested on July, used for September
RA = [275.072, 275.074]
DEC = [7.199, 7.201]       

while(filestdve_num < og_filestdve_num+100): #do ten iterations #added 0
    if filestdve_num < 1000: #added 0
        bash_copy = bash_copy_start+"0"+str(filestdve_num)+bash_copy_end
    else:
        bash_copy = bash_copy_start+str(filestdve_num)+bash_copy_end 
    subprocess.run(bash_empty,shell=True)
    subprocess.run(bash_copy,shell=True)
    subprocess.run(bash_gunzip,shell=True)
    filestdve_num = filestdve_num+1

    data = PESTO_lib.raw_PESTO_data(location,type,name,target) 
    data.flush()
    data.produce_lists(NA=True) 
    data.make_working_directory() 
    data.pyraf_reduction()

    reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1830/"+name, "minired") 
    if whichdate == '180709':
        reduced_data.WCS_preparation(42.0934)
        reduced_data.photometry(2.0, RA, DEC) #changed from 7.0 for 1k stacks, 3.0 for 100 stacks 
    elif whichdate == '180928':
        reduced_data.WCS_preparation(40.3927)
        reduced_data.photometry(2.0, RA, DEC) #changed from 7.0 for 1k, 3.0 for 100


