"""
@authors: Nicholas Vieira & Valérie Desharnais
@maxi_overnight.py

This script was used specifically to pipe data on the transient binary black hole system MAXI J1820+070 in 2018 and 2019. 
It was used on the irulan server of the physics department of McGill University. It allowed the users to obtain 
lightcurves of this object using either stacks of 100 or 1000 images for several different epochs. 
The script runs for either 100 stacks (iters) of 100 images or 10 stacks (iters) of 1000 images. 

A similar script which can be run indefinitely (until forced to stop), maxi_overnight.py, was also used. 
"""

import PESTO_lib 
import subprocess

whichdate = 'x'
whomstdve = 'x'
filestdve = 'x'
stackstdve = 'x'

while not(whichdate in ['180709','180928','190312','190317','190318', '190326', '190404']):
    whichdate = input("What date do you want? Enter 180709, 180928, 190312, 190317, 190318, 190326 or 190404 \n> ") 
while not(whomstdve in ['N','n','V','v']):
    whomstdve = input("\nAre you Nick [N/n] or Val [V/v]? \n> ")
while not (stackstdve in ['100', '1000']):
    stackstdve = input("\nWhat is the stack size? Pick 100 or 1000 \n> ")
while (filestdve == 'x'):
    filestdve = input("\nWhat is the integer starting point? Pick using the table below:\n\
    +--------+------------+----------+\n\
    | Date   | Stack size | Range    |\n\
    +--------+------------+----------+\n\
    | 180709 | 100        | 1070-5284|\n\
    +        +------------+----------+\n\
    |        | 1000       | 107-527  |\n\
    +--------+------------+----------+\n\
    | 180928 | 100        | 310-4634 |\n\
    +        +------------+----------+\n\
    |        | 1000       | 31-462   |\n\
    +--------+------------+----------+\n\
    | 190312 | 100        | 83-2560  |\n\
    +        +------------+----------+\n\
    |        | 1000       | 9-255    |\n\
    +--------+------------+----------+\n\
    | 190317 | 100        | 428-5840 |\n\
    +        +------------+----------+\n\
    |        | 1000       | 43-584   |\n\
    +--------+------------+----------+\n\
    | 190318 | 100        | 99-5659  |\n\
    +        +------------+----------+\n\
    |        | 1000       | 9-565    |\n\
    +--------+------------+----------+\n\
    | 190326 | 100        | 220-4373 |\n\
    +        +------------+----------+\n\
    |        | 1000       | 22-437   |\n\
    +--------+------------+----------+\n\
    | 190404 | 100        | 83-6134  |\n\
    +        +------------+----------+\n\
    |        | 1000       | 8-613    |\n\
    +--------+------------+----------+\n\n> ")

filestdve_num = int(filestdve)

if whomstdve in ["V","v"]:
    name = "mini_reduced_val"
elif whomstdve in ["N","n"]:
    name = "mini_reduced_data" 

###### FORMAT
# JULY: 096931 to 528554
# SEPT: 030054 to 463776
# MARCH 12: 008241 to 334055
# MARCH 17: 042805 to 584457
# MARCH 18: 009878 to 565912
# MARCH 26: 022000 to 437399
# APRIL 04: 008235 to 613404

if whichdate == '180709':
    location = ["/exports/scratch/MAXIJ1820/180709_works","/exports/scratch/MAXIJ1820/180709_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/180709/Target/MAXIJ1820+070/180709_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1820/180709_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/180709_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/180709_works/*"
elif whichdate == '180928':
    location = ["/exports/scratch/MAXIJ1820/180928_works","/exports/scratch/MAXIJ1820/180928_calibs"] 
    bash_copy_start = "cp /exports/scratch/MAXIJ1820/180928/MAXI1820+070-180928-OMM/Target/180928_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1820/180928_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/180928_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/180928_works/*"
elif whichdate == '190312':
    location = ["/exports/scratch/MAXIJ1820/190312_works","/exports/scratch/MAXIJ1820/190312_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190312/MAXIJ1820+070/190312_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1820/190312_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190312_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190312_works/*"
elif whichdate == '190317':
    location = ["/exports/scratch/MAXIJ1820/190317_works","/exports/scratch/MAXIJ1820/190317_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190317/MAXIJ1820+070/190317_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1820/190317_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190317_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190317_works/*"
elif whichdate == '190318':
    location = ["/exports/scratch/MAXIJ1820/190318_works","/exports/scratch/MAXIJ1820/190318_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190318/MAXIJ1820+070/190318_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1820/190318_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190318_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190318_works/*"
elif whichdate == '190326':
    location = ["/exports/scratch/MAXIJ1820/190326_works","/exports/scratch/MAXIJ1820/190326_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190326/MAXIJ1820+070/190326_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1820/190326_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190326_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190326_works/*"
elif whichdate == '190404':
    location = ["/exports/scratch/MAXIJ1820/190404_works","/exports/scratch/MAXIJ1820/190404_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190404/MAXIJ1820+070/190404_0000"
    bash_copy_end = "* /exports/scratch/MAXIJ1820/190404_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190404_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190404_works/*"

og_filestdve_num = filestdve_num
target = "/exports/scratch/MAXIJ1820"
type = ["object","calibration"]

# Here, select the RA and DEC bounds, depending on whether we're looking at MAXI or an alt source
# Uncomment the arrays which you want to use

# FOR MAXI
RA = [275.087, 275.095]
DEC = [7.184, 7.187]

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
#RA = [275.072, 275.074]
#DEC = [7.199, 7.201]

if (int(stackstdve) == 100):   # if making stacks of 100
    filestdve_iters  = 100     # 100 iterations
    filestdve_lim1 = 100      
    filestdve_lim2 = 1000      
elif (int(stackstdve) == 1000): # if making stacks of 1000
    filestdve_iters  = 10      # 10 iterations
    filestdve_lim1 = 10
    filestdve_lim2 = 100 

while(filestdve_num < og_filestdve_num+filestdve_iters): 
    if filestdve_num < filestdve_lim1:  # check if file is less than ...10000.fits
        bash_copy = bash_copy_start+"00"+str(filestdve_num)+bash_copy_end
    elif filestdve_num < filestdve_lim2: # check if file is less than ..100000.fits
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

    if whichdate == '180709':
        data.pyraf_reduction("results_july.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        reduced_data.WCS_preparation(42.0934)     # old WCS technique
        reduced_data.photometry(2.0, RA, DEC, "results_july.txt")
  
    elif whichdate == '180928':
        data.pyraf_reduction("results_sept.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        reduced_data.WCS_preparation(40.3927)     # old WCS technique
        reduced_data.photometry(2.0, RA, DEC, "results_sept.txt")

    elif whichdate == '190312':
        data.pyraf_reduction("results_190312.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190312_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(2.0, RA, DEC, "results_190312.txt") 

    elif whichdate == '190317':
        data.pyraf_reduction("results_190317.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190317_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(2.0, RA, DEC, "results_190317.txt") 

    elif whichdate == '190318':
        data.pyraf_reduction("results_190318.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190318_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(2.0, RA, DEC, "results_190318.txt") 

    elif whichdate == '190326':
        data.pyraf_reduction("results_190326.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190326_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location, delta_x=-45.0) # telescope shifts after fullframes
        reduced_data.photometry(2.0, RA, DEC, "results_190326.txt") 

    elif whichdate == '190404':
        data.pyraf_reduction("results_190404.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190404_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(2.0, RA, DEC, "results_190404.txt") 

