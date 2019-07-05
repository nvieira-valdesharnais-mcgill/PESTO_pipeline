#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 15:21:54 2019

@authors: Nicholas Vieira & Valérie Desharnais
@generalpurpose.py

"""
import PESTO_lib 
import subprocess
import re

whichdate = 'x'
file = 'x'
stack = 'x'

while not(whichdate in ['180709','180928','190312','190317','190318', 
                        '190326', '190404']):
    whichdate = input("Which date do you want? Enter 180709, 180928, "+
                      "190312, 190317, 190318, 190326 or 190404 \n> ")
while (file == 'x'):
    file = input("\nWhat is the integer starting point? Pick using the "+
                 "table below:\n\
    +----------+--------+--------+\n\
    | JULY 09  |  96931 | 528554 |\n\
    +----------+--------+--------+\n\
    | SEPT 28  |  30054 | 463776 |\n\
    +----------+--------+--------+\n\
    | MARCH 12 |   8241 | 334055 |\n\
    +----------+--------+--------+\n\
    | MARCH 17 |  42805 | 584457 |\n\
    +----------+--------+--------+\n\
    | MARCH 18 |   9878 | 565912 |\n\
    +----------+--------+--------+\n\
    | MARCH 26 |   2200 | 437399 |\n\
    +----------+--------+--------+\n\
    | APRIL 04 |   8235 | 613404 |\n\
    +----------+--------+--------+\n\n> ")                 
    
stack = int(input("\nWhat is the stack size? \n> "))

file = re.sub(whichdate+"*_", "", file) # get rid of date at filname start
#while file[0] == "0":
#    file = file[1:] # get rid of unecessary leading 0s
file_num = int(file)

name = "mini_reduced_data" 

###### FORMAT
# JULY: 096931 to 528554
# SEPT: 030054 to 463776
# MARCH 12: 008241 to 334055
# MARCH 17: 042805 to 584457
# MARCH 18: 009878 to 565912
# MARCH 26: 002200 to 437399
# APRIL 04: 008235 to 613404


if whichdate == '180709':
    location = ["/exports/scratch/MAXIJ1820/180709_works","/exports/scratch/MAXIJ1820/180709_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/180709/Target/MAXIJ1820+070/180709"
    bash_copy_end = " /exports/scratch/MAXIJ1820/180709_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/180709_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/180709_works/*.fits.gz"
    limit = 528554
elif whichdate == '180928':
    location = ["/exports/scratch/MAXIJ1820/180928_works","/exports/scratch/MAXIJ1820/180928_calibs"] 
    bash_copy_start = "cp /exports/scratch/MAXIJ1820/180928/MAXI1820+070-180928-OMM/Target/180928"
    bash_copy_end = " /exports/scratch/MAXIJ1820/180928_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/180928_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/180928_works/*.fits.gz"
    limit = 463776
elif whichdate == '190312':
    location = ["/exports/scratch/MAXIJ1820/190312_works","/exports/scratch/MAXIJ1820/190312_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190312/MAXIJ1820+070/190312"
    bash_copy_end = " /exports/scratch/MAXIJ1820/190312_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190312_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190312_works/*.fits.gz"
    limit = 334055
elif whichdate == '190317':
    location = ["/exports/scratch/MAXIJ1820/190317_works","/exports/scratch/MAXIJ1820/190317_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190317/MAXIJ1820+070/190317"
    bash_copy_end = " /exports/scratch/MAXIJ1820/190317_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190317_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190317_works/*.fits.gz"
    limit = 584457
elif whichdate == '190318':
    location = ["/exports/scratch/MAXIJ1820/190318_works","/exports/scratch/MAXIJ1820/190318_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190318/MAXIJ1820+070/190318_0000"
    bash_copy_end = " /exports/scratch/MAXIJ1820/190318_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190318_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190318_works/*.fits.gz"
    limit = 565912
elif whichdate == '190326':
    location = ["/exports/scratch/MAXIJ1820/190326_works","/exports/scratch/MAXIJ1820/190326_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190326/MAXIJ1820+070/190326"
    bash_copy_end = " /exports/scratch/MAXIJ1820/190326_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190326_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190326_works/*.fits.gz"
    limit = 437399
elif whichdate == '190404':
    location = ["/exports/scratch/MAXIJ1820/190404_works","/exports/scratch/MAXIJ1820/190404_calibs"] 
    bash_copy_start = "cp /data/irulan/omm_transients/MAXIJ1820/190404/MAXIJ1820+070/190404"
    bash_copy_end = " /exports/scratch/MAXIJ1820/190404_works/"
    bash_empty = "rm -rf /exports/scratch/MAXIJ1820/190404_works/*"
    bash_gunzip = "gunzip /exports/scratch/MAXIJ1820/190404_works/*.fits.gz"
    limit = 613404

target = "/exports/scratch/MAXIJ1820"
imtype = ["object","calibration"]

# Here, select the RA and DEC bounds, depending on whether we're looking at 
# MAXI or an alt source
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

start = file_num
end = start+stack

while start+stack < limit: # will run until it can't anymore 
    subprocess.run(bash_empty,shell=True)
    i = start 
    while i < end:     
        cp_cmd = bash_copy_start+"*000"+str(i)+".fits* "+bash_copy_end
        subprocess.run(cp_cmd, shell=True)
        #print(i)
        i += 1
    
    subprocess.run(bash_gunzip,shell=True)
    
    # MAXI STUFF
    data = PESTO_lib.raw_PESTO_data(location,imtype,name,target) 
    data.flush()
    data.produce_lists() 
    data.make_working_directory() 

    if whichdate == '180709':
        data.pyraf_reduction("results_july.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        reduced_data.WCS_preparation(42.0934)     # old WCS technique
        reduced_data.photometry(RA, DEC, 2.0, "results_july.txt")
  
    elif whichdate == '180928':
        data.pyraf_reduction("results_sept.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        reduced_data.WCS_preparation(40.3927)     # old WCS technique
        reduced_data.photometry(RA, DEC, 2.0, "results_sept.txt")

    elif whichdate == '190312':
        data.pyraf_reduction("results_190312.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190312_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(RA, DEC, 2.0, "results_190312.txt") 

    elif whichdate == '190317':
        data.pyraf_reduction("results_190317.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190317_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(RA, DEC, 2.0, "results_190317.txt") 

    elif whichdate == '190318':
        data.pyraf_reduction("results_190318.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190318_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(RA, DEC, 2.0, "results_190318.txt") 

    elif whichdate == '190326':
        data.pyraf_reduction("results_190326.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190326_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location, delta_x=45.0) # telescope shifts after taking fullframes
        reduced_data.photometry(RA, DEC, 2.0, "results_190326.txt") 

    elif whichdate == '190404':
        data.pyraf_reduction("results_190404.txt") 
        reduced_data = data.extract_reduced_images("/exports/scratch/MAXIJ1820/"+name, "minired") 
        wcs_location = "/data/irulan/omm_transients/wcs_solutions/190404_soln.fits"    # new WCS technique
        reduced_data.WCS_merge(wcs_location)
        reduced_data.photometry(RA, DEC, 2.0, "results_190404.txt")

    
    start = end # new beginning
    end += stack # new end 
    

