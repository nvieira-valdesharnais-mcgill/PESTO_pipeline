def photometry(header, data, name, im, thresh_factor):
    # Added the thresh_factor to make function more useful. Argument im is currently not used.    
    #Estimates background in the picture
    from astropy.io import fits
    from astropy.stats import SigmaClip
    from photutils import Background2D, MedianBackground
    from photutils.utils import calc_total_error
    from astropy.table import Table
    sigma_clip = SigmaClip(sigma=3.0, iters=20)
    bkg_estimator = MedianBackground()
    mask = (data == 0) 
    bkg = Background2D(data, (50,50), filter_size=(3,3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask)

    #Finds sources using image segmentation
    threshold = bkg.background + (thresh_factor*bkg.background_rms) #originally 5.0
    from astropy.convolution import Gaussian2DKernel
    from astropy.stats import gaussian_fwhm_to_sigma
    from photutils import detect_sources
    sigma = 3.0*gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3.0, y_size=3.0)
    kernel.normalize()
    segm = detect_sources(data, threshold, npixels=7, filter_kernel=kernel)
    segm.remove_masked_labels(mask)
    try: 
        segm.remove_border_labels(10, partial_overlap=True, relabel=True)
    except: 
        print("The background threshold factor is too large; sources are being ignored during image segmentation.\nPlease try a smaller value.\n")
        return
   
 
    #Pictures to see what's going on
    if(im):
        import numpy as np
        import matplotlib.pyplot as plt
        plt.switch_backend('agg')
        from astropy.visualization import SqrtStretch
        from astropy.visualization.mpl_normalize import ImageNormalize
        norm = ImageNormalize(stretch=SqrtStretch()) #added this
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(90, 90))
        ax1.imshow(data-bkg.background, origin='lower', cmap='Greys_r', norm=norm) 
        #vmin = 0, vmax = 0.50*np.amax(data))
        ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345)) 
        plt.savefig('segmentationtest_'+name+'.png')
    

    #Find source properties and convert their position to WCS
    from photutils import source_properties
    from astropy.wcs import WCS
    
    #calculate the error on the photon count first before building table
    tf = open('/data/irulan/omm_transients/results.txt','r')
    contents = tf.readlines()
    tf.close()
    tf_last = contents[len(contents)-1]
    tf_data = tf_last.split("\t")   
    effective_gain = 13.522*(float(tf_data[1])/1000.0)*float(tf_data[0])  #gain*exposure*stack
    error = calc_total_error(data, bkg.background_rms, effective_gain)
    
    cat = source_properties(data-bkg.background, segm, wcs='all_pix2world', error=error)
    tbl = cat.to_table()

    pix_x = tbl['xcentroid'].value
    pix_y = tbl['ycentroid'].value
    
    #print("\n\nThe pixel coords of the sources found:\n",tbl['xcentroid'],tbl['ycentroid'])   # Print the pixel coords of the sources
    coord = WCS(header)

    #print("\n\nThe WCS parameters to be used in WCS conversion:",coord)                      # Print WCS objects's qualities
    xposition, yposition = coord.all_pix2world(tbl['xcentroid'], tbl['ycentroid'],1)
    #print("\n\nThe converted WCS coords:\n",xposition,yposition)                              # Print the converted WCS coords
    tbl['xcentroid'] = xposition
    tbl['ycentroid'] = yposition
    tbl.write('source_table_'+name+'.csv', format = 'csv', overwrite=True)

    maxi_RA_min = 275.085 # minimum RA we will allow for MAXI
    maxi_RA_max = 275.097 # maximum RA we will allow for MAXI
    maxi_DEC_min = 7.1849 # minimum DEC we will allow for MAXI
    maxi_DEC_max = 7.1860 # maximum DEC we will allow for MAXI

    # This source worked for both, so leaving it in
    # Source 3 in Nick's finderchart, used with threshfactor 3
    #RA_min = 275.08 # minimum RA for alternate source
    #RA_max = 275.09 # maximum RA for alternate source
    #DEC_min = 7.180 # minimum DEC
    #DEC_max = 7.184 # maximum DEC
    
    # Source 1 (ref pix) in Nick's finderchart, used with threshfactor 2
    # Worked in July
    RA_min = 275.128
    RA_max = 275.138
    DEC_min = 7.1398
    DEC_max = 7.1410

    # Source 2 in Nick's finderchart, used with threshfactor 2   
    #Worked in July
    #RA_min = 275.1098
    #RA_max = 275.1106
    #DEC_min = 7.1655
    #DEC_max = 7.1725

    #print("\nNow, parsing the objects found by image segmentation for MAXI:\n")
    for i in range(len(tbl['id'])):
        #print("RA="+str(tbl[i]['xcentroid'])+", DEC="+str(tbl[i]['ycentroid'])) # prints RA, DEC of objects found by segmenting
        #if (maxi_RA_min <= tbl[i]['xcentroid'] <= maxi_RA_max) and (maxi_DEC_min <= tbl[i]['ycentroid'] <= maxi_DEC_max):  # MAXI
        if (RA_min <= tbl[i]['xcentroid'] <= RA_max) and (DEC_min <= tbl[i]['ycentroid'] <= DEC_max):                       # alt source
            print("\nFound a source.\n")
            new_tbl = Table(rows=tbl[i])
            #print(new_tbl) # prints just the info on MAXI/the alt source
        
            #append the pixel x-min, pixel y-min, photon count, photon count error and pixel area to results.txt
            #only photon count and photon count error are essential; others are just for debugging
            import pandas as pd
            csv_file = pd.read_csv('source_table_'+name+'.csv')
            x_min = csv_file.iloc[i,10] # if changes drastically from set to set, WCS coord error
            y_min = csv_file.iloc[i,12] # if changes drastically from set to set, WCS coord error
            pc = csv_file.iloc[i,5]
            pc_err = csv_file.iloc[i,6] 
            area = csv_file.iloc[i,20]
            line = str(x_min)+"\t"+str(y_min)+"\t"+str(area)+"\t"+str(pc)+"\t"+str(pc_err)+"\n" # write everything
            #line = str(pc)+"\t"+str(pc_err)+"\n"                                           # write the essentials
            tf = open('/data/irulan/omm_transients/results.txt','a')
            tf.write(line)
            tf.close()
            return None

    # if no MAXI is found
    line = "NO SOURCE FOUND.\n" 
    print("No source found.\n")
    print("\a") # make a beep to alert the user 

    tf = open('/data/irulan/omm_transients/results.txt','a')
    tf.write(line)
    tf.close()

    return None



    
