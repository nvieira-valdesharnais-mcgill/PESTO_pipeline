"""
@authors: Valérie Desharnais & Nicholas Vieira 
@aperturephotometry.py

Perform image segmentation to detect sources in the field and then search for a specific source. 
"""


def photometry(header, data, name, thresh_factor, RA_bound, DEC_bound, results_file, im=True):
    """
    Input: the header of a reduced object's .fits file, the image data of the file, the name 
    to be used when creating the segmented image and a csv containing all detected sources, a 
    threshold factor to be used in image segmentation, 2 arrays giving the RA and DEC (in degrees)
    boundaries on the desired source, the name of the results textfile to which the photometry will be 
    appended, and a boolean indicating whether or not to save the image 
    of the segmentation test (optional; default is True)
    Output: None
    
    Obtains a stack of images in the form of a header and data from a .fits file. Estimates the background photon count of this 
    image. Sets a threshold above which we declare a cluster of pixels to be a source: this threshold is defined as 
    the background + the input thresh_factor*the RMS deviation of the background. 
    [Should go into more detail here, or in the README.md]
    
    Scans the input image and looks for sources which are at least 7 pixels in area and above this threshold. Saves an image 
    of the sources found in the original field to the working directory. Uses a pixel coordinate to WCS coordinate 
    transformation, via a previously known reference pixel (see PESTO_lib.WCS_preparation()) to obtain WCS coords of all 
    detected sources. Outputs a .csv containing the properties of all detected sources.
    
    A tab-delimited results file is appended to. The number of images used in the stack, the total exposure 
    time of the stack and its error, the timestamp (averaged across all images) and its error are already included.    
    IF a source is found, this script appends the x and y minima of the source's centroid, the pixel area of the source, 
    the photon count, and the error on the photon count. If not, the flag NO SOURCE FOUND is appended to the image. 
    """
       
    # estimate background in the picture
    from astropy.io import fits
    from astropy.stats import SigmaClip
    from photutils import Background2D, MedianBackground
    from photutils.utils import calc_total_error
    from astropy.table import Table
    sigma_clip = SigmaClip(sigma=3.0, iters=20) # perform 20 iterations of sigma clipping where needed 
    bkg_estimator = MedianBackground() # background is estimated as the median of the entire image 
    mask = (data == 0) # mask all pixels where the ADU is 0 
    bkg = Background2D(data, (50,50), filter_size=(3,3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask)

    # find sources using image segmentation
    threshold = bkg.background + (thresh_factor*bkg.background_rms) # set the threshold for source detection 
    from astropy.convolution import Gaussian2DKernel
    from astropy.stats import gaussian_fwhm_to_sigma # fwhm ~ 2.355sigma for a Gaussian
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
 
    # pictures to see what's going on
    if(im):
        import numpy as np
        import matplotlib.pyplot as plt
        plt.switch_backend('agg') # to prevent matplotlib from trying to display images
        from astropy.visualization import SqrtStretch
        from astropy.visualization.mpl_normalize import ImageNormalize
        norm = ImageNormalize(stretch=SqrtStretch()) # normalize the image 
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(90, 90)) # 2 plots 

	# show the (data-background) of the image with a greyscale colourmap using the 
        # above normalization
        ax1.imshow(data-bkg.background, origin='lower', cmap='Greys_r', norm=norm) 
        ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345)) 
        plt.savefig('segmentationtest_'+name+'.png')
    
    # find source properties (centroid, source pixel area, etc.) 
    from photutils import source_properties
    from astropy.wcs import WCS
    
    # calculate the error on the photon counts
    tf = open('/data/irulan/omm_transients/'+results_file,'r')
    contents = tf.readlines()
    tf.close()
    tf_last = contents[len(contents)-1]
    tf_data = tf_last.split("\t")   
    effective_gain = 13.522*(float(tf_data[1])/1000.0)*float(tf_data[0])  # gain*exposure*stack
    error = calc_total_error(data, bkg.background_rms, effective_gain) # compute photon count error 
    
    cat = source_properties(data-bkg.background, segm, wcs='all_pix2world', error=error)
    tbl = cat.to_table() # contruct a table of source properties 
    
    #print("\n\nThe pixel coords of the sources found:\n",tbl['xcentroid'],tbl['ycentroid'])   # print the pixel coords of the sources
    coord = WCS(header)

    #print("\n\nThe WCS parameters to be used in WCS conversion:",coord)                      # print WCS objects's qualities
    xposition, yposition = coord.all_pix2world(tbl['xcentroid'], tbl['ycentroid'],1)
    #print("\n\nThe converted WCS coords:\n",xposition,yposition)                             # print the converted WCS coords
    tbl['xcentroid'] = xposition
    tbl['ycentroid'] = yposition
    tbl.write('source_table_'+name+'.csv', format = 'csv', overwrite=True)

    # boundaries on the desired source
    RA_min = RA_bound[0]
    RA_max = RA_bound[1]
    DEC_min = DEC_bound[0]
    DEC_max = DEC_bound[1] 

    # parse a list of all sources for a source within the RA, Dec bounds 
    for i in range(len(tbl['id'])):
        if (RA_min <= tbl[i]['xcentroid'] <= RA_max) and (DEC_min <= tbl[i]['ycentroid'] <= DEC_max):                      
            print("\nFound a source.\n")
            new_tbl = Table(rows=tbl[i])
            # print(new_tbl) # prints source properties of the one source 

            # write the lower left x and y bounds, the pixel area of the source, the photon count, and the 
            # photon count error to the results file
            import pandas as pd
            csv_file = pd.read_csv('source_table_'+name+'.csv') # read
            x_min = csv_file.iloc[i,10] # if changes drastically from set to set, WCS coord error
            y_min = csv_file.iloc[i,12] # if changes drastically from set to set, WCS coord error
            pc = csv_file.iloc[i,5]
            pc_err = csv_file.iloc[i,6] 
            area = csv_file.iloc[i,20]
            # line to write to the file: 
            line = str(x_min)+"\t"+str(y_min)+"\t"+str(area)+"\t"+str(pc)+"\t"+str(pc_err)+"\n" 
            tf = open('/data/irulan/omm_transients/'+results_file,'a')
            tf.write(line)
            tf.close()
            return None

    # if no source is found
    # line to write to the file 
    line = "NO SOURCE FOUND.\n" 
    print("No source found.\n")
    tf = open('/data/irulan/omm_transients/'+results_file,'a')
    tf.write(line)
    tf.close()

    return None



    
