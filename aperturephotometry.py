"""
@authors: Valérie Desharnais & Nicholas Vieira 
@aperturephotometry.py

Perform image segmentation to detect sources in the field and then search for 
a specific source. 
"""

def photometry(header, data, name, thresh_factor, RA_bound, DEC_bound, 
               results_file, im=True):
    """
    Input: the header of a reduced object's .fits file, the image data of the 
    file, the name to be used when creating the segmented image and a csv 
    containing all detected sources, a threshold factor to be used in image 
    segmentation, 2 arrays giving the RA and DEC (in degrees) boundaries on the 
    desired source, the name of the results textfile to which the photometry 
    will be appended, and a boolean indicating whether or not to save the image 
    of the segmentation test (optional; defaultTrue)
    Output: None
    
    Obtains a stack of images in the form of a header and data from a .fits 
    file. Estimates the background photon count of this image. Sets a threshold 
    above which we declare a cluster of pixels to be a source: this threshold 
    is defined as the background + the input thresh_factor*the RMS deviation of 
    the background. 
    
    Scans the input image and looks for sources which are at least 7 pixels in 
    area and above this threshold. Saves an image of the sources found in the 
    original field to the working directory. Uses a pixel coordinate to WCS 
    coordinate transformation, via a previously known reference pixel (see 
    PESTO_lib.WCS_merge()) to obtain WCS coords of all detected sources. 
    Outputs a .csv containing the properties of all detected sources.
    
    A tab-delimited results file is appended to. The number of images used in 
    the stack, the total exposure time of the stack and its error, the 
    timestamp (averaged across all images) and its error are already included.    
    IF a source is found, this script appends the x and y minima of the 
    source's centroid, the pixel area of the source, the photon count, and the 
    error on the photon count. If not, the flag NO SOURCE FOUND is appended to 
    the image. 
    """
    
    import numpy as np   
    import numpy.ma as ma 
    import os
    from astropy.stats import (SigmaClip, gaussian_fwhm_to_sigma, 
                               sigma_clipped_stats)
    from astropy.convolution import Gaussian2DKernel
    from astropy.table import Table, Column
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from photutils import Background2D, MedianBackground
    from photutils.utils import calc_total_error
    from photutils import detect_sources
    
    # perform 20 iterations of sigma clipping where needed 
    sigma_clip = SigmaClip(sigma=3.0, iters=20) 
    # background is estimated as the median of the entire image 
    bkg_estimator = MedianBackground() 
    mask = (data == 0) # mask all pixels where the ADU is 0  
    bkg = Background2D(data, (50,50), filter_size=(3,3), sigma_clip=sigma_clip, 
                       bkg_estimator=bkg_estimator, mask=mask)

    ### find sources using image segmentation
    
    # set the threshold for source detection 
    threshold = bkg.background + (thresh_factor*bkg.background_rms)
    sigma = 3.0*gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3.0, y_size=3.0)
    kernel.normalize()
    segm = detect_sources(data, threshold, npixels=7, filter_kernel=kernel)
    segm.remove_masked_labels(mask)
    
    try: 
        segm.remove_border_labels(10, partial_overlap=True, relabel=True)
    except: 
        print("The background threshold factor is too large; sources are "+
              "being ignored during image segmentation.\nPlease try a smaller"+
              " value.\n")
        return
 
    # pictures to see what's going on
    if(im):
        import matplotlib.pyplot as plt
        plt.switch_backend('agg') # stop matplotlib from trying to show image
        from astropy.visualization import SqrtStretch
        from astropy.visualization.mpl_normalize import ImageNormalize
        
        norm = ImageNormalize(stretch=SqrtStretch()) # normalize the image 
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(90, 90)) # 2 plots 

	     # show the (data-background) of the image with a greyscale colourmap 
        # using the above normalization
        ax1.imshow(data-bkg.background, origin='lower', cmap='Greys_r', 
                   norm=norm) 
        ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345)) 
        plt.savefig('segmentationtest_'+name+'.png')
    
    # find source properties (centroid, source pixel area, etc.) 
    from photutils import source_properties
    from astropy.wcs import WCS
    
    # calculate the error on the photon counts
#    tf = open('/data/irulan/omm_transients/'+results_file,'r')
#    contents = tf.readlines()
#    tf.close()
#    tf_last = contents[len(contents)-1]
#    tf_data = tf_last.split("\t")  
#    # gain*exposure*stack: 
#    effective_gain = 13.522*(float(tf_data[1])/1000.0)*float(tf_data[0])  
    effective_gain = 13.522
    # compute photon count error :
    error = calc_total_error(data, bkg.background_rms, effective_gain) 
    
    cat = source_properties(data-bkg.background, segm, wcs='all_pix2world', 
                            error=error)
    segm_tbl = cat.to_table() # contruct a table of source properties 
    
    # WCS object
    w = WCS(header)
    # get WCS of all sources, add to segm_tbl, and write the table
    ra, dec = w.all_pix2world(segm_tbl['xcentroid'], segm_tbl['ycentroid'],1)
    segm_tbl["ra"] = ra
    segm_tbl["dec"] = dec
    segm_tbl.write('segmentation_table_'+name+'.csv', format = 'csv', 
              overwrite=True)
    
    # build a new table with only the parameters we care about 
    tbl = Table()
    tbl["id"] = segm_tbl["id"] # id 
    tbl["xcentroid"] = segm_tbl["xcentroid"] # x coord
    tbl["ycentroid"] = segm_tbl["ycentroid"] # y coord 
    tbl["area"] = segm_tbl["area"] # area in pixels
    tbl["ra"] = ra # ra 
    tbl["dec"] = dec # dec
    tbl["pc"] = segm_tbl["source_sum"] # flux 
    tbl["pc_err"] = segm_tbl["source_sum_err"] # error on flux 
    tbl["mag_fit"] = -2.5*np.log10(tbl["pc"]) # instrumental magnitude
    tbl["mag_fit_unc"] = 2.5/(tbl["pc"]*np.log(10)) # error on magnitude 
    
    ### query Vizier to match sources and do aperture photometry
    
    # set the catalogue and filter of the image
    ref_catalog = "II/349/ps1"
    ref_catalog_name = "PS1" # PanStarrs 1
    filt = header["filtre"][0]
    # get the centre of the image and its RA, Dec
    x_size = data.shape[1]
    y_size = data.shape[0]
    ra_centre, dec_centre = np.array(w.all_pix2world(x_size/2.0, 
                                                     y_size/2.0, 1))
    # set radius to search in, minimum, maximum magnitudes
    minmag = 10.0
    maxmag = 22.0 
    max_emag = 0.3 # maximum error on magnitude 
    pixscale = np.mean(np.abs([header["CDELT1"], header["CDELT2"]])) # in deg
    pixscale = pixscale*3600.0 # in arcsec
    radius = pixscale*y_size/60.0 # radius in arcmin 
    
    # query print statement
    #print('Querying Vizier %s around RA %.4f, Dec %.4f with a radius of %.4f arcmin\n'%(
    #        ref_catalog, ra_centre, dec_centre, radius))
    
    # querying 
    v = Vizier(columns=["*"], column_filters={
            filt+"mag":str(minmag)+".."+str(maxmag),
            "e_"+filt+"mag":"<"+str(max_emag)}, row_limit=-1) # no row limit
    
    Q = v.query_region(SkyCoord(ra=ra_centre, dec=dec_centre, 
                    unit = (u.deg, u.deg)), radius = str(radius)+'m', 
                    catalog=ref_catalog, cache=False)
    cat_coords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
    # mask out edge sources
    x_lims = [int(0.05*x_size), int(0.95*x_size)] 
    y_lims = [int(0.05*y_size), int(0.95*y_size)]
    mask = (cat_coords[0] > x_lims[0]) & (
            cat_coords[0] < x_lims[1]) & (
            cat_coords[1] > y_lims[0]) & (
            cat_coords[1] < y_lims[1])
    good_cat_sources = Q[0][mask] # sources in catalogue 
    
    # cross-matching coords of sources found by astrometry
    source_coords = SkyCoord(ra=tbl['ra'], dec=tbl['dec'], frame='fk5', 
                             unit='degree')
    # and coords of valid sources in the queried catalog 
    cat_source_coords = SkyCoord(ra=good_cat_sources['RAJ2000'], 
                                     dec=good_cat_sources['DEJ2000'], 
                                     frame='fk5', unit='degree')
        
    # indices of matching sources (within 5.0 pix of each other) 
    idx_image, idx_cat, d2d, d3d = cat_source_coords.search_around_sky(
            source_coords, 5.0*pixscale*u.arcsec)
    
    # compute magnitude offsets and zero point
    mag_offsets = ma.array(good_cat_sources[filt+'mag'][idx_cat] - 
                      tbl['mag_fit'][idx_image])
    zp_mean, zp_med, zp_std = sigma_clipped_stats(mag_offsets) # zero point
    
    mag_calib = tbl['mag_fit'] + zp_mean # compute magnitudes 
    mag_calib.name = 'mag_calib'
    mag_calib_unc = np.sqrt(tbl['mag_fit_unc']**2 + zp_std**2) # propagate errs
    mag_calib_unc.name = 'mag_calib_unc'
    tbl['mag_calib'] = mag_calib
    tbl['mag_calib_unc'] = mag_calib_unc
    
    #print(zp_mean)
    #print(zp_std)
    
    # add flag indicating if source is in catalog
    #in_cat = []
    #for i in range(len(tbl)):
    #    if i in idx_image:
    #        in_cat.append(True)
    #    else:
    #        in_cat.append(False)
    #in_cat_col = Column(data=in_cat, name="in "+ref_catalog_name)
    #tbl["in "+ref_catalog_name] = in_cat_col
    
    #return tbl

    # boundaries on the desired source
    RA_min, RA_max = RA_bound
    DEC_min, DEC_max = DEC_bound

    # parse a list of all sources for a source within the RA, Dec bounds
    cwd = os.getcwd() # current working dir
    for i in range(len(tbl['id'])):
        # if source is found:
        if (RA_min <= tbl[i]['ra'] <= RA_max) and (
                DEC_min <= tbl[i]['dec'] <= DEC_max):                      
            print("\nFound a source.\n")

            # Write xcentroid and ycentroid, the pixel area of the source, 
            # the photon count, photon count error, calibrated magnitude, 
            # calibrated magnitude error, and filter used to the file.
            # If xcentroid and ycentroid change drastically from one stack to 
            # another, the sources are not the same or astrometric calibration 
            # may have failed. 
            xcentroid = tbl["xcentroid"].data[i]
            ycentroid = tbl["ycentroid"].data[i]
            area = tbl["area"].data[i]
            pc = tbl["pc"].data[i]
            pc_err = tbl["pc_err"].data[i] 
            mag = tbl["mag_calib"].data[i]
            mag_err = tbl["mag_calib_unc"][i]
            # line to write to the file: 
            line = str(xcentroid)+"\t"+str(ycentroid)+"\t"+str(area)
            line += "\t"+str(pc)+"\t"+str(pc_err) 
            line += "\t"+str(mag)+"\t"+str(mag_err)+"\t"+filt+"\n"
            tf = open(cwd+"/"+results_file,'a')
            tf.write(line)
            tf.close()
            return tbl

    # if no source is found:
    line = "NO SOURCE FOUND.\n" 
    print("No source found.\n")
    tf = open(cwd+"/"+results_file,'a')
    tf.write(line)
    tf.close()

    return tbl



    
