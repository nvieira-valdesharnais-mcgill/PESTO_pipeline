def clean_no_source(tf_path,output_path=-1):
    """
    Removes lines from the data file containing the NO SOURCE error message.
    """
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    #if no output location given, store as home/smooth.txt
    if output_path == -1:
            from os.path import expanduser
            output_path = expanduser("~")
            output_path = output_path+"clean.txt"    
    if tf_path == output_path:
        output = open(output_path,"a")
        output.write("Cleaned data (NO SOURCES error removed) appended below:"+"\n")
        output.close()    
    for line in contents:
        if "SOURCE" not in line:
            output = open(output_path,"a")
            output.write(line)
            output.close()

def clean_broken_lines(tf_path,output_path):
    """
    Removes lines from the data file containing less than 10 entries or more than 11 entries.
    In other words, only full lines (10 tab entries) with one potential comment (1 tab entry) are left.
    """
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    for line in contents:
        data=line.split("\t")
        #print(len(data))
        if len(data) > 9 and len(data) < 12:
            output = open(output_path,"a")
            output.write(line)
            output.close()
            
def standardize_stacking(tf_path,output_path,stack):
    """
    Removes lines from the data file containing an incorrect number of stacked images.
    """
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    for line in contents:
        data=line.split("\t")
        if abs(float(data[0])-stack)<0.0001:
            output = open(output_path,"a")
            output.write(line)
            output.close()    

def set_initial_time_zero(tf_path,output_path):
    """
    Sets the first timestamp to be t=0 and adjusts all subsequent timestamps.
    Will warn you and break if a day transition is detected
    """    
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()        
    initial_line = contents[0].split("\t")
    initial_time = float(initial_line[3])
    old_timestamp = initial_time
    for line in contents:
        data = line.split("\t")
        if abs(old_timestamp - float(data[3])) > 85000:
            print("Warning: day transition in observational data detected.")
            print("Run correct_day_transition() on this data before modifying.")
            return
        else:
            old_timestamp = float(data[3])        
        time=float(data[3])-initial_time
        new_line = data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+str(time)+"\t"+data[4]+"\t"+data[5]+"\t"+data[6]+"\t"+data[7]+"\t"+data[8]+"\t"+data[9]#+"\n"
        output = open(output_path,"a")
        output.write(new_line)
        output.close()
            
def correct_day_transition(tf_path,output_path,lower_limit=200,upper_limit=86000):
    """
    Fixes the timestamps if a day transition occured during observation time.
    Can adjust sensitivity if needed to detect the transition by modifying the limits.
    Note: will not automatically set the first timestamp to be t=0.
    """
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    
    times, starts, ends = [],[],[]
    for line in contents:
        data = line.split("\t")
        times.append(float(data[3]))
        start = data[0]+"\t"+data[1]+"\t"+data[2]
        end = data[4]+"\t"+data[5]+"\t"+data[6]+"\t"+data[7]+"\t"+data[8]+"\t"+data[9] 
        starts.append(start)
        ends.append(end)
        
    index=-1
    for i in range(len(times)-1):
        if times[i]>upper_limit and times[i+1]<lower_limit:
            index = i+1
            break
    if index < 0:
        print("No transition detected in this dataset. Exiting.")
    else:
        for i in range(index+1,len(times)):
            times[i]=times[i]+86400
    
    for i in range(len(contents)):
        line = starts[i]+"\t"+str(times[i])+"\t"+ends[i]+"\n"
        output = open(output_path,"a")
        output.write(line)
        output.close()
        
def sort_by_timestamp(tf_path,output_path):
    """
    Sort the data by timestamp. Warns if a day transition occurs in the data. 
    """
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()    
    contents_to_sort = []
    initial_line = contents[0].split("\t")
    old_timestamp = float(initial_line[3])
    for line in contents:
        data = line.split("\t")
        if abs(old_timestamp - float(data[3])) > 85000:
            print("Warning: day transition in observational data detected.")
            print("Run correct_day_transition() on this data before sorting.")
            return
        else:
            old_timestamp = float(data[3])
        entry= [float(data[3]),line]
        contents_to_sort.append(entry)    
    contents_sorted = quicksort(contents_to_sort)    
    for entry in contents_sorted:
        output = open(output_path,"a")
        output.write(entry[1])
        output.close()
        
def quicksort(contents_to_sort):
    """
    Implementation of quicksort for sort_by_timestamp()
    """
    if len(contents_to_sort) < 2:
        return contents_to_sort
    pivot_index = int(len(contents_to_sort)/2.0)
    pivot = contents_to_sort[pivot_index][0]
    less, more = [], []
    for i in range(len(contents_to_sort)):
        if i is pivot_index:
            continue
        if contents_to_sort[i][0] <= pivot:
            less.append(contents_to_sort[i])
        else:
            more.append(contents_to_sort[i])
    contents_sorted = []  
    if len(less) > 0:
        less_sorted = quicksort(less)
        contents_sorted.extend(less_sorted)
    contents_sorted.append(contents_to_sort[pivot_index])
    if len(more) > 0:
        more_sorted = quicksort(more)
        contents_sorted.extend(more_sorted) 
    return contents_sorted

def correct_outliers(tf_path,output_path,sig=3.0):
    
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    photon_counts = []
    for line in contents:
        data = line.split("\t")
        photon_counts.append(float(data[8]))
    
    photon_count_differences = []
    for i in range(len(photon_counts)-1):
        photon_count_differences.append(photon_counts[i+1]-photon_counts[i])
        
    import statistics 
    av = statistics.mean(photon_count_differences)
    stdev = statistics.stdev(photon_count_differences)
    lower_lim = av-sig*stdev
    upper_lim = av+sig*stdev
    
    output = open(output_path,"a")
    output.write(contents[0])
    output.close() 
    counter = 0
    for i in range(len(photon_counts)-1):
        if photon_count_differences[i] > lower_lim and photon_count_differences[i] < upper_lim:
            output = open(output_path,"a")
            output.write(contents[i+1])
            output.close()
        else:
            counter = counter+1
    print("Lines removed: "+str(counter))

def strip_above(tf_path,output_path,strip=100000,tmin=0,tmax=100000):
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    for line in contents:
        data = line.split("\t")
        if float(data[3]) < tmin or float(data[3]) > tmax:
            output = open(output_path,"a")
            output.write(line)
            output.close()
        elif float(data[8]) < strip:
            output = open(output_path,"a")
            output.write(line)
            output.close()
        
def smooth(tf_path,output_path=-1,threshold=-1,factor=1.025):
    """
    Smooths data by averaging a bin by both its predecessor and successor.
    Removes lines where the time gap is above the threshold value.
    Threshold default set to 102.5% of the average smoothed timestamp differences.
    """
    
    #acquiring the data in the text file to average
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    d0,d1,d2,d3,d4,d5,d6,d7,d8,d9 = ([] for i in range(10))
    for line in contents:
        data = line.split("\t")
        d0.append(float(data[0]))
        d1.append(float(data[1]))
        d2.append(float(data[2]))
        d3.append(float(data[3]))
        d4.append(float(data[4]))
        d5.append(float(data[5]))
        d6.append(float(data[6]))
        d7.append(float(data[7]))
        d8.append(float(data[8]))
        d9.append(float(data[9]))

    #naively smoothing the data by averaging the predecessor, the bin and the succesor together    
    nsd0,nsd1,nsd2,nsd3,nsd4,nsd5,nsd6,nsd7,nsd8,nsd9 = ([] for i in range(10))
    for i in range(len(d0)-2):
        nsd0.append((d0[i]+d0[i+1]+d0[i+2])/3.0)
        nsd1.append((d1[i]+d1[i+1]+d1[i+2])/3.0)
        nsd2.append((d2[i]+d2[i+1]+d2[i+2])/3.0)
        nsd3.append((d3[i]+d3[i+1]+d3[i+2])/3.0)
        nsd4.append((d4[i]+d4[i+1]+d4[i+2])/3.0)
        nsd5.append((d5[i]+d5[i+1]+d5[i+2])/3.0)
        nsd6.append((d6[i]+d6[i+1]+d6[i+2])/3.0)
        nsd7.append((d7[i]+d7[i+1]+d7[i+2])/3.0)
        nsd8.append((d8[i]+d8[i+1]+d8[i+2])/3.0)
        nsd9.append((d9[i]+d9[i+1]+d9[i+2])/3.0)
        
    #acquiring the time differences between each subsequent pair of smoothed data points    
    t_dif = []
    for i in range(len(nsd3)-1):
        t_dif.append(nsd3[i+1]-nsd3[i])

    #defining a threshold if no argument was passed    
    if threshold == -1:
        threshold=factor*sum(t_dif)/float(len(t_dif))
        print("The threshold for gaps was set to about "+("%.2f"%threshold)+", about "+("%.1f"%(100*factor))+"% of the average time spacing."+"\n")

    #if the output location = input location, append a delimitation.
    if tf_path == output_path:
        output = open(output_path,"a")
        output.write("Smoothed data of threshold "+("%.2f"%threshold)+" appended below:"+"\n")
        output.close()

    #correcting the naively smoothed data by appending only entries with time differences less than the threshold
    #if no output location given, store as home/smooth.txt
    if output_path == -1:
            from os.path import expanduser
            output_path = expanduser("~")
            output_path = output_path+"smooth.txt"
    line = str(nsd0[0])+"\t"+str(nsd1[0])+"\t"+str(nsd2[0])+"\t"+str(nsd3[0])+"\t"+str(nsd4[0])+"\t"+str(nsd5[0])+"\t"+str(nsd6[0])+"\t"+str(nsd7[0])+"\t"+str(nsd8[0])+"\t"+str(nsd9[0])+"\n"
    output = open(output_path,"a")
    output.write(line)
    output.close() 
    counter = 0    
    for i in range(len(t_dif)):
        if t_dif[i] <= threshold:
            counter = counter+1
            output = open(output_path,"a")
            line = str(nsd0[i+1])+"\t"+str(nsd1[i+1])+"\t"+str(nsd2[i+1])+"\t"+str(nsd3[i+1])+"\t"+str(nsd4[i+1])+"\t"+str(nsd5[i+1])+"\t"+str(nsd6[i+1])+"\t"+str(nsd7[i+1])+"\t"+str(nsd8[i+1])+"\t"+str(nsd9[i+1])+"\n"
            output.write(line)
            output.close()            
    print("The original file contained "+str(len(d0))+" entries. The smooth file contains "+str(counter)+" entries, a reduction of about "+("%.0f"%(100*(1-(counter/float(len(d0))))))+"%."+"\n")
            
def build_weather_database(db_path,tf_paths):
    """
    Builds a weather database for the correct_weather() function.
    """    
    
    for i in range(len(tf_paths)):
        # read in the data to be used to build the database
        tf = open(tf_paths[i],"r")
        contents_alt = tf.readlines()
        tf.close()
        # build the entry set of the database
        t, x, y, pc = [], [], [], []
        pc_largest = 0
        for line in contents_alt:
            data = line.split("\t")
            t.append(float(data[3]))
            x.append(float(data[5]))
            y.append(float(data[6]))
            pc.append(float(data[8]))
            if float(data[8]) > pc_largest:
                pc_largest = float(data[8])
        f = [pc_largest/p for p in pc]
        # write the entry set into the database file
        for j in range(len(t)):
            output = open(db_path,"a")
            line = str(t[j])+"\t"+str(x[j])+"\t"+str(y[j])+"\t"+str(pc[j])+"\t"+str(f[j])+"\n"
            output.write(line)
            output.close() 
        # indicate switching to next entry set in database
        if i is not (len(tf_paths)-1):
            output = open(db_path,"a")
            output.write("ALT_CHANGE\n")
            output.close()    

def correct_weather(db_path, tf_path, output_path, weighted=False, gap_threshold=3.0):     
    """
    Corrects a source for weather effects via factor averaging from alternative sources.
    Requires a weather database containing at least one alternative source to function.
    """           
    
    #load in the weather database, check the number of alts present
    tf = open(db_path,"r")
    contents_db = tf.readlines()
    tf.close()  
    alts = []
    ts,xs,ys,pcs,fs = [],[],[],[],[]
    for line in contents_db:        
        #read the line
        data = line.split("\t")        
        #if ALT is there, skip this line, append the completed alt, reset alt
        if 'ALT' in data[0]:
            alts.append([ts,xs,ys,pcs,fs])
            ts,xs,ys,pcs,fs = [],[],[],[],[]
        #else, attach this entry to the current alt
        else:
            ts.append(float(data[0]))
            xs.append(float(data[1]))
            ys.append(float(data[2]))
            pcs.append(float(data[3]))
            fs.append(float(data[4]))
    alts.append([ts,xs,ys,pcs,fs])            
    print(str(len(alts))+" alternate sources were found within the database.")
    # results in a 3D array: alts[alt number][data type][ith entry]
    
    #load in the source we are interested in correcting using the weather database
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    ts,xs,ys,pcs,pces = [],[],[],[],[]    
    line_segments = []
    for line in contents:
        data = line.split("\t")
        ts.append(float(data[3]))
        xs.append(float(data[5]))
        ys.append(float(data[6]))
        pcs.append(float(data[8]))
        pces.append(float(data[9]))
        line_segments.append(data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+data[3]+"\t"+data[4]+"\t"+data[5]+"\t"+data[6]+"\t"+data[7])
    source = [ts,xs,ys,pcs,pces]
    
    #we now have the alts and the source structured similarly
    #the source however has no data type <factors>, which we now need to interpolate    
    #we will save time by avoiding checking the weighted condition inside the loop
    
    if weighted is False:
 
        for j in range(len(source[0])):
            #print(j)
            alt_factors = []            
            for alt in alts:
                for i in range(len(alt[0])-1):                    
                    if alt[0][i] <= source[0][j] and alt[0][i+1] > source[0][j]:
                        #must be in <= and > state, which we can interpolate with if the gap is not extreme                       
                        if abs(alt[0][i+1]-alt[0][i]) <= gap_threshold:
                            m = (alt[4][i+1]-alt[4][i])/(alt[0][i+1]-alt[0][i])
                            b = alt[4][i]-m*alt[0][i]
                            alt_factors.append(m*source[0][j]+b) 
                            #print("successful interpolation")
                            break
                        else:   
                            #print("unsuccessful interpolation")
                            break                            
                    elif alt[0][i] > source[0][j] and alt[0][i+1] > source[0][j]:
                        #no relevant data from this alt for this source timestamp
                        #print("no factor for this alt at this timestamp")
                        break
            if len(alt_factors) < 1:
                #no alts could be used for correcting, skip this source time
                #print("no factors collected for this timestamp, correction failed")
                continue
            
            #correct the source photon count and photon count error
            source_factor = sum(alt_factors)/len(alt_factors)
            corrected_pc = source_factor*source[3][j]
            corrected_pc_error = source_factor*source[4][j]
            
            #output the new corrected line to the file
            #print("successful correction")
            line = line_segments[j]+"\t"+str(corrected_pc)+"\t"+str(corrected_pc_error)+"\n"
            output = open(output_path,"a")
            output.write(line)
            output.close()

    elif weighted is True:
        
        for j in range(len(source[0])):
            #print(j)
            alt_distances, alt_factors = [],[] 
            for alt in alts:
                for i in range(len(alt[0])-1):                    
                    if alt[0][i] <= source[0][j] and alt[0][i+1] > source[0][j]:
                        #must be in <= and > state, which we can interpolate with if the gap is not extreme                       
                        if abs(alt[0][i+1]-alt[0][i]) <= gap_threshold:
                            x = (alt[2][i]+alt[2][i+1])/2.0
                            y = (alt[3][i]+alt[3][i+1])/2.0
                            alt_distances.append(((source[2][j]-x)**2 + (source[3][j]-y)**2)**0.5)
                            m = (alt[4][i+1]-alt[4][i])/(alt[0][i+1]-alt[0][i])
                            b = alt[4][i]-m*alt[0][i]
                            alt_factors.append(m*source[0][j]+b) 
                            #print("successful interpolation")
                            break
                        else:
                            #print("unsuccessful interpolation")
                            break                            
                    elif alt[0][i] > source[0][j] and alt[0][i+1] > source[0][j]:
                        #no relevant data from this alt for this source timestamp
                        #print("no factor for this alt at this timestamp")
                        break

            if len(alt_factors) < 1:
                #no alts could be used for correcting, skip this source time
                #print("no factors collected for this timestamp, correction failed")
                continue
            
            #calculate the weighted average
            sum_weights = 0
            alt_weights = []
            for d in alt_distances:
                alt_weights.append(1/d)
                sum_weights = sum_weights+1/d            
            source_factor = 0
            for i in range(len(alt_weights)):
                source_factor = source_factor+(alt_weights[i]/sum_weights)*alt_factors[i]
                
            #correct the source photon count and photon count error
            corrected_pc = source_factor*source[3][j]
            corrected_pc_error = source_factor*source[4][j]
            
            #output the new corrected line            
            #print("successful correction")
            line = line_segments[j]+"\t"+str(corrected_pc)+"\t"+str(corrected_pc_error)+"\n"
            output = open(output_path,"a")
            output.write(line)
            output.close()
            
    print("Weather correction completed.")

def relative_photometry(threshold,source_path,alt_path,output_path):    
    """
    Creates a new file of a (main) source's photon count adjusted relative to another (alt) source's counts.
    """
    
    #read in source data
    sf = open(source_path,"r")
    contents_source = sf.readlines()
    sf.close()
    ts_s, pcs_s, pces_s = [],[],[]
    line_start, line_mid = [],[]
    for line in contents_source:
        data = line.split("\t")
        ts_s.append(float(data[3]))
        pcs_s.append(float(data[8]))
        pces_s.append(float(data[9]))
        line_start.append(str(float(data[0]))+"\t"+str(float(data[1]))+"\t"+str(float(data[2])))
        line_mid.append(str(float(data[4]))+"\t"+str(float(data[5]))+"\t"+str(float(data[6]))+"\t"+str(float(data[7])))
    
    #read in constant data
    af = open(alt_path,"r")
    contents_alt = af.readlines()
    af.close()
    ts_a, pcs_a, pces_a = [], [], []
    for line in contents_alt:
        data = line.split("\t")
        ts_a.append(float(data[3]))
        pcs_a.append(float(data[8]))
        pces_a.append(float(data[9]))
    
    #do the relative photometry
    cs, ca = 0, 0
    for i in range(len(ts_s)+len(ts_a)): 
        if abs(ts_s[cs]-ts_a[ca]) < threshold: #if within threshold, calculate relative photon count and error
            #error_on_div = (pcs_s[cs]/pcs_a[ca])*((pces_s[cs]/pcs_s[cs])**2+(pces_a[ca]/pcs_a[ca])**2)**0.5
            #for now, we do source_count - alt_count and source_error + alt_error
            line = line_start[cs]+"\t"+str(ts_s[cs])+"\t"+line_mid[cs]+"\t"+str(pcs_s[cs]-pcs_a[ca])+"\t"+str(pces_s[cs]+pces_a[ca])+"\n"
            output = open(output_path,"a")
            output.write(line)
            output.close()
            if ts_s[cs] > ts_a[ca] or (cs == len(ts_s)-1):
                ca = ca+1
            else:
                cs = cs+1
        else: #else, increment lowest counter
            if ts_s[cs] > ts_a[ca] or (cs == len(ts_s)-1):
                ca = ca+1
            else:
                cs = cs+1
        if (ca == len(ts_a)-1) and (cs == len(ts_s)-1): #break if end of source or alt reached
            break
            
def correct_poisson(freqs, powers, fit_lows=True):
    """
    Correct Poisson noise from a periodogram by fitting a decaying exponential to the data.
    Default is to use only local minima during the fitting, breaks if less than two are found.
    """
    # find local minima in powers if fit_lows set to True (default)
    min_freqs, min_powers = [],[]
    if fit_lows is True:
        for i in range(len(freqs)-2):
            if powers[i+2] > powers[i+1] and powers[i] > powers[i+1]:
                min_freqs.append(freqs[i+1])
                min_powers.append(powers[i+1])
        if len(min_freqs) < 2:
            print("One or no local minima were found, please set fit_lows to False during correction.")
            return
    else:
        min_freqs,min_powers = freqs,powers
            
    #fit decaying exponential to power minima
    import numpy as np
    from scipy.optimize import curve_fit
    def func(x, a, b):
        return a*np.exp(b*x)
    popt, pcov = curve_fit(func, xdata=min_freqs, ydata=min_powers,maxfev=1000,bounds=((0,-np.inf),(np.inf,0)))
    #add a catch if it does not converge?
    
    #print the curve fit for reference
    print("Best a+exp(bt) fit:")
    print("\na: "+str(popt[0])+" +/- "+str(pcov[0,0]))
    print("b: "+str(popt[1])+" +/- "+str(pcov[1,1])+"\n")
    
    #catch if b param is non-negative as Poisson noise must be decaying
    if popt[1] > 0:
        print("\nPositive exponent parameter, decaying exponential fit failed.\n")
        return
        
    #apply correction and write output file
    corrected_powers = []
    corrections = []
    for i in range(len(freqs)):
            correction = func(freqs[i],popt[0],popt[1])#,popt[2]
            corrected_powers.append(powers[i]-correction)
            corrections.append(correction)
    
    #plot correction for debugging purposes        
    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots()
    #ax.plot(freqs, powers, color='red',zorder=1)
    #ax.plot(freqs, corrections, color='blue',zorder=2)
    #plt.xlabel("Frequency [mHz]")
    #plt.ylabel("Power")
    #plt.savefig("C:/Users/Val/Desktop/ls_poisson.png")
    #plt.gcf().clear()
            
    return corrected_powers
    
def ls_fft_compare(tf_path,output=-1,spacing=-1):
    """
    Creates a PSD-normalized Fourier transform and a PSD-normalized Lomb-Scargle transform of the data.
    Outputs a comparative graph. Assumes a uniform time spacing in the data set.
    Spacing default set to the average of the timestamp differences.
    """
    from astropy.stats import LombScargle
    import numpy as np
    import matplotlib.pyplot as plt
    
    #acquiring the data for the transforms
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    x = []
    y = []
    for line in contents:
        data = line.split("\t")
        x.append(float(data[3]))
        y.append(float(data[8]))
        
    #defining a uniform time spacing if no argument was passed
    if spacing == -1:
        t_dif = []
        for i in range(len(x)-1):
            t_dif.append(x[i+1]-x[i])
        spacing = sum(t_dif)/float(len(t_dif)) 
        print("Spacing is set to "+str(spacing)+" ms based off of average timestamp differences.")

    #the fourier transform of the data with PSD normalization    
    frequency = np.fft.fftfreq(len(x), spacing)
    y_fft = np.fft.fft(y)
    positive = (frequency > 0)
    frequency = frequency[positive]
    PSD_fourier = (1./len(x))*abs(y_fft[positive])**2
    
    #the lomb-scargle transform of the data with PSD normalization
    PSD_LS = LombScargle(x, y).power(frequency, normalization='psd')
    
    #plotting the comparison figure
    #if no output location given, store as home/ls_fft.png
    if output == -1:
            from os.path import expanduser
            output = expanduser("~")
            output = output+"ls_fft.png"
    plt.switch_backend('agg')
    fig, (ax0,ax1)  = plt.subplots(nrows=2,sharex=False)
    ax0.plot(frequency, PSD_fourier, color='#8f1402')
    ax0.set_title('Fourier (Top), Lomb-Scargle (Bottom)')
    ax0.set_ylabel('Power')
    ax1.plot(frequency, PSD_LS, color='#8f1402')
    ax1.set_xlabel('Frequency')
    ax1.set_ylabel('Power')
    plt.savefig(output,bbox_inches='tight')
    
def lightcurve(tf_path,output=-1,plotcolor='black',saveimage=True): 
    """
    Assumes data reduction is complete and that the textfile at tf_path is populated 
    as one row per reduced image and columns delimited by \t in order of:
    (0) stack, (1) exp avg (ms), (2) exp stdev (ms), (3) time avg (s), 
    (4) time stdev (s), (5) x_min, (6) y_min, (7) area, (8) count, (9) count error
    
    If saveimage is False, no .pngs are saved
    Returns 5 lists: time, time errors, photon counts, photon count errors, 
    normalized time and normalized photon counts
    """

    import matplotlib.pyplot as plt

    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    time = []
    time_err = []
    photon = []
    photon_err = []
    for line in contents:
        data = line.split("\t")
        #time axis
        time.append(float(data[3]))
        #time axis error
        time_err.append(float(data[4]))
        #photon axis
        photon.append(float(data[8]))
        #photon axis error
        photon_err.append(float(data[9]))

    #normalization of photon count with first count = 0
    #first_count = photon[0]
    first_count = 0
    photon_normed = [p - first_count for p in photon]
    photon_normed_kilo = [p/1000.0 for p in photon_normed]
    photon_err_kilo = [perr/1000.0 for perr in photon_err]
    
    #normalization of time with t1 = 0
    first_time = time[0]
    time_normed = [t - first_time for t in time]
    
    #creating the normalized lightcurve plot
    #if no output, store as home/lightcurve.png
    if (saveimage == True):
        if output == -1:
            from os.path import expanduser
            output = expanduser("~")
            output = output+"lightcurve.png"
        plt.switch_backend('agg')
        plt.rc('text',usetex=False)
        fig, ax0  = plt.subplots(figsize=(6,2),nrows=1,sharex=False)
        ax0.set_title('Light Curve: 28 September 2018')#HERE
        ax0.errorbar(time_normed,photon_normed_kilo,xerr=time_err,yerr=photon_err_kilo,fmt='.',
                     markerfacecolor=plotcolor,markeredgecolor=plotcolor)
        ax0.set_xlabel('Time (s)')
        ax0.set_ylabel('Photon Count $[10^3]$')
        plt.savefig(output,bbox_inches='tight')
        
    return time, time_err, photon, photon_err, time_normed, photon_normed
    
def QPO_detect(tf_path,output=-1,df_path=-1,minfreq=0,maxfreq=1,sinterms=1,
               saveimage=True,savetext=False,is_window=False,norm='standard',
               renorm=False,poisson=True,FALs=True,plotcolor='black',
               probs=[0.95,0.5,0.05]):
    """
    Assumes data reduction is complete and that the textfile output is populated 
    as one row per reduced image and columns delimited by \t in order of:
    (0) stack, (1) exp avg (ms), (2) exp stdev (ms), (3) time avg (s), 
    (4) time stdev (s), (5) x_min, (6) y_min, (7) area, (8) count, (9) count error
    
    minfreq and maxfreq are the set of frequencies to restrict the Lomb-Scargle
    (LS) autopower method; sinterms is the number of sinusoidal terms to be used in 
    the LS model fit 
    
    If saveimage is False (non-default), no .png is saved
    
    df_path is the data file which contains the reslults of the LS
    If savetext=False (default), this text file is not populated 
    If savetext=True, the text file is populated in one of two ways:
    * No frequency restriction -> frequency and power (2 columns)
    * Frequency restriction -> frequency, power, restricted frequency and 
    restructed power (4 columns) 
    
    If is_window is True (non-default), take the LS of the window function
    
    If renorm is True (non-default), every power p obtained from the LS is divided 
    by the maximum power before plotting
    
    If poisson is True (default), a decaying exponential is fit to the local minima 
    of the powers obtained from the LS, which is then used to correct all powers.
    
    If FALs is True (default), FALs are plotted in the background using the
    probabilities listed in probs parameter. probs is an optional argument.
    
    Output: 
    - Lomb-Scargle object and 4 lists: restricted frequency, error on restricted 
    frequency, restricted power, and heights for FALs
    ** All frequencies returned in Hz, NOT mHz 
    """
    
    from astropy.stats import LombScargle
    import matplotlib.pyplot as plt
    import numpy as np

    #Acquire the data for the Lomb-Scargle periodogram
    tf = open(tf_path,"r")
    contents = tf.readlines()
    tf.close()
    time = []
    time_err = []
    photon = []
    photon_err = []
    for line in contents:
        data = line.split("\t")
        time.append(float(data[3]))
        time_err.append(float(data[4]))
        photon.append(float(data[8]))
        photon_err.append(float(data[9]))

    # Normalization of photon count with first count = 0
    first_count = photon[0]
    photon_normed = [p - first_count for p in photon]
        
    # Normalization of time with t1 = 0
    first_time = time[0]
    time_normed = [t - first_time for t in time]
    
    # Edit the center data and fit mean depending on whether the data is to be viewed in window mode.
    cd_opt = True
    fm_opt = True
    if is_window is True:
        photon_normed = [1.0] * len(photon_normed)
        photon_err = None
        cd_opt = False
        fm_opt = False
    # If fit_mean is true, then (from LS documentation):
    # "Include a constant offset as part of the model at each frequency. 
    # This can lead to more accurate results, especially in the case of incomplete phase coverage."
    # If center_data is True, then (from LS documentation):
    # "Pre-center the data by subtracting the weighted mean of the input data. 
    # This is especially important if fit_mean = False"

    # Create the Lomb-Scargle object using autopower to generate a frequency sweep.
    limbo = LombScargle(time_normed, photon_normed, photon_err,center_data=cd_opt,
                        fit_mean=fm_opt,nterms=sinterms)    
    
    # Determine the error on the general frequency sweep to use in the restricted frequencies error
    frequency, power = limbo.autopower(normalization=norm)
    frequency_rebin = np.histogram(frequency, len(time))[1].tolist()
    while (len(time) < len(frequency_rebin)):
        frequency_rebin.pop()
    frequency_err = [0]
    for i in range(1, len(frequency_rebin)): # start at 1 to avoid t=0
        frequency_err.append(frequency_rebin[i]*time_err[i]/time[i])
    frequency_err[0] = frequency_err[1]
    frequency_err = np.histogram(frequency_err, len(frequency))[1].tolist()
    while (len(frequency) < len(frequency_err)):
        frequency_err.pop()

    # Set the frequencies and powers using autopower to generate a restricted frequency sweep.
    frequency_strict, power_strict = limbo.autopower(minimum_frequency=minfreq, 
                                                     maximum_frequency=maxfreq,
                                                     normalization=norm)
    
    # Determine the error on the resticted frequencies sweeped
    frequency_strict_rebin = np.histogram(frequency_strict, len(time))[1].tolist()
    while (len(time) < len(frequency_strict_rebin)):
        frequency_strict_rebin.pop()
    frequency_strict_err = [0]
    for i in range(1, len(frequency_strict_rebin)): # start at 1 to avoid t=0
        frequency_strict_err.append(frequency_strict_rebin[i]*time_err[i]/time[i])
    frequency_strict_err[0] = frequency_err[1]
    frequency_strict_err = np.histogram(frequency_strict_err, len(frequency_strict))[1].tolist()
    while (len(frequency_strict) < len(frequency_strict_err)):
        frequency_strict_err.pop()
    
    # Enforce a renormalization of min power = 0, max power = 1 if desired
    if renorm is True:
        m = max(power)
        power = [p/m for p in power]
        ms = max(power_strict)
        power_strict = [ps/ms for ps in power_strict]
        
    # Create the plot if image is to be saved, with mHz frequency units.
    # If FALs are to be added, calculate those and add to background. 
    # If no output location is given, store as home/lomb_scargle.png
    frequency_strict_milli = [f*1000.0 for f in frequency_strict]
    frequency_strict_milli_err = [f*1000.0 for f in frequency_strict_err]
    
    heights = limbo.false_alarm_level(probs) # moved outside conditional 
    
    if poisson is True:
        power_strict = correct_poisson(frequency_strict_milli,power_strict,fit_lows=True)
        if power_strict is None:
            print("\nPoisson correction failed, breaking QPO detection.\n")
            return
        print("\nNote: FALs not plotted as Poisson noise is corrected for.\n")
        FALs = False
    
    if (saveimage==True):
        if output == -1:
            from os.path import expanduser
            output = expanduser("~")
            output = output+"lomb_scargle.png"
        if (FALs==False):
            plt.switch_backend('agg')
            fig, ax0  = plt.subplots(figsize=(6,2),nrows=1,sharex=False)
            ax0.set_title('Lomb-Scargle Periodogram: 28 September 2018')#HERE
            ax0.set_ylabel('Power')   
            ax0.plot(frequency_strict_milli, power_strict, color=plotcolor) #mHz
            ax0.set_xlabel('Frequency [mHz]')
            plt.savefig(output,bbox_inches='tight')
        else:
            #heights = limbo.false_alarm_level(probs) #moved outside conditonal
            plt.switch_backend('agg')
            fig, ax0  = plt.subplots(figsize=(6,2),nrows=1,sharex=False)
            for i in range(len(heights)):
                ax0.axhline(y=heights[i],color='#d3d3d3')
            ax0.plot(frequency_strict_milli, power_strict, color=plotcolor) #mHz
            ax0.set_title('Lomb-Scargle Periodogram')
            ax0.set_ylabel('Power')
            ax0.set_xlabel('Frequency [mHz]')
            plt.savefig(output,bbox_inches='tight')
            
    # Create the text file if data is to be saved. 
    # If no output location is given, store as home/lomb_scargle.txt
    if (savetext == True):
        if df_path == -1:
            from os.path import expanduser
            df_path = expanduser("~")
            df_path = df_path+"lomb_scargle.txt"
        df = open(df_path,'w+')
        for i in range(len(frequency)): #Hz
            if (i < len(frequency_strict)):
                line = str(frequency[i])+"\t"+str(power[i])+"\t"+str(frequency_strict[i])+"\t"+str(power_strict[i])+"\n"
                df.write(line)
            else:
                line = str(frequency[i])+"\t"+str(power[i])+"\t \t \n"
                df.write(line) 
        df.close()
                        
    return limbo, frequency_strict, frequency_strict_err, power_strict, heights

def spinmob_fit(freq, freq_err, pows, num_peaks, width):
    """
    Input: the frequency, powers, and frequency errors output by QPO_detect, 
    as well as the expected number of peaks to be fit and the estimated peak 
    width (in mHz), assuming all peaks have about the same width. 
    
    Fit num_peaks Lorentzians to the periodogram. 
    
    *** Currently broken.  
    """
    import spinmob as s 
    f = s.data.fitter() # create fitter object
    
    # fitting function and parameters
    function = 'd +' # start with constant offset from power=0
    parameters = 'd' # fit parameters
    
    # add a lorentzian for each expected peak
    # lorentzian of form (g/2) / ((x-f0)**2 + (g/2)**2) 
    # where g is the FWHM and f0 is the centre
    for i in range(num_peaks):
        function += ' (g%d/2)/((x-f%d)**2 + (g%d/2)**2)'%(i,i,i)
        parameters += ',g%d,f%d'%(i,i)
    f.set_functions(f=function,p=parameters)
    
    print(function)
    print(parameters)
    
    # loading, setting data (where it currently breaks)
    freq_mhz = [f*1000.0 for f in freq]
    powers = pows.tolist()

    f.set_data(xdata=freq_mhz,ydata=powers)
    # says here that ndarray is not callable
    # but freqs_mhz and powers are NOT ndarrays, they are lists (verified)
       
    print("(Single) click the estimated baseline.")
    base_x, base_y = f.ginput()[0] # obtain the clicked (x,y)
    
    print("Next, (single) click the estimated centre of all 6 peaks, one at a time.")
    click_x0, click_y0 = f.ginput()[0] # obtain the clicked (x,y)
    click_x1, click_y1 = f.ginput()[0]
    click_x2, click_y2 = f.ginput()[0]
    click_x3, click_y3 = f.ginput()[0]
    click_x4, click_y4 = f.ginput()[0]
    click_x5, click_y5 = f.ginput()[0]
    
    # set guesses for fit parameters
    f.set(d=base_y, f0=click_x0, f1=click_x1, f2=click_x2, 
          f3=click_x3, f4=click_x4, f5=click_x5)
    f.set(g0=width, g1=width, g2=width, g3=width, g4=width, g5=width)
    f.set(xlabel='Frequency [Hz]', ylabel='Power')
    
    f.fit()
    
    print(f)

def lorentz(x, *p):
    """
    Input: an array of frequencies [Hz] and the parameters which describe the sum 
    of six Lorentzians. f_i represent centres of peaks, g_i their width, i_i
    their height.
    Output: the value of the six-Lorentzian sum.
    """
    import numpy as np 
    
    d, f0, f1, f2, f3, f4, f5, g0, g1, g2, g3, g4, g5, i0, i1, i2, i3, i4, i5 = p
    ret = d
    ret += (i0*(g0)**2)/((x-f0)**2 + (g0)**2) 
    ret += (i1*(g1)**2)/((x-f1)**2 + (g1)**2) 
    ret += (i2*(g2)**2)/((x-f2)**2 + (g2)**2) 
    ret += (i3*(g3)**2)/((x-f3)**2 + (g3)**2) 
    ret += (i4*(g4)**2)/((x-f4)**2 + (g4)**2) 
    ret += (i5*(g5)**2)/((x-f5)**2 + (g5)**2)
    return ret/(np.pi)

def six_lorentzian(freq, freq_err, pows, fit_params, FAL_heights, 
                   make_plot=True, bars=True, FALs=False):
    """
    Input: 
    - Frequencies, frequency errors [Hz] and powers of a periodogram
      to be fit;
    - Array of guesses for the (19) fit parameters in the order
      [d, f0, f1, ..., f5, g0, g1, ..., g5, i0, i1, ..., i5]
      where d is some baseline power, and f_i, g_i and i_i are the centres, widths
      and heights, respectively, of individual Lorentzians;
    - The false alarm probability (FALs) heights for the source being fit;
    - make_plot bool (default True): whether or not to plot & save figure;
    - bars bool (default True): whether or not to plot vertical bars at peaks
    - FALs bool (default False): whether or not to plot FALs lines
    
    Saves a figure showing the data to be fit and the fit function plotted over
    the entire range, and vertical bars showing where the peaks are located 
    Bars may be disabled. Fitting done using the least-squares method.
    
    Prints the peak locations, with errors, in a readable format. 
    
    TO IMPLEMENT: Errors on the input data. Very important. !!!!!!!!
    
    TO IMPLEMENT: Not just vertical bars, but bars of confidence. Maybe a 
    gradient or something else pretty.
    
    TO IMPLEMENT: A legend
    
    Output: Fit parameters; errors on the fit parameters; covariance matrix
    """
    
    from scipy import optimize as opt
    import numpy as np
    import matplotlib.pyplot as plt
    
    freq_mhz = [f*1000.0 for f in freq] # convert to mHz
    freq_err_mhz = [ferr*1000.0 for ferr in freq_err]
    powers = [p for p in pows]
 
    # use scipy to optimize
    # obtain parameters and covariance matrix
    popt, pcov = opt.curve_fit(lorentz, freq_mhz, powers, fit_params, maxfev=100000)
    # obtain errors on parameters
    perr = np.sqrt(np.diag(pcov)) 
    
    # inform user of peak locations
    print("\nPeaks located at:")
    for i in range(6):
        print("%f ± %f mHz"%(popt[i+1],perr[i+1]))
    #print("\nWith widths:")
    #for i in range(6):
    #    print("%f ± %f mHz"%(popt[i+7],perr[i+7]))
    
    if make_plot == True: 
        fig, ax = plt.subplots()
        fit_visual = lorentz(freq_mhz, *popt) # obtain fit over entire range
        
        ax.errorbar(freq_mhz,powers, color='cornflowerblue', zorder=1,
                    label="Data") # data
        ax.plot(freq_mhz, fit_visual, color='green',linestyle='--',zorder=2,
                label="Fit") # fit 
        
        if bars == True: # if bars desired at peak locations
            peak_locations = [popt[i] for i in range(1,7)]
            for p in peak_locations:
                ax.axvline(p, color='grey') # vertical lines at peak locations 

        if FALs == True: # if FALs desired 
            for h in FAL_heights:
                ax.axhline(h,color='#d3d3d3')
                            
        plt.title("6-Lorentzian Fit")
        plt.xlabel("Frequency [mHz]")
        plt.ylabel("Power")
        plt.legend()
        plt.grid(which="both")
        #plt.show()
        plt.savefig("data+fit+peaks.pdf")
        plt.gcf().clear()
        
    
    return popt, perr, pcov

def source_versus_alt(freq1, freqerr1, power1, freq2, freqerr2, power2,
                      fit_params, FAL_heights, output, 
                      plotcolor, compare_method="vertical lines", FALs=False):
    """
    Input:  
    - Frequency, frequency errors [Hz] and powers of one periodogram 
      and those of another periodogram to which we compare;
    - Parameters to be used as guesses to the fit of the second periodogram [see
      six_lorentzian for formatting]; 
    - Heights for the false alarm probabilities (FALs) of the second periodogram, obtained 
      from the QPO_detect() method; 
    - Comparison method (vertical lines for the peaks of the second periodogram, by default); 
    - Name of the plot to be saved;
    - Color for the plot of the first periodogram
    
    Extracts 6 peaks from the second periodogram by fitting a 6-Lorentzian 
    curve. Plots these peaks as vertical lines by default, and simply plots
    the fit of the second periodogram over the entire range otherwise
    
    Usage:
    Allows the user to see if peaks present in a source which is not of interest
    are present in the periodogram of a source of interest. e.g., are peaks in 
    the periodogram of MAXI J1820+070 (input as source #2) unique, or are these same peaks 
    located in periodograms of other sources (input as source #1) in the field? 
    Do sources in the field have peaks that MAXI doesn't?
    
    TO IMPLEMENT: Not just vertical bars, but bars of confidence. Maybe a 
    gradient or something else pretty.
    
    TO IMPLEMENT: A legend.
    
    Output: None
    """
    
    import matplotlib.pyplot as plt
    
    # Hz -> mHz conversion for plotting first/second periodogram
    freq2_mhz = [f*1000.0 for f in freq2]
    freq1_mhz = [f*1000.0 for f in freq1]
    freqerr1_mhz = [ferr*1000.0 for ferr in freqerr1]
    power1 = [p for p in power1]
    
    # fit second periodogram
    fitopt, fiterr, fitcov = six_lorentzian(freq2, 
                            freqerr2, power2, fit_params, FAL_heights,
                            make_plot=False, FALs=True) # mHz
    
    # obtain peak locations 
    peak_locations = [fitopt[i] for i in range(1,7)] # mHz
    
    fig, ax = plt.subplots()
    # plot first periodogram
    ax.plot(freq1_mhz, power1, color=plotcolor, 
            label="Comparison source #3", zorder=3) #mHz
    
    if compare_method == "vertical lines":
        # plot vertical lines representing peaks from alternate source
        for p in peak_locations:
            ax.axvline(p, color='blue', label="MAXI J1820+070 peak", zorder=2)
    else:
        # plot the fit obtained for the second periodogram
        fit_visual = lorentz(freq2_mhz, *fitopt)
        ax.plot(freq2_mhz, fit_visual, color='#1fa774',linestyle='-.', 
                 label="MAXI J1820+070 fit",zorder=4,linewidth=1.5) # fit 
    
    if FALs == True:
        for h in FAL_heights:
                plt.axhline(h,color='#d3d3d3', zorder=1)
    
    #plt.title("Two-source Comparison") # more general 
    plt.title("MAXI J1820+070 versus alternate source #3: 28 September 2018") # for us
    plt.xlabel("Frequency [mHz]")
    plt.ylabel("Power")
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[-2:], labels[-2:])
    #plt.show()
    plt.savefig(output)
    plt.gcf().clear()
    
    
    
    
