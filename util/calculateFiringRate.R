#Creates angle bins given the number of bins wanted
getAngleBins = function(angle_data, num_bins){
    angle_data = unlist(angle_data)
    angle_bins=rep(0,num_bins)
    for(i in 1:num_bins){
        lower_bound = pi*2/num_bins*(i-1)
        upper_bound = pi*2/num_bins*i
        angle_bins[i] = 0.5*(lower_bound+upper_bound)*360/(2*pi)
    }
    
    return(angle_bins)
}

calculateBinnedFiringRate = function(spiketimes, angle_data, num_bins=40){
  
    angle_data = unlist(angle_data)
    spiketimes = unlist(spiketimes)

    ## Find Which Angles Correspond to Firing Times
    # Make a row vector with zeros to contain the angles corresponding to the time of firing.
    index_angle = spiketimes
    index_angle[index_angle==0]=1
    angles_at_firing = angle_data[index_angle]
    ## Use angles_at_firing to Make Firing Rates For Given Angle Bins
    # Now that the angles at which a given cell is firing have been computed,
    # angles_at_firing, they can be used to compute the firing rates necessary to plot the tuning cuves.
    
    # Make a row vector firing_rate of zeros with the total number of bins (nrpoints)
    firing_rate=rep(0, num_bins)
    
    # BINNING to find how many spikes are within each radian interval.
    
    # For each bin:
    # The firing rate for the particular bin is
    # computed by dividing the number of spikes in the interval by the
    # occupancy. Recall that each element in occupancy is the total count
    # of milliseconds the animal spends in a particular bin. Multiply by
    # 1,000 because occupancy is in milliseconds.
    
    # Occupancy is  defined for each bin as the number of elements in
    # binnedAwakeangle_data between two given angles AA and BB. The values of
    # AA and BB is given by dividing the complete circle (2*pi) into the
    # number of bins (nrpoints) and then multiply by the previous bin number
    # or the current bin number (to find an angle interval).
    
    # As such, in occupancy you get a "measure" of how much time was spent in
    # each interval it is the total count of milliseconds the animal spends
    # in a particular bin. This follows as binnedAwakeangle_data expand the
    # awakeangle_data vector to contain a measure for each millisecond. Thus, if
    # an animal spends several milliseconds in a particular angle (e.g. when
    # staring at a cue), the subsequent elements in the vector will be equal.
    
    for(i in 1:num_bins){
        lower_bound = pi*2/num_bins*(i-1)
        upper_bound = pi*2/num_bins*i
        numspikes=length(which(lower_bound<angles_at_firing & angles_at_firing<upper_bound))
        occupancy=length(which(lower_bound<angle_data & angle_data<upper_bound))
        firing_rate[i]=1000* numspikes/occupancy# Convert to firing rates
    }
    return(list(firing_rate=firing_rate,
                angles_at_firing=angles_at_firing))
}