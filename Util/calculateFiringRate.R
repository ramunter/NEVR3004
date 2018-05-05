calculateFiringRate = function(spiketimes, angle_data){
  
    angle_data = unlist(angle_data)
    spiketimes = unlist(spiketimes)
  
    ## Define Occupancy
    
    # Occupancy is needed to compute the tuning curve and is defined for each
    # bin as the number of elements in binnedAwakeangle_data between two given
    # angles AA and BB. The values of AA and BB is given by dividing the
    # complete circle (2*pi) into the number of bins (nrpoints) and then
    # multiply by the previous bin number or the current bin number (to
    # find an angle interval).
    
    # As such, in occupancy you get a "measure" of how much time was spent in
    # each interval it is the total count of milliseconds the animal spends
    # in a particular bin. This follows as binnedAwakeangle_data expand the
    # awakeangle_data vector to contain a measure for each millisecond. Thus, if
    # an animal spends several milliseconds in a particular angle (e.g. when
    # staring at a cue), the subsequent elements in the vector will be equal.
    
    
    nrpoints=40
    occupancy=rep(0,nrpoints)
    plotting_angles=rep(0,nrpoints)
    for(i in 1:nrpoints){
        AA = pi*2/nrpoints*(i-1)
        BB = pi*2/nrpoints*i
        plotting_angles[i] = 0.5*(AA+BB)*360/(2*pi)
        occupancy[i]=length(which(AA<angle_data & angle_data<BB))
    }
    
    ## Find Which Angles Correspond to Firing Times
    
    # Make a row vector with zeros to contain the angles corresponding to the time of firing.
    angles_at_firing=rep(0, length(spiketimes))
    ####find the angle that corresponds to each spike time.
    for(i in 1:length(spiketimes)){
        if(spiketimes[i] == 0){
            index_angle=1
        }
        else{
            index_angle=spiketimes[i]
        }
    
        angles_at_firing[i]=angle_data[index_angle]
    }
    ## Use angles_at_firing to Make Firing Rates For Given Angle Bins
    # Now that the angles at which a given cell is firing have been computed,
    # angles_at_firing, they can be used to compute the firing rates necessary to plot the tuning cuves.
    
    # Make a row vector firing_rate of zeros with the total number of bins (nrpoints)
    firing_rate=rep(0,nrpoints)
    
    # BINNING to find how many spikes are within each radian interval.
    
    # For each bin:
    # AA will contain the lower limit of the given bin interval, BB will
    # contain the upper limit of the same bin. The limits are given in
    # radians (because we are to bin the angles at which the cell fired).
    # The number of spikes, numspikes, contained in the given interval is
    # found. Then the firing rate for the particular bin is
    # computed, dividing the number of spikes in the interval by the
    # occupancy. Recall that each element in occupancy is the total count
    # of milliseconds the animal sp}s in a particular bin. Multiply by
    # 1,000 because occupancy is in milliseconds. Hence, the unit of the
    # tuning curve's firing rate becomes (spikes * 1000) / (1000 * s) = spikes/s.
    
    for(i in 1:nrpoints){
        AA = pi*2/nrpoints*(i-1)
        BB = pi*2/nrpoints*i
        numspikes=length(which(AA<angles_at_firing & angles_at_firing<BB))
        firing_rate[i]=1000* numspikes/occupancy[i]# Convert to firing rates
    }
    return(list(firing_rate=firing_rate,
                plotting_angles=plotting_angles,
                angles_at_firing=angles_at_firing))
}