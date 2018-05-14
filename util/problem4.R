#Calculates the activity over a timeframe
calculateActivity = function(data, cellname, timeframe=200){

    spiketimes = data$spiketimes[[cellname]]
    total_time = length(data$awake_angle)
    
    activity = rep(0, total_time)
    sum_activity = activity
    
    activity[spiketimes] = 1
    
    for(i in 1:timeframe){
        sum_activity = sum_activity + c(rep(0,i),activity[1:(total_time-i)])
    }

    # Remove activity that doesnt have enough history to include timeframe ms of data 
    result = sum_activity[(timeframe+1):total_time]
    return(result)
}

#Bins the awake_angles vector
binAwakeAngles = function(awake_angle, bins=40){
    angle_bins = getAngleBins(40)
    bin_length = angle_bins[1]*2
    which_bin = ceiling(awake_angle*360/(2*pi)/bin_length)
    which_bin[which_bin==0] = 1
    binned_angles = angle_bins[which_bin]
    return(binned_angles)
}

# Set the threshold of activity at 1/3 of the max activity
thresholdedActivity = function(binned_angles, activity, cellname, bins=40){

    # activity = activity[which(!is.na(binned_angles))]
    # binned_angles = na.omit(binned_angles)
    # 
    df=data.frame(binned_angles=binned_angles, activity=activity, cellname=cellname)
    df=df[order(df$activity),]
    
    end = length(df$activity)

    df$category                = factor('medium', levels=c('low', 'medium', 'high'), ordered=TRUE)
    df$category[1:(end/3)]     = factor('low',  levels=c('low', 'medium', 'high'), ordered=TRUE)
    df$category[(2/3*end):end] = factor('high', levels=c('low', 'medium', 'high'), ordered=TRUE)
    
    return(df)
}

# Calculate the firing rate for each angle bin for the different thresholds of activity
thresholdedFiringRate = function(activity_df, data, cellname, bins=40){
    
    spiketimes = data$spiketimes[[cellname]]
    angle_data = data$awake_angle
    categories = as.ordered(c('low', 'medium', 'high'))

    angle_bins = getAngleBins(bins)
    firing_rate_df = data.frame()
    
    for(category in categories){
        
        spike_index = spiketimes
        spike_index[spike_index==0]=1
        spike_index = spike_index[which(activity_df$category[spike_index]==category)] # Filter angles in the correct category
        angles_at_firing = angle_data[spike_index]
        
        firing_rate = rep(0, bins)
        df = data.frame()
        
        for(i in 1:bins){
            lower_bound = pi*2/bins*(i-1)
            upper_bound = pi*2/bins*i
            numspikes=length(which(lower_bound<angles_at_firing & angles_at_firing<upper_bound))
            if(numspikes < 1){
                firing_rate[i]=NA
            }
            else{
                occupancy=length(which(lower_bound<angle_data & angle_data<upper_bound))
                firing_rate[i]=1000*numspikes/occupancy # Convert to firing rates
            }
        }
        
        df = data.frame(firing_rate=firing_rate, angle_bin=angle_bins, category = category, cellname=cellname)
        firing_rate_df = rbind(firing_rate_df, df)
    }
    return(firing_rate_df)
}

# Calculate the activity and calculate the firing_rate at different angles for thresholds based on the activity
calculateTresholdedValues = function(data, cell){
    activity = calculateActivity(data, cell)
    binned_awake_angles = binAwakeAngles(data$awake_angle[201:length(data$awake_angle)])
    activity_df = thresholdedActivity(binned_awake_angles, activity, cell)
    firing_rate_df = thresholdedFiringRate(activity_df, data, cell)
    return(list(activity_df=activity_df, firing_rate_df=firing_rate_df))
}