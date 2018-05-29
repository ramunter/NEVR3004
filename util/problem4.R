calculateActivity = function(data, cellname, timeframe=200){
    library(Matrix)
    
    spiketimes = data$spiketimes[[cellname]]
    total_time = length(data$awake_angle)
    
    activity=rep(0, total_time)
    activity[spiketimes] = 1
    
    # Create matrix where we will apply the smoothing window
    sum_activity = activity
    
    shifted = activity

    for(i in 1:timeframe){
        shifted = lag(shifted)
        sum_activity = shifted + sum_activity
    }
    # sum_activity = sum_activity[!is.na(sum_activity)]
    return(sum_activity)
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
thresholdedActivity = function(data, binned_angles, activity, cellname, bins=40){

    # activity = activity[which(!is.na(binned_angles))]
    # binned_angles = na.omit(binned_angles)
    usable_indexs = which(!is.na(activity) & activity!=0)
    
    binned_angles = binned_angles[usable_indexs]
    activity = activity[usable_indexs]

    df=data.frame(index=usable_indexs, binned_angles=binned_angles, activity=activity,
                  cellname=cellname)
    
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
    categories = factor(c('low', 'medium', 'high'), levels=c('low', 'medium', 'high'), ordered=TRUE)

    angle_bins = getAngleBins(bins)
    firing_rate_df = data.frame()
    
    for(j in 1:3){
        
        spiketimes_filtered = activity_df %>% filter(category==categories[j], index %in% spiketimes) %>% select(index)
        angles_at_firing = angle_data[unlist(spiketimes_filtered)]

        bin_size = 2*pi/bins
        
        angle_bins = seq(bin_size/2,
                         2*pi-bin_size/2, bin_size)
        angle_bins = 360*angle_bins/(2*pi)
        
        firing_rate = rep(0, bins)
        
        for(i in 1:bins){
            lower_bound = bin_size*(i-1)
            upper_bound = bin_size*i
            
            numspikes=length(which(lower_bound<angles_at_firing & angles_at_firing<upper_bound))
            occupancy=length(which(lower_bound<angle_data & angle_data<upper_bound))
            firing_rate[i]=1000*numspikes/occupancy # Convert to firing rates
        }
        
        df = data.frame(firing_rate=firing_rate, angle_bin=angle_bins, category = categories[j], cellname=cellname)
        firing_rate_df = rbind(firing_rate_df, df)
    }
    return(firing_rate_df)
}

# Calculate the activity and calculate the firing_rate at different angles for thresholds based on the activity
calculateTresholdedValues = function(data, cell){
    activity = calculateActivity(data, cell)
    binned_awake_angles = binAwakeAngles(data$awake_angle)
    activity_df = thresholdedActivity(data, binned_awake_angles, activity, cell)
    firing_rate_df = thresholdedFiringRate(activity_df, data, cell)
    
    return(list(activity_df=activity_df, firing_rate_df=firing_rate_df))
}