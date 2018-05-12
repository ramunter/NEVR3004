calculateCoefficentOfVariation = function(cellname, data, tuning_curve){
    
    source('./util/calculateFiringRate.R')
    spiketimes = data$spiketimes[[cellname]]
    angles_at_firing = calculateBinnedFiringRate(spiketimes, data$awake_angle)$angles_at_firing
    
    # Find Inter Spike Intervals and their corresponding angle
    ISI = rep(0, length(spiketimes)-1)
    angle_ISI = rep(0, length(spiketimes)-1)
    
    for(i in 1:(length(spiketimes)-1)){
        ISI[i] = spiketimes[i+1]-spiketimes[i]
        angle_ISI[i] = angles_at_firing[i]
        
        if (ISI[i]>250){ # Remove ISIs above 250 ms
            ISI[i] = NaN
            angle_ISI[i] = NaN
        }
    }
    
    # Collect data in data frame
    df_isi <- data.frame(isi = ISI, rad = angle_ISI, deg = angle_ISI*360/(2*pi), cellname=cellname)
    
    bin_width = tuning_curve$angle_bins[1]*2 # Width of each bin
    bin_index = ceiling(df_isi$deg/bin_width) # Says which bin each angle belongs to
    binned_isi_count = rep(0, length(tuning_curve$angle_bins))
    
    for(index in bin_index){
        if(!is.nan(index)){ 
            binned_isi_count[index] = binned_isi_count[index] + 1 # Counts number of ISIs in each bin
        }
    }
    
    # Append data to data frame
    df_isi = cbind(df_isi, index=bin_index, binned_angle = tuning_curve$angle_bins[bin_index])
    df_isi = df_isi[!is.nan(df_isi$deg),]
    
    
    CV = rep(NaN, length(tuning_curve$angle_bins))

    
    # Compute coefficient of variation
    for(j in 1:length(tuning_curve$angle_bins)){
        if(binned_isi_count[j] > 50){
            binned_isi = na.omit(df_isi[df_isi$index == j,])$isi
            CV[j] = sd(binned_isi)/mean(binned_isi) # Coefficient of variation
        }
    }
    
    # CV = 1 means totally random
    return(list(CV=CV,
                isi=df_isi))
    
} # function