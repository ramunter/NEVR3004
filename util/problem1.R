# Creates a dataframe with the cell firing rates at different angles.
createTuningCurveDataFrame = function(data, bins=40, min_firing_rate=10){
    source('./util/calculateFiringRate.R')
    tuning_curves = data.frame()
    angle_bins = getAngleBins(bins)
    
    # Loop through cells
    for(cellname in unique(data$cellnames)){
        
        # If we have measurements
        if(length(data$spiketimes[[cellname]]) != 0){
            firing_rate = calculateBinnedFiringRate(data$spiketimes[cellname], data$awake_angle)$firing_rate
            
            # If the minimum firing rate is high, the cell is not tuned
            if(min(firing_rate)<min_firing_rate){
                name = unlist(strsplit(cellname, 'C'))
                result = data.frame(angle_bins, firing_rate, cellname, tetrode=name[1], cell=name[2])
                tuning_curves = rbind(tuning_curves, result) 
            }
        }
    }
    return(tuning_curves)
}