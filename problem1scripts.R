# Creates a dataframe with the cell firing rates at different angles.
createTuningCurveDataFrame = function(data){
    source('./Util/calculateFiringRate.R')
    tuning_curves = data.frame()
    angle_bins = getAngleBins(data$awake_angle, 40)
    
    for(cellname in data$cellnames){
        if(length(data$spiketimes[[cellname]]) != 0){
            firing_rate = calculateBinnedFiringRate(data$spiketimes[cellname], data$awake_angle, angle_bins)$firing_rate
            name = unlist(strsplit(cellname, 'C'))
            result = data.frame(angle_bins, firing_rate, cellname, tetrode=name[1], cell=name[2])
            tuning_curves = rbind(tuning_curves, result)
        }
    }
    return(tuning_curves)
}