# Shuffles and splits the dataset into a training and testset.
# In the dataset we have cellnames, spiketimes and awake_angle.
train_test_split = function(data){
    
    # Create a vector of timestamps
    timespan = length(data$awake_angle)
    time_seq = seq(1,timespan)
    
    # Shuffle the timesequence
    # shuffled_time_seq = sample(time_seq, timespan)
    # Dont shuffle the timesequence
    shuffled_time_seq = time_seq
    
    # Split and order the timesequence
    training_times = order(shuffled_time_seq[1:timespan/2])
    test_times = order(shuffled_time_seq[-(1:timespan/2)])
    
    training_data = data
    test_data = data
    
    # Extract relevant awake angles
    training_data$awake_angle = data$awake_angle[training_times]
    test_data$awake_angle = data$awake_angle[test_times]
    
    for(cellname in data$cellnames){
        #Extract relevant spiketimes
        training_data$spiketimes[[cellname]] = data$spiketimes[[cellname]][data$spiketimes[[cellname]] %in% training_times]
        test_data$spiketimes[[cellname]] = data$spiketimes[[cellname]][data$spiketimes[[cellname]] %in% test_times]
    }
    return(list(training_data = training_data,
                test_data = test_data))
}

# Removes indexs and values where we dont have angle data 
processMomentaryPopulationVector = function(momentary_population_matrix, smoothing_window, total_time){
    
    # Convert to firing rate
    momentary_population_matrix = momentary_population_matrix/(2*smoothing_window/1000)
    
    # Remove bins that are smaller than 200ms or where we don't have the relevant
    # angle data for the mouse
    usable_indexs = smoothing_window:(total_time-smoothing_window)
    usable_indexs = usable_indexs[which(!is.nan(test_data$awake_angle[usable_indexs]))]
    
    momentary_population_matrix = as.matrix(momentary_population_matrix[,usable_indexs])
    
    return(list(usable_indexs=usable_indexs,
                momentary_population_matrix=momentary_population_matrix))
}

# Calculates the momentary population vector for every timestamp
calculateMomentaryPopulationVectors = function(data, smoothing_windows){
    
    momentary_population_matrix_list = list()
    
    num_cells = length(data$cellnames)
    total_time = length(data$awake_angle)
    
    # Matrix with cells as rows and timestamps as columns
    activity_matrix=matrix(0, nrow=num_cells, ncol=total_time)
    
    # Set activity to one for every timestamp there is a spike in the
    # corresponding cell.
    for(i in 1:num_cells){
       activity_matrix[i,data$spiketimes[[data$cellnames[[i]]]]] = 1
    }
    
    activity_matrix = Matrix(activity_matrix, sparse=TRUE)
    
    # Create matrix where we will apply the smoothing window
    momentary_population_matrix = activity_matrix
    
    left = activity_matrix
    right = activity_matrix
    zero_vec = Matrix(0, nrow=num_cells, ncol=1, sparse=TRUE)
    
    for(i in 1:max(smoothing_windows)){
        left = cbind(zero_vec,left[,-dim(left)[2]])
        right = cbind(right[,-1], zero_vec)
        momentary_population_matrix = momentary_population_matrix + left + right
        if(i %in% smoothing_windows){
            temp_result = processMomentaryPopulationVector(momentary_population_matrix, i, total_time)
            momentary_population_matrix_list = cbind(momentary_population_matrix_list, list(temp_result))
        }
    }
    names(momentary_population_matrix_list) = smoothing_windows
    colnames(momentary_population_matrix_list) = smoothing_windows
    return(momentary_population_matrix_list)
}

# Returns the reference population vector based on an input angle
calculateReferencePopulationVector = function(angle, tuning_curves){
    
    #Place the input angle into the corresponding bin
    bin_length = 2*min(tuning_curves$angle_bins)
    bin_index = ceiling(angle/bin_length)
    binned_input_angle = tuning_curves$angle_bin[bin_index]
    
    #Fetch the reference population vector for the binned angle
    reference_vector = tuning_curves[tuning_curves$angle_bin == binned_input_angle,]$firing_rate
    return(reference_vector)
}

# Finds the angle with the lowest pearson coefficient
findHighestCorAngle = function(momentary_population, tuning_curves){

    population_vectors = lapply(unique(tuning_curve$angle_bins), calculateReferencePopulationVector, tuning_curves=tuning_curves)
    correlations = lapply(population_vectors, cor, momentary_population, use="pairwise.complete.obs", method="pearson")

    best_angle=unique(tuning_curves$angle_bins)[which.max(correlations)]
    
    return(best_angle)
}

angleDiff = function(angle1, angle2){
    angle1 = pi/180*angle1
    angle2 = pi/180*angle2
    diff = abs(atan2(sin(angle1-angle2), cos(angle1-angle2)))
    return(diff*180/pi)
}


## Mutual information

calculateMutualInformation = function(binned_awake_angles, momentary_population_matrix, cellnumber, smoothing_window){
    
    library(entropy)
    
    response = momentary_population_matrix[cellnumber,]
    
    # Convert from firing rate to number of spikes
    response = response*(2*smoothing_window/1000)
    

    stimulus_dist = estimateDistribution(binned_awake_angles)
    p_conditional_response = calculateConditionalResponse(response, binned_awake_angles)
    response_dist = estimateDistribution(response)

    ## I tried implementing this myself first, the trends look correct, but my values are wrong
    ## probably some normalization issues.
    # information = 0
    # for(s in 1:length(stimulus_dist)){
    #     for(r in 1:length(response_dist)){
    #         log_term = log2(p_conditional_response[r,s]/response_dist[r])
    #         if(is.infinite(log_term)){
    #             next
    #         }
    #         information = information + stimulus_dist[s]*p_conditional_response[r,s]*log_term
    #         # print(paste('i', s, 'j', r))
    #         # print(information)
    # 
    #     }
    # }
    
    # count_stim = freqs.empirical(table(binned_awake_angles))
    # count_resp = freqs.empirical(table(response))
    return(mi.empirical(p_conditional_response))
    # return(information)
}

# I.e P(r|s)
calculateConditionalResponse = function(response, binned_awake_angles){

    # Stimulus response dataframe
    df = data.frame(stimulus=binned_awake_angles, response=as.factor(response))
    result = rep(0,length(unique(response)))
    for(stim in unique(df$stimulus)){
        filtered = df %>% filter(stimulus==stim) %>% select(response)
        response_given_stim = as.data.frame(table(filtered))
        result = cbind(result, response_given_stim$Freq)
    }
    result = result[,-1]
    result = result/sum(rowSums(result))
    return(result)
}

# Plug in principle
estimateDistribution = function(data){
    return(as.data.frame(table(data))$Freq/length(data))
}
