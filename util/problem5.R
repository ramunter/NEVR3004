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

# Calculates the momentary population vector for every timestamp
calculateMomentaryPopulationVector = function(data, smoothing_window=100){
    num_cells = length(data$cellnames)
    total_time = length(data$awake_angle)
    
    # Matrix with cells as rows and timestamps as columns
    activity_matrix=matrix(0, nrow=num_cells, ncol=total_time)
    
    # Set activity to one for every timestamp there is a spike in the
    # corresponding cell.
    for(i in 1:num_cells){
       activity_matrix[i,data$spiketimes[[data$cellnames[[i]]]]] = 1
    }

    # Create matrix where we will apply the smoothing window
    momentary_population_matrix = activity_matrix

    for(i in 1:smoothing_window){
        momentary_population_matrix = momentary_population_matrix + 
            cbind(matrix(0, nrow=num_cells, ncol=i), activity_matrix[,1:(total_time-i)]) + 
            cbind(activity_matrix[,(i+1):total_time], matrix(0, nrow=num_cells, ncol=i))
    }
    momentary_population_matrix = momentary_population_matrix/(2*smoothing_window/1000)
    
    # Remove bins that are smaller than 200ms or where we don't have the relevant
    # angle data for the mouse
    
    usable_indexs = smoothing_window:(total_time-smoothing_window)
    usable_indexs = usable_indexs[which(!is.nan(test_data$awake_angle[usable_indexs]))]
    
    momentary_population_matrix = momentary_population_matrix[,usable_indexs]

    return(list(momentary_population_matrix = momentary_population_matrix,
                usable_indexs = usable_indexs))
}

# Finds the angle with the lowest pearson coefficient
findHighestCorAngle = function(momentary_population, tuning_curves){
    best_correlation = -1
    for(angle in unique(tuning_curves$angle_bins)){
        reference_vector = calculateReferencePopulationVector(angle, tuning_curves)
        correlation = cor(momentary_population, reference_vector, use="pairwise.complete.obs", method="pearson")
        if(correlation > best_correlation){
            best_correlation = correlation
            best_angle = angle
        }
    }
    return(best_angle)
}

angleDiff = function(angle1, angle2){
    angle1 = pi/180*angle1
    angle2 = pi/180*angle2
    diff = abs(atan2(sin(angle1-angle2), cos(angle1-angle2)))
    return(diff*180/pi)
}
