interpolateBetweenPoints = function(xx,ind0,ind1,y0,y1,angle_freq){
  # if this looks stange, google "linear interpolation" in wikipedia
  yd=y1-y0
  if(is.nan(y0) || is.nan(y1)){   #if firingTime is close to a nan angle, we do not interpolate 
    return(NaN)
  }
  if(yd>pi){               
    y1=y1-2*pi
  }
  else if(yd<(-pi)){
    y1=y1+2*pi
  }
  x0=1000*ind0/angle_freq            
  x1=1000*ind1/angle_freq
  yy=y0 + (xx-x0)*(y1-y0)/(x1-x0) #linear interpolation
  return(yy%%(2*pi))
}

importData = function(directory_name, numfiles){
    ## Importing raw data
     
    # Create empty cell array to store res and clu data from each tetrode
    # tetrode is a multiple electron measurement device.
    resclu = matrix(list(), numfiles, 2) #
    # res data stores the timestamps corresponding to when a cell fired.
    # clu data stores the cell numbers corresponding to each timestamp.
    
    for(i in 1:numfiles){
        # Set file names:
        clufile = paste(directory_name, 'clu', i, sep=".")
        resfile = paste(directory_name, 'res', i, sep=".")
        clusters = as.matrix(read.table(clufile))
        timestamps = as.matrix(read.table(resfile))
        
        # Assign values to the previously defined cell array resclu. Cluster
        # arrays are assigned to the ith row and 1st column and timestamps
        # are assigned to the ith row and second column. Curly brackets are
        # used because resclu is a cell array.
        resclu[[i,1]] = clusters
        resclu[[i,2]] = timestamps
    }
    
    
    ## Cell Names and Timestamps
    
    awake = as.matrix(read.table(paste(directory_name,'states','Wake', sep='.')))
    awake_start = awake[1,1]
    awake_end = awake[dim(awake)[1],dim(awake)[2]]
    
    electrode_freq=20000;  #Electrophysiologal rate (i.e. the rate at which the raw data were recorded; Hz)
    
    ### Create empty cell array to store spike times.
    spiketimes=list();
    for(tetrode in 1:dim(resclu)[1]){ ### Length of resclu variable will be the same as numfiles
        
        timestamps_ms=1000*resclu[[tetrode,2]]/electrode_freq # Time of firings in all cells measured by electrode. 
        #Dividing timestamps by electrode freq rescales into seconds, multiplying by 1000 then rescales to milliseconds.
        tetrode_number=resclu[[tetrode,1]]
        tetrode_number = tetrode_number + 1 # The first cell in the tetrode is marked as 0, but matlab cannot handle an index of 0, so we add 1 to all cell numbers.
        number_cells_in_tetrode=tetrode_number[1,]-1 #The first element in the cluster is an overview of nr of Cells measured. We minus 1 because we previously added one and to get the true number of cells we must minus that again.
        cellnames_per_tetrode = matrix(list(), 1, number_cells_in_tetrode)
        spiketimes_this_tetrode = matrix(list(), 1, number_cells_in_tetrode)
        for(cellnum in 1:number_cells_in_tetrode){
            cellnames_per_tetrode[,cellnum] = paste('T',tetrode,'C',cellnum,sep='')
            cell_index = which(tetrode_number == cellnum)-1 #### Must be -1 because first index of Tetrode Number does not correspond to a real cell number
            timestamps_cell=timestamps_ms[cell_index]
            #goes and finds the cell timestamp at the input index, as decided by the previous code
            #The -1 makes sure it doesn't find index one because that is not the name of a cell.
            #Since it's a for loop, this will create an array for each cell
            #individually (and iterates for each cell in each cluster because
            #it's a for loop).
            index_awake_timestamps_cell = which((timestamps_cell>(awake_start*1000) & (timestamps_cell<(awake_end*1000)))) #find returns the Index of all times when the mouse is awake and therefore has usable firing data.
            awake_cell_timestamps = timestamps_cell[index_awake_timestamps_cell] #now only firings when the mouse is awake.
            awake_cell_timestamps = floor(awake_cell_timestamps-awake_start*1000) #### Rescale awake data to starting time of 0 with whole integers.
            spiketimes_this_tetrode[[1,cellnum]]= awake_cell_timestamps
        }
        if(tetrode == 1){
            cellnames=cellnames_per_tetrode
        }
        else{
            cellnames=cbind(cellnames, cellnames_per_tetrode)
        }
        
        spiketimes = cbind(spiketimes, spiketimes_this_tetrode)
        spiketimes_this_tetrode=NULL; #variable must be cleared on every iteration of tetrode for loop.
    }
    
    names(spiketimes) = cellnames[1,]
    ####
    ## Angle Data
    angle_freq=1250/32;    # Angle measurement rate (1,250 is the EEG sampling rate and 32 the
                           #number of LFP samples; Hz/sample)* 
    
    angle_data=paste(directory_name, '-HD-filt.ang', sep="")
    angle_data=as.matrix(read.table(angle_data))
    
    ### -1 was used as a code for missing data. In order for Matlab to
    ### understand that these values are missing, we must recode them to NaN
    ### (Not a Number)
    tmpIndex1=which(angle_data==-1)
    angle_data[tmpIndex1,]=NaN
    
    # Find the index for a given time to find the corresponding angle in
    # angle_data. round() rounds to the nearest integer (away from zero for
    # > 0.5). Units: s * (Hz/sample) = 1/sample (i.e. index of a given sample).
    # Then, use the indices to select the angle_data recorded from the awake
    # mouse.
    indexAwakeStart=round(awake_start*angle_freq)
    indexAwakeEnd= round(awake_end*angle_freq)
    awake_angle_data=angle_data[indexAwakeStart:indexAwakeEnd,]
    
    ## BIN THE ANGLE DATA
    # Creating binnedangle_data, a vector containing one angle measurement
    # for each awake ms, because the length of the binnedangle_data is
    # given by the number of milliseconds between awakeStart and awakeEnd.
    # ceil() rounds a decimal number upwards, regardless of the value of the decimal.
    
    # For each index i of binnedangle_data, find the corresponding index in
    # the awake_angle_data vector by dividing the index by angle_freq in mHz/sample.
    # Then, find the index of the next angle in angle_data (ind1), find the
    # corresponding angles from the angle_data vector (y0 and y1) and use the
    # interpolate function to find yy, which will be the angle of firing for a
    # given bin (the function is described in more detail below).
    
    binned_awake_angle_data=rep(0, ceiling((awake_end-awake_start)*1000));
    for(i in 1:length(binned_awake_angle_data)){
        ind0=ceiling(i/(1000/angle_freq)) ### SHOULD THIS BE FLOOR???
        ind1=ind0+1;
        if(ind1>length(awake_angle_data)){
            next
        }
        y0=awake_angle_data[ind0]
        y1=awake_angle_data[ind1]
        binned_awake_angle_data[i]=interpolateBetweenPoints(i,ind0,ind1,y0,y1,angle_freq)
    }
    
    return(list(cellnames=cellnames,
                spiketimes=spiketimes,
                awake_angle=binned_awake_angle_data))
}


