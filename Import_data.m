function[cellnames,spiketimes,binnedAwakeAngleData,resclu]=Import_data(dir_nm,numfiles)
%% Importing raw data
 
% Create empty cell array to store res and clu data from each tetrode
% tetrode is a multiple electron measurement device.
resclu = cell(numfiles,2);
% res data stores the timestamps corresponding to when a cell fired.
% clu data stores the cell numbers corresponding to each timestamp.
for i = 1:numfiles
    % Set file names:
    clufile = sprintf('%s.clu.%d',dir_nm,i); % The %s uses the dir_nm input and the %d uses the i input
    resfile = sprintf('%s.res.%d',dir_nm,i);
    clusters = load(clufile);
    timestamps = load(resfile);
% Assign values to the previously defined cell array resclu. Cluster
% arrays are assigned to the ith row and 1st column and timestamps
% are assigned to the ith row and second column. Curly brackets are
% used because resclu is a cell array.
    resclu{i,1} = clusters;resclu{i,2} = timestamps;
end


%% Cell Names and Timestamps

awake=load(sprintf('%s.states.Wake',dir_nm));
awakeStart = awake(1,1);
awakeEnd = awake(end,end);

electrodeFreq=20000;  %Electrophysiologal rate (i.e. the rate at which the raw data were recorded; Hz)

% Create empty string array to store names of cells, defined by tetrode
% and cell number
cellnames=string();
%%% Create empty cell array to store spike times.
spiketimes_All_Cells={};
for tetrode=1:length(resclu) %%% Length of resclu variable will be the same as numfiles
    timestamps_ms=1000*resclu{tetrode,2}/electrodeFreq; % Time of firings in all cells measured by electrode. 
    %Dividing timestamps by electrode freq rescales into seconds, multiplying by 1000 then rescales to milliseconds.
    TetrodeNumber=resclu{tetrode,1};
    TetrodeNumber = TetrodeNumber + 1; % The first cell in the tetrode is marked as 0, but matlab cannot handle an index of 0, so we add 1 to all cell numbers.
    NuofCellsInTetrode=TetrodeNumber(1)-1; %The first element in the cluster is an overview of nr of Cells measured. We minus 1 because we previously added one and to get the true number of cells we must minus that again.
    cellnames_per_tetrode=string();
    for cellnum=1:NuofCellsInTetrode
        cellnames_per_tetrode{cellnum,1} = sprintf('T%dC%d', tetrode, cellnum);
        thisCellIndex=find(TetrodeNumber==cellnum)-1; %%%% Must be -1 because first index of Tetrode Number does not correspond to a real cell number
        TimestampsThisCell=timestamps_ms(thisCellIndex);
        %goes and finds the cell timestamp at the input index, as decided by the previous code
        %The -1 makes sure it doesn't find index one because that is not the name of a cell.
        %Since it's a for loop, this will create an array for each cell
        %individually (and iterates for each cell in each cluster because
        %it's a for loop).
        IndexAwakeTimestampsThisCell=find((TimestampsThisCell>(awakeStart*1000)) & (TimestampsThisCell<(awakeEnd*1000))); %find returns the Index of all times when the mouse is awake and therefore has usable firing data.
        awakeCellTimestampsThisCell=TimestampsThisCell(IndexAwakeTimestampsThisCell); %now only firings when the mouse is awake.
        awakeCellTimestampsThisCell = floor(awakeCellTimestampsThisCell-awakeStart*1000.); %%%% Rescale awake data to starting time of 0 with whole integers.
        spiketimesThisTetrode{cellnum,1}=awakeCellTimestampsThisCell;
    end
    if tetrode == 1
        cellnames=cellnames_per_tetrode;
    else
        cellnames=[cellnames;cellnames_per_tetrode];
    end
    
    spiketimes_All_Cells=[spiketimes_All_Cells;spiketimesThisTetrode];
    spiketimesThisTetrode={}; %variable must be cleared on every iteration of tetrode for loop.
end
spiketimes=containers.Map({cellnames{1:end}},{spiketimes_All_Cells{1:end}});
%% Angle Data
angleFreq=1250/32;    % Angle measurement rate (1,250 is the EEG sampling rate and 32 the
                       %number of LFP samples; Hz/sample)* 

angleData=sprintf('%s-HD-filt.ang',dir_nm); 
angleData=load(angleData);

%%% -1 was used as a code for missing data. In order for Matlab to
%%% understand that these values are missing, we must recode them to NaN
%%% (Not a Number)
tmpIndex1=find(angleData==-1);
angleData(tmpIndex1)=NaN;

% Find the index for a given time to find the corresponding angle in
% angleData. round() rounds to the nearest integer (away from zero for
% > 0.5). Units: s * (Hz/sample) = 1/sample (i.e. index of a given sample).
% Then, use the indices to select the angleData recorded from the awake
% mouse.
indexAwakeStart=round(awakeStart*angleFreq);
indexAwakeEnd= round(awakeEnd*angleFreq);
awakeAngleData=angleData(indexAwakeStart:indexAwakeEnd);

%% BIN THE ANGLE DATA
% Creating binnedAngleData, a vector containing one angle measurement
% for each awake ms, because the length of the binnedAngleData is
% given by the number of milliseconds between awakeStart and awakeEnd.
% ceil() rounds a decimal number upwards, regardless of the value of the decimal.

% For each index i of binnedAngleData, find the corresponding index in
% the awakeAngleData vector by dividing the index by angleFreq in mHz/sample.
% Then, find the index of the next angle in angleData (ind1), find the
% corresponding angles from the angleData vector (y0 and y1) and use the
% interpolate function to find yy, which will be the angle of firing for a
% given bin (the function is described in more detail below).

binnedAwakeAngleData=zeros(ceil((awakeEnd-awakeStart)*1000),1);
for i=1:length(binnedAwakeAngleData)
    ind0=ceil(i/(1000/angleFreq)); %%% SHOULD THIS BE FLOOR???
    ind1=ind0+1;
    if ind1>length(awakeAngleData)
        continue
    end
    y0=awakeAngleData(ind0);
    y1=awakeAngleData(ind1);
    binnedAwakeAngleData(i)=interpolateBetweenPoints(i,ind0,ind1,y0,y1,angleFreq);
end

end