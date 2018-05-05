function [firingRate,plottingangles,anglesAtFiring] = Calculate_Firing_Rate(spiketimes,angleData)

%% Define Occupancy

% Occupancy is needed to compute the tuning curve and is defined for each
% bin as the number of elements in binnedAwakeAngleData between two given
% angles AA and BB. The values of AA and BB is given by dividing the
% complete circle (2*pi) into the number of bins (nrpoints) and then
% multiply by the previous bin number or the current bin number (to
% find an angle interval).

% As such, in occupancy you get a "measure" of how much time was spent in
% each interval; it is the total count of milliseconds the animal spends
% in a particular bin. This follows as binnedAwakeAngleData expand the
% awakeAngleData vector to contain a measure for each millisecond. Thus, if
% an animal spends several milliseconds in a particular angle (e.g. when
% staring at a cue), the subsequent elements in the vector will be equal.
nrpoints=40;
occupancy=zeros(1,nrpoints);
plottingangles=zeros(1,nrpoints);
for i=1:nrpoints
    AA = pi*2/nrpoints*(i-1);
    BB = pi*2/nrpoints*i;
    plottingangles(i) = 0.5*(AA+BB)*360./(2.*pi);
    occupancy(i)=length(find(AA<angleData & angleData<BB));
end

%% Find Which Angles Correspond to Firing Times

% Make a row vector with zeros to contain the angles corresponding to the time of firing.
anglesAtFiring=zeros(length(spiketimes),1);
%%%%find the angle that corresponds to each spike time.
for i=1:length(spiketimes)
    if spiketimes(i) == 0
        indexAngle=1;
    else
        indexAngle=spiketimes(i);
    end
anglesAtFiring(i)=angleData(indexAngle);
end
%% Use anglesAtFiring to Make Firing Rates For Given Angle Bins
% Now that the angles at which a given cell is firing have been computed,
% anglesAtFiring, they can be used to compute the firing rates necessary to plot the tuning cuves.

% Make a row vector firingRate of zeros with the total number of bins (nrpoints)
firingRate=zeros(1,nrpoints);

% BINNING to find how many spikes are within each radian interval.

% For each bin:
% AA will contain the lower limit of the given bin interval, BB will
% contain the upper limit of the same bin. The limits are given in
% radians (because we are to bin the angles at which the cell fired).
% The number of spikes, numspikes, contained in the given interval is
% found. Then the firing rate for the particular bin is
% computed, dividing the number of spikes in the interval by the
% occupancy. Recall that each element in occupancy is the total count
% of milliseconds the animal spends in a particular bin. Multiply by
% 1,000 because occupancy is in milliseconds. Hence, the unit of the
% tuning curve's firing rate becomes (spikes * 1000) / (1000 * s) = spikes/s.

for i=1:nrpoints
    AA = pi*2/nrpoints*(i-1);
    BB = pi*2/nrpoints*i;
    numspikes=length(find(AA<anglesAtFiring & anglesAtFiring<BB));
    firingRate(i)=1000. * numspikes/occupancy(i); % Convert to firing rates
end


end