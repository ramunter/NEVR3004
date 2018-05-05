# From raw data to neurons

*Topics*

1. Recording the brain.
2. Processing the data.
3. Spike sorting.
4. A quick overview of imaging.

## Recording the brain

### What do we use to record the brain?

Tetrode $\rightarrow$ cheap but only 4 channels
Silicon probe $\rightarrow$ very expensive but up to several hundred channels.  

Tetrode is used for long time data collection, probes are used for more short term experiments.

Spike sorting consits of finding what neuron is being recorded by comparing the potentials on different channels.

## Processing the data

### Data analysis pipeline

1. Filter data.
2. Spike detection.
3. Feature extraction.
4. Clustering

Scale: Microvolt, picoampere.

### Filtering

Data recordings have high-frequency noise, and by low-frequency  variations due to network activity.

* A low-pass filter* is a filter that passes signals with a frequency lower than a certain frequency.

* A high-pass filter* is a filter that passes signals with a frequency higher than a certain frequency.

* Band-pass filter* keeps a certain range of frequencies

It is always better to record raw data, leads to less distortion and can use the raw data as a reference.

## Spike detection

### Detecting the spikes

1. Identify the peaks that cross a threshold(usually 3-5 sd from baseline) on the filtered data.
2. Extract the values before and after the peak (1ms) to extract the wave.

Will still have some noise. Usually then use the template method. The selected peaks are compared to some wavesforms' template fitted on the data.

### Spike sorting

A *potential unit* is what a neuron is called before we are sure it is a neuron.

We now have a bunch of spikes, but these come from different cells. What properties should the spikes from a _potential unit_ share? 

* Waveforms
    - Similar ampllitude
    - Similar valley
    - Similar width
* Temporal autocorrelation
    - Must respect the refractory period.
* Properties of the spikes

Neurons tend to have very stable waveforms though they can vary a lot from neuron to neuron. Plotting the spiking(firing rate) can help aswell.

## Clustering

Different manipulations of the data can help identify the clusters.

* PCA to reduce dimensionality and highlight the most variable components.
* Valley, length of the waveforms.
* *Temporal Correlation*

Different algorithms:

* Mountainsort
* Kilosort
* Spiking circus

Can be done by hand with little data.

### Cluster quality

After grouping the similar spikes into clusters, wee need to check their temporal correlation. We compute the autocorrelation. This is apparently counting the time between activation firings for every activation firing. Plotting a histogram of these times, one should have 0 for the 2ms firing window(not physically possible). If 1% of the firings have this contamination, don't include it in analysis.  

For quantitive measures we can use is the average energy of the the spikes.  

We can then use the L-ratio test or the isolation distance to see how different your clusters are from the noise. 

*Add these measures from the slides*

## Imaging

Neuronal actiivty results to an increase in intracellular calcium. You can inject things which change color with this calcium variation, or use GM for the same effect. Then you can image the calcium levels to find neuron spikes.

# Stimulus-response

Information is processed with *action potentials*. At first we will not consider the content of the AP, and simply use them as _timestamps_ of ecnoded information.

## Quantifying the stimilus response process

*Step 1:* calculate the firing rate, i.e. the number of spikes emitted in a pre-defined time window, normalized by the length of the time window.

$$r = \frac{n}{T}$$

where $r$ is the FR, $n$, is the number of spikes and, $T$ is the time window.  

*In a S-R relationship* the response of a neuron varies with the stimulus. We want to find this link between stimulus and firing rate.  

*Step 2:* Plot the firing rate as a function of the stimulus attribute. THe analyses strongly depens on teh observed stimulus.  

It can be good to fit a curve:

* Characterizes the cell type
* Allows to talk about precision, accuracy
* Allows a better qualification and quantification of the stimulus-response behavioural
* Permits potential clustering of cells  with similar properties

# Information theory 1

*The neural code* is the minimum set of symbols [most compact description of the neural activity] capable of representing all of the bilogically significant information.

We use information theory to identify the neural code.

* single neuron vs population code
* Spike count coding hypothesis
* Spike time coding hypothesis 
    - Spike pattern
    - Synchrony code

Also

* Speed of information transmission
* Finding the stimulus-response relationship
* Infering effective connectivity

## Stochastic framework

$\text{Stimulus}\:S\:\rightarrow\:\text{population of neurons}\:\rightarrow\:\text{Response}\:R$

This pipeline is stochastic because the stimulus is stochastic and there is noise in the neuron stage, leading to a distribution of different responses.  

The encoding problem is giving a response given a stimulus($p(r|s)=\frac{p(s,r)}{p(s)}$). The decoding problem is doing the opposite, and uses the theorem to calculate stimulus based on the response.

$$p(s|r) = \frac{p(r|s)p(s)}{p(r)}$$

## Entropy

*Entropy* is the uncertainty of the random variable $s$ and thus the information gained by measuring it.

$$H(S) = -\sum_sp(s)log_2p(s)$$

I.e. a random variable that is not random(a coin that always returns heads) has no entropy as we know before flipping it what it is going to be. A fair coin has entropy equal to 1 as we cannot predict its result so we gain "1 bit information" for flipping it.

Correlation is a linear measure of relations between two variables. But if two variables are related through a nonlinear relationship, the correlation can be zero. Mutual information captures nonlinear relationships.  

*Mutual information* is a measure of how large a reduction in uncertainty there is for one variable given the value of the other variable.

$$I(S;R)=\sum_{s,r}p(s,r)log_2\frac{p(s,r)}{p(s)p(r)}$$

$$I(S;R) = H(S) - H(S|R)$$

*Note*

* MI is symetric, we cannot say anything about causality based on it.
* $I(S;R) \ge I(S;f(R))$

# Spike Train Analysis

* Firing rate isn't everything, the temporal relation is important
* Therefore we prefer looking at the firing rate as a function of time
    - Bin firing rate over small $dt$
    - Sliding bin firing rate over small $dt$
    - Sliding bin with gaussian waiting over small $dt$
* When recording neuronal actiivty, closer neurons cause a low frequency potential oscilation. These are called waves(i.e. theta wave)

* Independent neuron firings should follow a Poisson process.
* *Inter-spike interval* - the time between action potentials.
* You can plot the inter-spike interval histogram and compare it to a Poisson process to check if the neuron firings are independent.
* *Burst neurons* will have a minimum between two peaks in the inter-spike histogram
* We can use the *coefficient of variation* to check the temporal variability
    - 1 for a poisson process
    - 0 for a constant firing at dt = some value
* *Fan factor* is the normalized standard deviation of the toal number of spikes from tial to trial. 
    - 1 for a poisson process
* *Autocorrelation* is the correlation between the firing rate and a shifted version of the spike train. 
* *Spike triggered average* is the average stimulus that causes a spike
* *Burst triggered average* is the average stimulus that causes a spike when considering spikes that occur in a *burst* to be a different response than a single spike.

# I dont know what this is about




