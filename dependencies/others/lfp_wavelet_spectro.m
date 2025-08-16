function [scalogram,frequency,t]=lfp_wavelet_spectro(lfp,frequencyLimits)
%Compute scalogram

%Parameters
sampleRate = 1000;
wavelet = 'amor';
voicesPerOctave = 10;
%frequencyLimits = [0,30];

%Compute time vector
t = 0:1/sampleRate:(length(lfp)*1/sampleRate)-1/sampleRate;

%Compute CWT
%If necessary, substitute workspace variable name for lfp as first input to cwt() function in code below
%Run the function call below without output arguments to plot the results
[waveletTransform,frequency] = cwt(lfp', sampleRate, wavelet,...
    VoicesPerOctave = voicesPerOctave,...
    FrequencyLimits = frequencyLimits);
scalogram = abs(waveletTransform);