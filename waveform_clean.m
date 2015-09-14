function w_clean = waveform_clean(w_raw, filter_object)
% waveform_clean  Clean waveform data
%   w_clean = waveform_clean(w_raw) cleans raw waveform data
%
% w_clean = waveform_clean(w_raw)
%
%   Inputs:
%       w_raw - object with raw waveform data
%       filter - filter object of the format filterobject(filter, values)
%   Output:
%       w_clean - cleaned waveform data that has been detrended, tapered,
%       gaps filled, filtered, and integrated
%
%   Author: Alexandra Farrell 2014/11/12 from work by Glenn Thompson
    w_cleaned = detrend(w_raw); %remve trend to center the data on zero
    %w_cleaned = taper(w_cleaned, 0.2); %taper data to join data so it can repeat indefinitely for FFT analysis
    w_cleaned = fillgaps(w_cleaned, 0); %fill gaps with a zero
    if nargin < 2
        w_clean = w_cleaned;
    else
        w_cleaned_filtered = filtfilt(filter_object, w_cleaned);
        w_clean = w_cleaned_filtered;
    end
    
    
    
%     f = filterobject('h', 0.1, 2); %usage f = filterobject(type, cutoff, poles), create a variable with filtering information
%     f = filterobject();
%     w_cleaned_filtered = filtfilt(f, w_cleaned); %apply the filter to the data
    %w_cleaned_filtered_integrate = integrate(w_cleaned_filtered); %integrate data to get displacement
    %w_clean = w_cleaned_filtered_integrate; %create returned variable
    %w_clean = w_cleaned;
    
    
    %w_clean = w_cleaned_filtered;
    
    
    
end