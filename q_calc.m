function [stations,Q] = q_calc(waveform_object, velocity, depth, frequency)
% q_calc  Calculate q_values using the amplitudes at stations. The largest
% amplitude is considered A0, so that station has an infinite q-value
%   q_vals = q_calc(waveform_object, velocity, frequency) calculates q-values from amplitudes
%
% q_vals = q_calc(waveform_object, velocity, frequency)
%
%   Inputs:
%       waveform_object - waveform object with amplitude values and stations from
%       velocity - average velocity of crust for the subsurface of
%       interest, in km/s
%       depth - length of subsurface area of interest, in km
%       frequency - frequency, in Hz, of data sample
%   Outputs:
%       Q - len(station) matrix of crustal Q values for the velocity and depth specified
%       stations - 4*station strings of individual letters in all of the
%       stastions
%
%   Author: Alexandra Farrell 2015/5/18 from formula in Jochen Braunmiller
%   comps question

    ln = @log;
    stations = char([get(waveform_object, 'station')]);
    amp_abs = [get(waveform_object, 'AMP_ABS')]
    max_val = abs(max(max(amp_abs), min(amp_abs)));
    Q = [];
    for i = 1:numel(amp_abs)
        %try
        %catch err
            Q(i) = (-pi*(depth/velocity)*frequency(i))/(ln(amp_abs(i)/max_val));
        %end
    end
