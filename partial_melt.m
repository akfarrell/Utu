function [ output_args ] = partial_melt( amps, delay)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
addpath('./matTaup/')
        tolerance = 0.0000001;
        for i=1:numel(dnum)
            t = find(dnum(i)>(time_value-tolerance) & dnum(i)<(time_value+tolerance));
            if t~=0
                index = i;
            end
        end
end

