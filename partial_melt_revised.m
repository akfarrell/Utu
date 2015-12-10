%function per_melt = partial_melt_revised(delay, aoi, model_number)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Read from file        
% fid = fopen('KTSZ1/text/KTSZ1_output_diffFreq_0.1875_0.3750.txt');
%         C = textscan(fid, '%s %10.5f %10.3f %10.4f');
%         fclose(fid);
%         f = {'Station','Q','Freq', 'Time_offset'};
%         values = cell2struct(C,f,2);
%         delay = values.Time_offset;
%         amps = max(values.Q)-values.Q;
        

        %tolerance for closest velocity to match slowness
delay = delay_corrected;
aoi = eq(earthquake_number).aoi;
model_number = 1;
percent_melt = linspace(1,100,100);

Vp_average = [5.26, 5.25, 5.24, 5.23, 5.22, 5.21, 5.20, 5.19, 5.17, 5.16...  Average P velocity in km/s
              5.15, 5.14, 5.13, 5.12, 5.11, 5.10, 5.08, 5.07, 5.06, 5.05...
              5.04, 5.03, 5.02, 5.01, 5.00, 4.98, 4.97, 4.96, 4.95, 4.94...
              4.93, 4.92, 4.91, 4.90, 4.88, 4.87, 4.86, 4.85, 4.84, 4.83...
              4.82, 4.81, 4.80, 4.78, 4.77, 4.76, 4.75, 4.74, 4.73, 4.72...
              4.71, 4.70, 4.68, 4.67, 4.66, 4.65, 4.64, 4.63, 4.62, 4.61...
              4.59, 4.58, 4.57, 4.56, 4.55, 4.54, 4.53, 4.52, 4.51, 4.49...
              4.48, 4.47, 4.46, 4.45, 4.44, 4.43, 4.42, 4.41, 4.39, 4.38...
              4.37, 4.36, 4.35, 4.34, 4.33, 4.32, 4.31, 4.29, 4.28, 4.27...
              4.26, 4.25, 4.24, 4.23, 4.22, 4.20, 4.19, 4.18, 4.17, 4.16];

TT_for_20_km = [3.80, 3.81, 3.82, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.87...   Travel time for 20 km at x percent partial melt
                3.88, 3.89, 3.90, 3.91, 3.92, 3.92, 3.93, 3.94, 3.95, 3.96...
                3.97, 3.98, 3.99, 3.99, 4.00, 4.01, 4.02, 4.03, 4.04, 4.05...
                4.06, 4.07, 4.08, 4.09, 4.09, 4.10, 4.11, 4.12, 4.13, 4.14...
                4.15, 4.16, 4.17, 4.18, 4.19, 4.20, 4.21, 4.22, 4.23, 4.24...
                4.25, 4.26, 4.27, 4.28, 4.29, 4.30, 4.31, 4.32, 4.33, 4.34...
                4.35, 4.36, 4.37, 4.38, 4.40, 4.41, 4.42, 4.43, 4.44, 4.45...
                4.46, 4.47, 4.48, 4.49, 4.51, 4.52, 4.53, 4.54, 4.55, 4.56...
                4.57, 4.59, 4.60, 4.61, 4.62, 4.63, 4.65, 4.66, 4.67, 4.68...
                4.69, 4.71, 4.72, 4.73, 4.74, 4.76, 4.77, 4.78, 4.79, 4.81];
tolerance = 0.005;   
Vp_solid = 5.27;
if model_number == 1 %Using just the velocities of the UTU velocity model, in same 20-km space as Ward anomaly
    thickness = 70; %thickness of partial melt body, in km
    dist_in_magma = thickness/cosd(aoi); %distance ray travels through magma, in km
    top_distance = (15-4)/cosd(aoi); top_velocity = 4.1; %distance in top part of velocity model in km and velocity of this space in km/s
    bottom_distance = (25-15)/cosd(aoi); bottom_velocity = 6.4; %distance in bottom part of velocity model in km and velocity of this space in km/s
    fastest_possible_time = top_distance/top_velocity + bottom_distance/bottom_velocity
    delay_through_anomaly = delay + dist_in_magma/Vp_solid%fastest_possible_time
    per_melt = [zeros(size(delay))];
    for s = 1:numel(delay)
        per_melt(s) = 999;
        for k = 1:numel(Vp_average)
            Vp_in = Vp_average(k);
            travel_time(k) = dist_in_magma/Vp_in;
            travel_time;
            delay_through_anomaly(s)-tolerance;
            del_val = find(travel_time(k)>(delay_through_anomaly(s)-tolerance) & travel_time(k)<(delay_through_anomaly(s)+tolerance));
            if del_val~=0
                per_melt(s) = k;
            end
        end
    end

end
per_melt