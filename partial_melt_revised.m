%function per_melt = partial_melt_revised(delay, aoi, model_number,thickness)
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
 %clc       
%%   
model_number = 4; %--------You can change dis!!!!!-----------

% ----- Determining average P velocity ------- %
K = 45; %bulk modulus in GPa
u = 26; %rigidity in GPa
rho = 2.7; %density in Kg/m3
%Vp_solid = sqrt((K+1.33+u)/rho);
%Vp_melt = sqrt(K/rho);
%percent_melt = 0:100;
%Vp_average = (Vp_solid.*(100-percent_melt)+Vp_melt.*percent_melt)/100;
%additions = linspace(1,300,300);
% for i = 2:numel(additions)
%       varz(1) = Vp_average(end); 
%       varz(i) = varz(i-1)-(Vp_average(1)-Vp_average(2));
% end
% %%
% varz(end) = varz(end-1)+(Vp_average(2)-Vp_average(1));
% Vp_average = horzcat(Vp_average, varz);

%Vp_solid = 5.0825;
%Vp_melt = 4.0089;
x = linspace(0,1,101);
%Vp = 1./((x./Vp_melt)+(1.-x)/Vp_solid);
%Vp_average = Vp;

    perc_water = 0.07;
    meltz = 1-perc_water;
    K_water = 2.2; %bulk modulus of water in GPa
    mu_min = 0.1;
    mu_melt = 0.5;
    K_water = 2.2;
    x2=1-x;
    
    Ktot = 1./(x2./(K+4/3*mu_min) + (1-x2).*perc_water./(2.2+4/3*mu_min) + (1-x2).*meltz/(K+4/3*mu_min))-(4/3)*mu_min;
    Ktot_nowater = K;
    mu_tot = mu_melt + x2./((1/(u-mu_melt)) + (2.*x*(K+2*mu_melt)/(5*mu_melt*(K+(4/3)*mu_melt))));
    clear Vp_average
    rho_melt = 2.3;
    rho_used = rho_melt;
    mu_upper_bound = u+(x./((1/(mu_melt-u)+2*x2.*(K+2*u)/(5*u*(K+(4/3*u))))));
    %Vp_average = sqrt((Ktot_nowater+1.33+mu_upper_bound)./rho_used)
    %%uncomment for no water upper bound results
    %Vp_average = sqrt((Ktot_nowater+1.33+mu_tot)./rho_used) %uncomment for
    %no water results
    Vp_average = sqrt((Ktot+1.33+mu_tot)./rho_used) %uncomment for water
    %results
    
    thickness = 20;
    Vp_solid = Vp_average(1);
%%        
delay = delay_corrected;
%aoi = eq(earthquake_number).aoi;

%percent_melt = linspace(1,100,100);
% Vp for 2600 kg/m^3 density
% Vp_average = [5.26, 5.25, 5.24, 5.23, 5.22, 5.21, 5.20, 5.19, 5.17, 5.16...  Average P velocity in km/s
%               5.15, 5.14, 5.13, 5.12, 5.11, 5.10, 5.08, 5.07, 5.06, 5.05...
%               5.04, 5.03, 5.02, 5.01, 5.00, 4.98, 4.97, 4.96, 4.95, 4.94...
%               4.93, 4.92, 4.91, 4.90, 4.88, 4.87, 4.86, 4.85, 4.84, 4.83...
%               4.82, 4.81, 4.80, 4.78, 4.77, 4.76, 4.75, 4.74, 4.73, 4.72...
%               4.71, 4.70, 4.68, 4.67, 4.66, 4.65, 4.64, 4.63, 4.62, 4.61...
%               4.59, 4.58, 4.57, 4.56, 4.55, 4.54, 4.53, 4.52, 4.51, 4.49...
%               4.48, 4.47, 4.46, 4.45, 4.44, 4.43, 4.42, 4.41, 4.39, 4.38...
%               4.37, 4.36, 4.35, 4.34, 4.33, 4.32, 4.31, 4.29, 4.28, 4.27...
%               4.26, 4.25, 4.24, 4.23, 4.22, 4.20, 4.19, 4.18, 4.17, 4.16];
            
Distances = 0:100;
%Vp_solid = 5.27;
if model_number == 1 %Using just the velocities of the UTU velocity model, in same 20-km space as Ward anomaly
    %thickness = 25.3; %thickness of partial melt body, in km
    %thickness = 10;
    %thickness = 1;
    %thickness = 15;
    tolerance = 0.007;
    dist_in_magma = thickness/cosd(aoi); %distance ray travels through magma, in km
    delay_through_anomaly = delay + dist_in_magma/Vp_solid;%fastest_possible_time
    per_melt = [zeros(size(delay))];
    for s = 1:numel(delay)
        per_melt(s) = 111;
        for k = 1:numel(Vp_average)
            Vp_in = Vp_average(k);
            travel_time(k) = dist_in_magma/Vp_in;
            travel_time;
            delay_through_anomaly(s)-tolerance;
            del_val = find(travel_time(k)>(delay_through_anomaly(s)-tolerance) & travel_time(k)<(delay_through_anomaly(s)+tolerance));
            if del_val~=0
                per_melt(s) = k-1;
            end
            if k == numel(Vp_average)
                if per_melt(s) == 111;
                    if delay_through_anomaly(s) > travel_time(1)
                        per_melt(s) = 999;
                    elseif delay_through_anomaly(s) < travel_time(numel(Vp_average))
                        per_melt(s) = -1;
                    end
                end
            end
        end
    end
    per_melt
elseif model_number == 2 %estimating thickness with given Vp
    %Vp_in = Vp_average(21);
    %Vp_in = Vp_average(18);
    Vp_in = Vp_average(16);
    thickness = [zeros(size(delay))];
    %tolerance = 0.009;
    for s = 1:numel(delay)
        thickness(s) = 999;
        for k = 1:numel(Distances)
            tolerance = 0.0075; %for Vp_average(20)
            Dist_try = Distances(k)*cosd(aoi);
            travel_time(k) = Dist_try/Vp_in;
            travel_time;
            delay_through_anomaly = delay + Dist_try/Vp_solid; %fastest_possible_time
            del_val = find(delay_through_anomaly(s)>(travel_time(k)-tolerance) & delay_through_anomaly(s)<(travel_time(k)+tolerance));
            if del_val~=0
                thickness(s) = k-1;
            end
            if k == numel(Distances)
              if thickness(s) == 999
                  tolerance = 0.03;
                  for k = 1:numel(Vp_average)
                      %Dist_try
                      Dist_try = Distances(k)*cosd(aoi);
                      travel_time(k) = Dist_try/Vp_in;
                      travel_time;
                      delay_through_anomaly = delay + Dist_try/Vp_solid; %fastest_possible_time
                      del_val = find(delay_through_anomaly(s)>(travel_time(k)-tolerance) & delay_through_anomaly(s)<(travel_time(k)+tolerance));
                      if del_val~=0
                          thickness(s) = k-1;
                      end
                  end
              end
            end
        end
    end
    mean(thickness)
    thickness
elseif model_number == 3 %Using Heather's depths
    thickness = abs(-25-depths_at_stas); %need to run get_data_fig_surf.m first to get depths at stas!!!!
    %tolerance = 0.006;
    for i = 2:numel(Vp_average)
        tolerance(i-1) = (Vp_average(i-1)-Vp_average(i))/2;
        if i == numel(Vp_average)
            tolerance(i) = tolerance(i-1);
        end
    end
    
    Vp_solid = Vp_average(1);
    dist_in_magma = thickness./cosd(aoi); %distance ray travels through magma, in km
    %delay_through_anomaly = delay + dist_in_magma./Vp_solid;%fastest_possible_time
    per_melt = [zeros(size(delay))];
    for s = 1:numel(delay)
        per_melt(s) = 111;
        for k = 1:numel(Vp_average)
%             if per_melt(s) ~= 111
%                 continue
%             end
            Vp_in = Vp_average(k);
            travel_time(k) = dist_in_magma(s)/Vp_in;
            travel_time;
            delay_through_anomaly = delay(s) + dist_in_magma(s)/Vp_solid;%fastest_possible_time
            delay_through_anomaly-tolerance;
            del_val = find(travel_time(k)>(delay_through_anomaly-tolerance(k)) & travel_time(k)<(delay_through_anomaly+tolerance(k)));
            if del_val~=0
                per_melt(s) = k-1;
            end
            if k == numel(Vp_average)
                if per_melt(s) == 111;
                    if delay_through_anomaly > travel_time(1)
                        per_melt(s) = 999;
                    elseif delay_through_anomaly < travel_time(numel(Vp_average))
                        per_melt(s) = -1;
                    end
                end
            end
        end
    end
    %per_melt
elseif model_number == 4 %Using Hashin-Shtrikman Bounds
    perc_water = 0.07;
    meltz = 1-perc_water;
    K_water = 2.2; %bulk modulus of water in GPa
    mu_min = 0.1;
    mu_melt = 0.5;
    K_water = 2.2;
    x2=1-x;
    
    Ktot = 1./(x2./(K+4/3*mu_min) + (1-x2).*perc_water./(2.2+4/3*mu_min) + (1-x2).*meltz/(K+4/3*mu_min))-(4/3)*mu_min;
    Ktot_nowater = K;
    mu_tot = mu_melt + x2./((1/(u-mu_melt)) + (2.*x*(K+2*mu_melt)/(5*mu_melt*(K+(4/3)*mu_melt))));
    clear Vp_average
    rho_melt = 2.3;
    rho_used = rho_melt;
    mu_upper_bound = u+(x./((1/(mu_melt-u)+2*x2.*(K+2*u)/(5*u*(K+(4/3*u))))));
    %Vp_average = sqrt((Ktot_nowater+1.33+mu_upper_bound)./rho_used)
    %%uncomment for no water upper bound results
    %Vp_average = sqrt((Ktot_nowater+1.33+mu_tot)./rho_used) %uncomment for
    %no water results
    Vp_average = sqrt((Ktot+1.33+mu_tot)./rho_used); %uncomment for water
    %results
    
    thickness = 20;
    dist_in_magma = thickness/cosd(aoi); %distance ray travels through magma, in km
    Vp_solid = Vp_average(1);
    delay_through_anomaly = delay + dist_in_magma/Vp_solid;%fastest_possible_time
    per_melt = [zeros(size(delay))];
    for s = 1:numel(delay)
        per_melt(s) = 111;
        for k = 1:numel(Vp_average)
            for i = 2:numel(Vp_average)
                tolerance(i-1) = (Vp_average(i-1)-Vp_average(i))/2+1/(i*10);
                if i == numel(Vp_average)
                    tolerance(i) = tolerance(i-1);
                end
            end
            Vp_in = Vp_average(k);
            travel_time(k) = dist_in_magma/Vp_in;
            travel_time;
            delay_through_anomaly(s)-tolerance(s);
            del_val = find(travel_time(k)>(delay_through_anomaly(s)-tolerance(s)) & travel_time(k)<(delay_through_anomaly(s)+tolerance(s)));
            if del_val~=0
                per_melt(s) = k-1;
            end
            if k == numel(Vp_average)
                if per_melt(s) == 111;
                    if delay_through_anomaly(s) > travel_time(1)
                        %per_melt(s) = 999;
                        for k = 1:numel(Vp_average)
                            Vp_in = Vp_average(k);
                            travel_time(k) = dist_in_magma/Vp_in;
                          tolerance = (Vp_average(1)-Vp_average(2))/2;
                          delay_through_anomaly(s)-tolerance;
                          del_val = find(travel_time(k)>(delay_through_anomaly(s)-tolerance) & travel_time(k)<(delay_through_anomaly(s)+tolerance));
                          if del_val~=0
                              per_melt(s) = k-1;
                          end
                        end
                    elseif delay_through_anomaly(s) < travel_time(numel(Vp_average))
                        per_melt(s) = -1;
                    end
                end
            end
        end
    end
end
    stationzz = {'PL07', 'PLAR', 'PLBR', 'PLLL', 'PLSS', 'PLTP'} %KTSZ2
    %stationzz = {'PLAR', 'PLBR', 'PLLL'} %JSZ3
    %stationzz = {'PLAR', 'PLLC', 'PLLL', 'PLSS'} %SSSZ4
    %per_melt
    %min(per_melt)
    %max(per_melt)
    
    for i = 1:numel(stationzz)
         per_melt(strcmp(sta,stationzz{i}))
         %thickness(strcmp(sta,stationzz{i}))
    end

