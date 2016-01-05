%%
% for j = 1:numel(sta)
%     %[sta_dist(j), az(j)] = distance(lat_ref, lon_ref, lat_sta(j), lon_sta(j));
%     wave_to_sta_dist(j) = (abs(A*(lat_sta(j))+B*(lon_sta(j))+C)/sqrt(A^2 + B^2));
%     if earthquake_number >=6 && earthquake_number<=9
%         y_vals = sind(360-az_volc_eq-90)*wave_to_sta_dist(j);
%         x_vals = cosd(360-az_volc_eq-90)*wave_to_sta_dist(j);
%         LAT_point_on_line(j) = lat_sta(j)-y_vals;
%         LON_point_on_line(j) = lon_sta(j)-x_vals;
%     elseif earthquake_number >=2 && earthquake_number<=5
%         x_vals = sind(90-(270-az_volc_eq))*wave_to_sta_dist(j);
%         y_vals = cosd(90-(270-az_volc_eq))*wave_to_sta_dist(j);
%         LAT_point_on_line(j) = lat_sta(j)+y_vals;
%         LON_point_on_line(j) = lon_sta(j)-x_vals;
%     elseif earthquake_number >=10 && earthquake_number<=14 %%%%%%%%%%CHANGE!
%         x_vals = sind(90-(270-az_volc_eq))*wave_to_sta_dist(j);
%         y_vals = cosd(90-(270-az_volc_eq))*wave_to_sta_dist(j);
%         LAT_point_on_line(j) = lat_sta(j)+y_vals;
%         LON_point_on_line(j) = lon_sta(j)-x_vals;
%     end
% 
% end
%distOut = distdim(sta_dist, 'degrees', 'km');
%wave_to_sta_dist_OUT = distdim(wave_to_sta_dist, 'degrees', 'km');
%partial_melt(m_values, time_vals_ref, distOut, eq(earthquake_number).aoi, elev)

ex = linspace(-22.75, -21.75, 100);
for i = 1:numel(ex)
    ye(i) = A*ex(i)+C;
end


% tolerance = 0.0005;
% for index = 1:numel(utu_lat)
%     for index2 = 1:numel(topo_data.lat)
%         Lat_for_z = find(topo_data.lat(index2)>(utu_lat(index)-tolerance) & topo_data.lat(index2)<(utu_lat(index)+tolerance));
%         Lon_for_z = find(topo_data.lon(index2)>(utu_lon(index)-tolerance) & topo_data.lon(index2)<(utu_lon(index)+tolerance));
%         if Lat_for_z~=0
%             Lat_index_for_z(index) = index2
%         end
%         if Lon_for_z~=0
%             Lon_index_for_z(index) = index2
%         end
%     end
% end

%floor((68-67.186035)/(lat(1,1)-lat(2,1)))-1     to determine matrix coordinates
%for given lat and lon

% for index=1:numel(sta)
%     elev_on_line(index) = interp2(lon, lat, topo_data.z, LON_point_on_line(index), LAT_point_on_line(index)); %Interpolate to find elevation at points on the line
% end


%Calculate distance from plane to station with topography
%elev_on_line, LAT_point_on_line, LON_point_on_line, lat_sta, lon_sta,
%elev

%----RUN AS A TEST-------%
%eq(earthquake_number).aoi = 0.0000001;
%c = (elev - elev(1))
%Compare c to Distance

close all
[x_on_line, y_on_line]= ll2utm(LAT_point_on_line, LON_point_on_line); %calculate utm coordinates of points on line, in m


diff_on_line = 1000/tand(eq(earthquake_number).aoi);
if earthquake_number >=6 && earthquake_number<=9
    third_point_x = x_of_sta(1)+diff_on_line*(sind(az_volc_eq-270));
    third_point_y = y_of_sta(1)-diff_on_line*(sind(az_volc_eq-270));
elseif earthquake_number >=2 && earthquake_number<=5
    third_point_x = x_of_sta(1)+diff_on_line*(sind(90-(270-az_volc_eq)));
    third_point_y = y_of_sta(1)-diff_on_line*(sind(90-(270-az_volc_eq)));
elseif earthquake_number >=10 && earthquake_number<=14 %%%CHANGE!
    third_point_x = x_of_sta(1)+diff_on_line*(sind(90-(270-az_volc_eq)));
    third_point_y = y_of_sta(1)-diff_on_line*(sind(90-(270-az_volc_eq)));
end

[third_point_lat, third_point_lon] = utm2ll(third_point_x, third_point_y, -19); %convert this third point from UTM to lat lon

P1 = [x_of_sta(1),y_of_sta(1), elev(1)]; %first station coordinates (UTM) and elevation
P2 = [x_on_line(2), y_on_line(2), elev(1)]; %point (UTM) on intersecting line at elevation of first station
P3 = [third_point_x, third_point_y, (elev(1)-10)]; %point along the plane under the surface of the earth
%determine distance between points on line and station
%line_to_sta_dist= sqrt((x_on_line-x_of_sta).^2+(y_on_line-y_of_sta).^2); %distance from line to station on surface

normal = cross(P1-P2, P1-P3);

x = [P1(1) P2(1) P3(1)];
y = [P1(2) P2(2) P3(2)];
z = [P1(3) P2(3) P3(3)];

A = normal(1); B = normal(2); C = normal(3);
D = -dot(normal, P1);
nnorm = normal/norm(normal);
p = D/norm(normal);
for i = 1:numel(sta)
    P4 = [x_of_sta(1), y_of_sta(i), elev(i)];
    Distance(i) = dot(nnorm, P4)+p;
end

if strcmp(eq(earthquake_number).name, 'JSZ4')
    Distance_corr = Distance - Distance(3);
elseif strcmp(eq(earthquake_number).name, 'KTSZ2')
    Distance_corr = Distance - Distance(3);
elseif strcmp(eq(earthquake_number).name, 'SSSZ1') %%%%%%%%%%CHANGE!!!
    Distance_corr = Distance - Distance(3);
end


h = figure;
latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.5]; %[western_limit eastern_limit]
set(h, 'Position', [1000 1000 1000 1000])
worldmap(latlim, lonlim);
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
hold on
invQ = [];
maxval = max(Q);
min2 = sort(Q(:));
minval = min2(2);
total = maxval+minval;
plotm(x, y, 'k')
hold on
scatterm(lat_sta, lon_sta, '^', 'k')
%textm(LAT_point_on_line, LON_point_on_line, sta)
scatterm(LAT_point_on_line, LON_point_on_line, 'o', 'k', 'filled')
scatterm(lat(884,977), lon(884,977), '^', 'filled')
scatterm(lat_ref, lon_ref, '^', 'k', 'filled')
scatterm(third_point_lat, third_point_lon, 'm', 'filled')
scatterm(lat_sta(1), lon_sta(1), 'c', 'filled')
scatterm(LAT_point_on_line(2), LON_point_on_line(2), 'c', 'filled')
textm(lat_sta, lon_sta, sta)
hypotenuse = 0.07;
if earthquake_number >=6 && earthquake_number<=9
    lat_az = -21.85;
    lon_az = -67.60;
    u = hypotenuse*sind(360-eq(earthquake_number).az-90); %vertical
    v = hypotenuse*cosd(360-eq(earthquake_number).az-90); %horizontal
elseif earthquake_number >=2 && earthquake_number<=5
    lat_az = -22.7;
    lon_az = -67.60;
    u = hypotenuse*cosd(eq(earthquake_number).az); %vertical
    v = hypotenuse*sind(eq(earthquake_number).az); %horizontal
elseif earthquake_number >=10 && earthquake_number<=14
    lat_az = -22.7;
    lon_az = -66.95;
    u = hypotenuse*sind(eq(earthquake_number).az-90); %vertical
    v = hypotenuse*-cosd(eq(earthquake_number).az-90); %horizontal
end
quiverm(lat_az, lon_az,u, v, 'k')
directory = sprintf('/home/a/akfarrell/Uturuncu/%s/figures', eq(earthquake_number).name);
filename = sprintf('%s_wavefront_%1.4f_%1.4f.png',eq(earthquake_number).name,fil(1),fil(2));
filename_wPath = fullfile(directory,filename);
hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
%%
velocity = 4100; %velocity between -6 km and 15 km, in m/s
travel_time_relativeToFirstStation = Distance/velocity;

delay = time_vals_ref - travel_time_relativeToFirstStation;
corr_factor = min(delay);
delay_corrected = delay-corr_factor;
% 
% %partial_melt_percent = partial_melt_revised(delay_corrected, eq(earthquake_number).aoi, 1);


%velocity = [4100, 2500, 7000];
%
% for v_val = 1:numel(velocity)
%     travel_time_relativeToFirstStation_tries = Distance_corr/velocity(v_val);
%     delay_corrected(v_val, :) = time_vals_ref - travel_time_relativeToFirstStation_tries;
%     %corr_factor = min(delay);
%     %delay_corrected(v_val,:) = delay-corr_factor;
% end


%%

g = figure;
set(g, 'Position', [1000 1000 1000 1000])
ind_var = linspace(0,max(Distance_corr),10);
zeroes = linspace(0,0,10);
for p = 1:numel(ind_var)
    dep_var(p) = ind_var(p);
end
plot(ind_var, dep_var, 'k')
scatter(Distance, delay_corrected, 'k')
%scatter(Distance_corr, delay_corrected(1,:), 'k')
hold on
%scatter(Distance_corr, delay_corrected(2,:), 'm')
%scatter(Distance_corr, delay_corrected(3,:))
plot(ind_var, zeroes, 'k-.')
text(Distance+0.01, delay_corrected(1,:), sta);
xlabel('Distance (m)')
ylabel('Time Values with Reference to Closest Station (s)')
filename = sprintf('%s_timeVsDist_%1.4f_%1.4f.png',eq(earthquake_number).name,fil(1),fil(2));
filename_wPath = fullfile(directory,filename);
hgexport(g, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

r = figure;
set(r, 'Position', [1000 1000 1000 1000])
zeroes = linspace(0,0,10);
hold on
scatter(delay_corrected, Q, 'k')
%scatter(delay_corrected(1,:), Q, 'k')
text(delay_corrected(1,:), Q, sta);
hold on
%scatter(delay_corrected(2,:), Q, 'm')
%scatter(delay_corrected(3,:), Q)
xlabel('Time Delay (s)')
ylabel('Apparent Q')
filename = sprintf('%s_QvsDelay_%1.4f_%1.4f.png',eq(earthquake_number).name,fil(1),fil(2));
filename_wPath = fullfile(directory,filename);
hgexport(r, filename_wPath, hgexport('factorystyle'), 'Format', 'png');








%%
close all
visualization(w_clean_sort, Q,eq(earthquake_number).name, fil, eq(earthquake_number).az, earthquake_number, delay_corrected(1,:));



%%


%%
% w_clean_tp = taper(w_clean, 0.2);
% 
% f = filterobject('b', [0.8 25], 2);
% w_clean_filt = filtfilt(f, w_clean_tp);
% 
% s = spectralobject(1024, [], 10, []);
% figure(2)
% specgram_iceweb(s, w_clean, 0.75)