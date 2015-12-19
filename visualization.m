function visualization(w_clean_sort, Q,name,fil,az,eq_number, slowness)
% visualization  Visualize Q data in several maps of Uturuncu area
%   visualization(w_clean_sort, Q) visualizes Q data with stations and the
%   volcano
%
% visualization(w_clean_sort, Q)
%
%   Inputs:
%       w_clean_sort - waveform object with latitude and longitude of
%       stations
%       Q - double containing Q values in order of sorted station
%       name - name of the earthquake, to name the figures
%       fil - 2x1 filter object, to name the figures
%       az - azimuth data for earthquake
%       eq_number - earthquake number, to determine where to put the
%       azimuth origin
%
%   Author: Alexandra Farrell 2015/6/1 using mapping toolbox and border
%   data from http://www.gadm.org/


lat = [get(w_clean_sort, 'LAT')];
lon = [get(w_clean_sort, 'LON')];
station = [get(w_clean_sort, 'station')];

latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.5]; %[western_limit eastern_limit]


%----------- First Figure - Outer Circle with Small and Big Q-Values
%Flipped, Filling outer circle ---------------

% h = figure;
% set(h, 'Position', [1000 1000 1000 1000])
% worldmap(latlim, lonlim);
% borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
% arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
% geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% hold on
% invQ = [];
% maxval = max(Q);
% min2 = sort(Q(:));
% minval = min2(2);
% total = maxval+minval;
% 
%     for i = 1:numel(Q)
%         invQ(i) = 4*(total - Q(i));
%         scatterm(lat(i), lon(i), maxval*4, 'k')
%         scatterm(lat(i),lon(i),invQ(i), 'k', 'filled')
%         dx = maxval*4/10000;
%         dy = maxval*4/10000;
%         c = cellstr(station(i));
%         textm(min(lat(i)+dx, lat(i)+0.020), min(lon(i)+dy, lon(i)+0.020), c);
%     end
%     scatterm(-22.27, -67.18, 100, '^', 'k')
  hypotenuse = 0.07;     
  if eq_number >=6 && eq_number<=9
    lat_az = -22.05;
    lon_az = -67.60;
    u = hypotenuse*sind(360-az-90); %vertical
    v = hypotenuse*cosd(360-az-90); %horizontal
  elseif eq_number >=2 && eq_number<=5
    lat_az = -22.7;
    lon_az = -67.60;
    u = hypotenuse*cosd(az); %vertical
    v = hypotenuse*sind(az); %horizontal
  elseif eq_number >=10 && eq_number<=14
    lat_az = -22.7;
    lon_az = -66.95;
    u = hypotenuse*sind(az-90); %vertical
    v = hypotenuse*-cosd(az-90); %horizontal
  end  


% 
% if eq_number == 1 %ESZ
%         lat_az = -21.80;
%         lon_az = -66.55;
%     elseif eq_number > 1 && eq_number < 6 %KTSZ
%         lat_az = -22.70;
%         lon_az = -67.70;
%     elseif eq_number > 5 && eq_number < 10 %JSZ
%         lat_az = -21.85;
%         lon_az = -67.60;
%     else %SSSZ
%         lat_az = -22.70;
%         lon_az = -66.55;
%     end

%     u = hypotenuse*sind(360-az-90); %vertical
%     v = hypotenuse*cosd(360-az-90); %horizontal
%     quiverm(lat_az, lon_az,u, v, 'k')
%     title(sprintf('%s',name))
% hold off
% directory = sprintf('/home/a/akfarrell/Uturuncu/%s/figures', name);
% filename = sprintf('%s_q_values1_%1.4f_%1.4f.png',name,fil(1),fil(2));
% filename_wPath = fullfile(directory,filename);
% hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
% 
% 
% %----------- Second Figure - Small and Big Q-Values Flipped -----------
% b = figure;
% set(b, 'Position', [1000 1000 1000 1000])
% worldmap(latlim, lonlim);
% %s = dcwdata('SOAMAFR', -21.75, -67.75, 'PO', 'line');
% borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
% arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
% geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% hold on
% for i = 1:numel(Q)
%         invQ(i) = 4*(total - Q(i));
%         scatterm(lat(i),lon(i),invQ(i), 'k', 'filled')
%         dx = invQ(i)/10000;
%         dy = invQ(i)/10000;
%         c = cellstr(station(i));
%         textm(min(lat(i)+dx, lat(i)+0.020), min(lon(i)+dy, lon(i)+0.020), c);
% end
% scatterm(-22.27, -67.18, 100, '^', 'k')
% quiverm(lat_az, lon_az,u, v, 'k')
% title(sprintf('%s',name))
% hold off
% directory = sprintf('/home/a/akfarrell/Uturuncu/%s/figures', name);
% filename = sprintf('%s_q_values2_%1.4f_%1.4f.png',name,fil(1),fil(2));
% filename_wPath = fullfile(directory,filename);
% hgexport(b, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
% 
% 
% %----------- Third Figure - Contour Map -----------
% % r = figure;
% % set(r, 'Position', [1000 1000 1000 1000])
% % worldmap(latlim, lonlim);
% % %s = dcwdata('SOAMAFR', -21.75, -67.75, 'PO', 'line');
% % borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
% % arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
% % geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% % geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% % hold on
% % for i = 1:numel(Q)
% %     scatterm(lat(i),lon(i), 'k', 'filled')
% %     c = cellstr(station(i));
% %     textm(lat(i)+0.020, lon(i)+0.020, c);
% % end
% % contourfm(lat, lon, Q)
% % contourcmap('jet', 'Colorbar', 'on', 'Location', 'horizontal', 'TitleString', 'Contour Intervals')
% % scatterm(-22.27, -67.18, 100, '^', 'k')
% % quiverm(lat_az, lon_az,u, v, 'k')
% % title(sprintf('%s',name))
% % hold off
% % filename = sprintf('%s_q_values_contour_%1.4f_%1.4f.png',name,fil(1),fil(2));
% % hgexport(r, filename, hgexport('factorystyle'), 'Format', 'png');
% 
% %------------- Fourth Figure - Colored Dots ----------
p=figure; hold on;
set(p, 'Position', [1000 1000 1000 1000])
worldmap(latlim, lonlim);
%s = dcwdata('SOAMAFR', -21.75, -67.75, 'PO', 'line');
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
colormap(flipud(colormap))
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
scatterm(lat,lon,20^2,Q,'filled');
textm(lat-0.04,lon-0.025,station); colorbar
scatterm(-22.27, -67.18, 100, '^', 'k')
quiverm(lat_az, lon_az,u, v, 'k')
title(sprintf('%s Q Values',name))
hold off
directory = sprintf('/home/a/akfarrell/Uturuncu/%s/figures', name);
filename = sprintf('%s_q_values_coloredDots_%1.4f_%1.4f.png',name,fil(1),fil(2));
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%---------------- Fifth Figure - Slowness Dots -------------
h=figure; hold on;
set(h, 'Position', [1000 1000 1000 1000])
worldmap(latlim, lonlim);
%s = dcwdata('SOAMAFR', -21.75, -67.75, 'PO', 'line');
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
scatterm(lat,lon,20^2,slowness,'filled');
textm(lat+0.006,lon+0.025,station); colorbar
scatterm(-22.27, -67.18, 100, '^', 'k')
quiverm(lat_az, lon_az,u, v, 'k')
title(sprintf('%s Slowness',name))
hold off
directory = sprintf('/home/a/akfarrell/Uturuncu/%s/figures', name);
filename = sprintf('%s_slowness1_%1.4f_%1.4f.png',name,fil(1),fil(2));
filename_wPath = fullfile(directory,filename);
hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
