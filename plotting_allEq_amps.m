function plotting_allEq_amps(siteStruct, siteSta, plot_variable, num_eqs, id, num_averaged)
h = figure;
latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.5]; %[western_limit eastern_limit]
set(h, 'Position', [1000 1000 1000 1000])
worldmap(latlim, lonlim);
% --- Define station characteristics ----%

siteStaRaw = siteStruct.sta(9:numel(siteStruct.sta)); %remove Lazufre stations
siteLatRaw = siteStruct.lat(9:numel(siteStruct.sta));
siteLonRaw = siteStruct.lon(9:numel(siteStruct.sta));
for i=1:numel(siteStaRaw)
    for k = 1:numel(siteSta) 
        if strcmp(siteStaRaw{i}, siteSta{k})
            lat_sta(k) = siteLatRaw(i);
            lon_sta(k) = siteLonRaw(i);
            sta_name(k) = siteStaRaw(i);
        end
    end
end

if exist('num_averaged', 'var')
    num_vals = 1:max(num_averaged);
end
colors = {'k', 'm', 'c', 'g', 'b', 'r', 'y'}; %need to add more when I have more events

% ----- Make Plot ------- %            
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
hold on
for i = 1:numel(sta_name)
    if isnan(plot_variable(i))
        continue
    else
        scatterm(lat_sta(i), lon_sta(i), plot_variable(i), 'k', 'filled')
        textm(lat_sta(i)+0.02, lon_sta(i), sta_name(i))
        if exist('num_averaged', 'var')
            for count = 1:max(num_averaged)
                if num_averaged(i) == num_vals(count)
                    scatterm(lat_sta(i), lon_sta(i), plot_variable(i), colors{count}, 'filled')
                end
            end
        end
    end
end
title(sprintf('%s %d',id, num_eqs))
directory = sprintf('/home/a/akfarrell/Uturuncu/Synth_data/');
if exist('num_averaged', 'var')
    filename = sprintf('averaged_amps_%d_%s_nums',num_eqs, id);
else
    filename = sprintf('averaged_amps_%d_%s',num_eqs, id);
end
filename_wPath = fullfile(directory,filename);
hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');