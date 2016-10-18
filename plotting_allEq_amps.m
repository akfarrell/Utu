function plotting_allEq_amps(siteStruct, siteSta, plot_variable, num_eqs, id, num_averaged, scale_factor)
%run from test_visuals_amps.m
h = figure;
latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.6]; %[western_limit eastern_limit]
set(h, 'Position', [1000 1000 1100 1000])

worldmap(latlim, lonlim);
% --- Define station characteristics ----%

setm(gca, 'fontsize',11);
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
grey = rgb('Grey');
thistle = rgb('Thistle');
colors = {'k', 'm', 'c', 'g', 'b', 'r', 'y', grey, thistle}; %need to add more when I have more events

% ----- Make Plot ------- %            
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
northarrow('latitude', -21.9, 'longitude', -67.6, 'scaleratio', 1/20);
scalebar('location', 'sw', 'FontSize', 11)
hold on
plot_variable = plot_variable * 2;
colormap(jet)
colormap(flipud(colormap))
set(gca,'CLim',[min(plot_variable),max(plot_variable)]);
for i = 1:numel(sta_name)
    if isnan(plot_variable(i))
        continue
    elseif strcmp(sta_name(i), 'PLRV') || strcmp(sta_name(i), 'PLMD') %remove stations with few samples
        continue
    else
        scatterm(lat_sta(i), lon_sta(i), plot_variable(i), plot_variable(i), 'filled')
        textm(lat_sta(i)+0.025, lon_sta(i)-0.02, sta_name(i), 'FontSize', 12)
        sta_name
        plot_variable
        if exist('num_averaged', 'var')
            for count = 1:max(num_averaged)
                if num_averaged(i) == num_vals(count)
                    textm(lat_sta(i)-0.01, lon_sta(i)+0.02, num2str(num_averaged(i)), 'FontSize', 12)
                end
            end
        end
    end
end
tick_labels = (60:20:220)./(scale_factor*2);
for i=1:numel(tick_labels)
    tick_vals(i) = {sprintf('%1.2f',tick_labels(i))}
end
tick_vals
c=colorbar('eastoutside', 'TickLabels',tick_vals)
c.Label.String='Average Normalized Amplitude'
c.Label.FontSize=11;
%title(sprintf('%s %d',id, num_eqs))
directory = sprintf('/home/a/akfarrell/Uturuncu/Synth_data/');
if exist('num_averaged', 'var')
    filename = sprintf('averaged_amps_%d_%s_nums.tiff',num_eqs, id);
else
    filename = sprintf('averaged_amps_%d_%s.tiff',num_eqs, id);
end
filename_wPath = fullfile(directory,filename);
hold off
%hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'pdf');
fp = fillPage(h, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
print(h, '-dtiff', '-r400',filename_wPath);