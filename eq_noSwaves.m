%function Plutons_waveform_analysis(earthquake_number, fil);
%load the waveform data from Antelope database dbmerged
clc
clear all
ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
scnl = scnlobject('*', '*', 'PL');

eq = struct('orid', 2183, 'snum', datenum(2011, 5, 13, 3, 34, 50), 'enum', datenum(2011, 5, 13, 3, 35, 30), 'lat', -21.8353, 'lon', -67.2705, 'depth', 179, 'evtime', datenum(2011, 5, 13, 3, 34, 27));

w_raw = waveform(ds, scnl, eq.snum, eq.enum);
%disp(w(1)) %check to make sure that the above worked



% %% Plot waveforms to see if the processing worked
% figure(1)
% len = numel(w_raw)
% for i=1:3
%     %subplot(len, 1, i)
%     figure
%     plot(w(i))
% end

%%
%w_clean = waveform_clean(w_raw);
fil=[0 0];
w_clean = waveform_clean(w_raw);


%get information for these waveforms
stations = [get(w_clean, 'station')];
stations;

%%
channels = [get(w_clean, 'channel')];
calibs = [get(w_clean, 'calib')];
len = numel(w_clean);
max_vals=[];
min_vals=[];
for i=1:len
    max_vals(i) = max(w_clean(i));
    min_vals(i) = min(w_clean(i));
end
absMax = max(max_vals);
absMin = min(min_vals);


siteStruct = loadSiteTable('/raid/data/antelope/databases/PLUTONS/dbmerged');
siteSta = siteStruct.sta;
staStruct = struct();
SECS2DAY = 60 * 60 * 24;
for i=1:len
    for k = 1:numel(siteSta) 
        if strcmp(stations{i}, siteSta{k})
            staStruct(i).sta = get(w_clean(i),'station');
            staStruct(i).lat = siteStruct.lat(k);
            staStruct(i).lon = siteStruct.lon(k);
            staStruct(i).dist = distance(eq.lat, eq.lon, siteStruct.lat(k), siteStruct.lon(k));
            [times, phasenames] = arrtimes(staStruct(i).dist, eq.depth);
            %staStruct(i).timeDiff = times[];
            staStruct(i);
            w_clean(i) = addfield(w_clean(i), 'LAT', staStruct(i).lat);
            w_clean(i) = addfield(w_clean(i), 'LON', staStruct(i).lon);
            w_clean(i) = addfield(w_clean(i), 'DISTANCE', distance(eq.lat, eq.lon, siteStruct.lat(k), siteStruct.lon(k)) );
            w_clean(i) = addfield(w_clean(i), 'EX_ARR_TIME', (eq.evtime + times(1)/SECS2DAY));
        end
    end
end

[Y,order] = sort([staStruct.dist]);


w_clean_sort = w_clean(order);


for cc=1:numel(w_clean_sort)
    get(w_clean_sort(cc),'station');
end


% ----------- Split up waveform object because earthquake shows up at too
% many stations ------------------ %
w_clean_sort1 = w_clean_sort(1:27);
w_clean_sort2 = w_clean_sort(28:54); 
w_clean_sort3 = w_clean_sort(55:78);

sta_array = [22 23 24 31 32 33 34 35 36 40 41 42 49 50 51 52 53 54 55 56 57 67 68 69 70 71 72];

w_clean_sort_paper1 = w_clean_sort(sta_array); %HH? - all 3 components
w_clean_sort_paper2 = w_clean_sort_paper1;
w_clean_sort_paper2(2:3:27) = [];


% figure(3)
% for i=1:7
%     subplot_tight(7, 1, i)
%     plot(w_clean(i), 'xunit', 'date')
%     ylim([absMin absMax])
%     ylabel('nm/s')
% end

close all;
edit_mulplt(w_clean_sort1, 0, absMax, absMin, eq.orid, fil, '1-27')
edit_mulplt(w_clean_sort2, 0, absMax, absMin, eq.orid, fil, '28-54')
edit_mulplt(w_clean_sort3, 0, absMax, absMin, eq.orid, fil, '55-78')

edit_mulplt(w_clean_sort_paper2, 0, absMax, absMin, eq.orid, fil, 'selected')

%%
lat = [get(w_clean_sort, 'LAT')];
lon = [get(w_clean_sort, 'LON')];
station = [get(w_clean_sort, 'station')];

latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.5]; %[western_limit eastern_limit]


%----------- First Figure - Outer Circle with Small and Big Q-Values
%Flipped, Filling outer circle ---------------

h = figure;
set(h, 'Position', [1000 1000 1000 1000])
worldmap(latlim, lonlim);
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
hold on

scatterm(eq.lat, eq.lon, '*', 'k')

    for i = 1:numel(station)
        ismember(i,sta_array)
        if ismember(i,sta_array) == 0
            scatterm(lat(i), lon(i), 'k')
        elseif ismember(i,sta_array) == 1
            scatterm(lat(i), lon(i), 'k', 'filled')
        end
        c = cellstr(station(i));
        textm(lat(i)+0.010, lon(i)+0.010, c);
    end
scatterm(-22.27, -67.18, 100, '^', 'k')
    
hold off

filename = sprintf('eq_loc_%4.0f.png',eq.orid);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');