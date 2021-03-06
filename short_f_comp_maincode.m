%function Plutons_waveform_analysis(earthquake_number, fil);
%load the waveform data from Antelope database dbmerged
clc
clear all
close all
startup_seismiclab
addpath(genpath('/raid/apps/src/GEOTOOLS/matlab_util'))
addpath('data_func_matTaup/')
addpath('lldistkm/')
addpath('latlonutm/Codes/matlab/')
addpath('readhgt/')
addpath('pathdist_v4')
addpath('scalebar_v3/scalebar')
NEWGISMODIR=fullfile('/raid/apps/matlab/toolbox/GISMO/startup_GISMO.m');
rmpath(genpath(NEWGISMODIR))
addpath('/raid/apps/src/gismotools/GISMO')
startup_GISMO

ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
earthquake_number = 6;
scnl = scnlobject('*', 'HHZ', 'PL');

utu_lat = -22.27;
utu_lon = -67.18;

%ESZ1
eq(1) = struct('name', 'ESZ1', 'snum', datenum(2011, 5, 19, 15, 02, 02), 'enum', datenum(2011, 5, 19, 15, 02, 31), 'lat', 39.21, 'lon', 14.96, 'depth', 292, 'mag', 4.7, 'evtime', datenum(2011, 5, 19, 14, 50, 54), 'freq', 1/0.9, 'az', 231, 'aoi', 13);
%KTSZ1
eq(2) = struct('name', 'KTSZ1', 'snum', datenum(2011, 4, 18, 13, 16, 17), 'enum', datenum(2011, 4, 18, 13, 16, 27), 'lat', -34.34, 'lon', 179.87, 'depth', 86, 'mag', 6.6, 'evtime', datenum(2011, 4, 18, 13, 3, 2), 'freq', 1/1.4, 'az', 50, 'aoi', 13);
%KTSZ2
eq(3) = struct('name', 'KTSZ2', 'snum', datenum(2011, 7, 29, 7, 55, 10), 'enum', datenum(2011, 7, 29, 7, 55, 25), 'lat', -23.78, 'lon', 179.76, 'depth', 523, 'mag', 6.7, 'evtime', datenum(2011, 7, 29, 7, 42, 22), 'freq', 1/2.1, 'az', 59, 'aoi', 13);
%KTSZ3
eq(4) = struct('name', 'KTSZ3', 'snum', datenum(2011, 9, 15, 19, 43, 45), 'enum', datenum(2011, 9, 15, 19, 44, 1), 'lat', -21.61, 'lon', -179.53, 'depth', 645, 'mag', 7.3, 'evtime', datenum(2011, 9, 15, 19, 31, 4), 'freq', 1/3.8, 'az', 61, 'aoi', 13);
%KTSZ4
eq(5) = struct('name', 'KTSZ4', 'snum', datenum(2012, 1, 24, 1, 4, 50), 'enum', datenum(2012, 1, 24, 1, 5, 5), 'lat', -24.98, 'lon', 178.52, 'depth', 580, 'mag', 6.3, 'evtime', datenum(2012, 1, 24, 0, 52, 5), 'freq', 1/3.4, 'az', 57, 'aoi', 13);
%JSZ1
eq(6) = struct('name', 'JSZ1', 'snum', datenum(2010, 11, 30, 3, 43, 35), 'enum', datenum(2010, 11, 30, 3, 43, 45), 'lat', 28.36, 'lon', 139.15, 'depth', 486, 'mag', 6.8, 'evtime', datenum(2010, 11, 30, 3, 24, 41), 'freq', 1/1.3, 'az', 290, 'aoi', 4);
%JSZ2
eq(7) = struct('name', 'JSZ2', 'snum', datenum(2011, 1, 12, 21, 51, 45), 'enum', datenum(2011, 1, 12, 21, 51, 55), 'lat', 26.98, 'lon', 139.87, 'depth', 527, 'mag', 6.4, 'evtime', datenum(2011, 1, 12, 21, 32, 55), 'freq', 1/1.5, 'az', 286, 'aoi', 4);
%JSZ3
eq(8) = struct('name', 'JSZ3', 'snum', datenum(2011, 5, 10, 15, 44, 52), 'enum', datenum(2011, 5, 10, 15, 45, 2), 'lat', 43.29, 'lon', 130.94, 'depth', 544, 'mag', 5.4, 'evtime', datenum(2011, 5, 10, 15, 26, 5), 'freq', 1/1.3, 'az', 329, 'aoi', 3);
%JSZ4
eq(9) = struct('name', 'JSZ4', 'snum', datenum(2011, 10, 4, 1, 56, 25), 'enum', datenum(2011, 10, 4, 1, 56, 40), 'lat', 26.77, 'lon', 140.43, 'depth', 455, 'mag', 5.6, 'evtime', datenum(2011, 10, 4, 1, 37, 29), 'freq', 1/1.3, 'az', 286, 'aoi', 4);
%SSSZ1
eq(10) = struct('name', 'SSSZ1', 'snum', datenum(2011, 1, 20, 3, 52, 35), 'enum', datenum(2011, 1, 20, 3, 52, 55), 'lat', -59.94, 'lon', -27.48, 'depth', 129, 'mag', 5.2, 'evtime', datenum(2011, 1, 20, 3, 44, 26), 'freq', 1/1.15, 'az', 154, 'aoi', 24);
%SSSZ2
eq(11) = struct('name', 'SSSZ2', 'snum', datenum(2011, 1, 23, 23, 0, 59), 'enum', datenum(2011, 1, 23, 23, 1, 15), 'lat', -56.47, 'lon', -26.96, 'depth', 124, 'mag', 5.3, 'evtime', datenum(2011, 1, 23, 22, 53, 2), 'freq', 1/0.75, 'az', 150, 'aoi', 24);
%SSSZ3
eq(12) = struct('name', 'SSSZ3', 'snum', datenum(2011, 6, 19, 8, 45, 40), 'enum', datenum(2011, 6, 19, 8, 45, 55), 'lat', -56.05, 'lon', -27.42, 'depth', 119, 'mag', 5.2, 'evtime', datenum(2011, 6, 19, 8, 37, 45), 'freq', 1/0.8, 'az', 150, 'aoi', 24);
%SSSZ4
eq(13) = struct('name', 'SSSZ4', 'snum', datenum(2011, 8, 21, 12, 46, 47), 'enum', datenum(2011, 8, 21, 12, 47, 5), 'lat', -56.43, 'lon', -27.49, 'depth', 130, 'mag', 5.6, 'evtime', datenum(2011, 8, 21, 12, 38, 54), 'freq', 1/1.05, 'az', 150, 'aoi', 24);
%SSSZ5
eq(14) = struct('name', 'SSSZ5', 'snum', datenum(2011, 8, 26, 7, 49, 18), 'enum', datenum(2011, 8, 26, 7, 49, 35), 'lat', -56.21, 'lon', -27.18, 'depth', 127, 'mag', 5.2, 'evtime', datenum(2011, 8, 26, 7, 41, 23), 'freq', 1/0.85, 'az', 150, 'aoi', 24);

w_raw = waveform(ds, scnl, eq(earthquake_number).snum, eq(earthquake_number).enum);
%disp(w(1)) %check to make sure that the above worked

% %% Plot waveforms to see if the processing worked
% figure(1)
% len = numel(w_raw)
% for i=1:3
%     %subplot(len, 1, i)
%     figure
%     plot(w(i))
% end

%w_clean = waveform_clean(w_raw);
fil=[0.375 0.75];
[tshift, C] = cross_corr(eq(earthquake_number), fil);
corr_vals = get(C, 'corr');

w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));

if earthquake_number == 4
    w_clean(12)=[]; %remove station PLLL - instrument noise KTSZ3
elseif earthquake_number == 5
    w_clean(1:7)=[]; %remove all LZ stations - KTSZ4
    w_clean(18)=[]; %remove PLSP - waveform unlike others
elseif earthquake_number == 7
    w_clean(6)=[]; %remove station PLQU - noise JSZ2
elseif earthquake_number == 8    
    w_clean(17)=[]; %remove station PLQU - noise JSZ3
    w_clean(20)=[]; %remove station PLSP - noise JSZ3
    w_clean(9)=[]; %remove station PLDK - noise JSZ3
elseif earthquake_number == 10
    w_clean(9)=[]; %remove station PLSM - noise SSSZ1
elseif earthquake_number == 11
    w_clean(6)=[]; %remove station PLQU - noise SSSZ2
elseif earthquake_number == 12
    w_clean(23)=[]; %remove station PLTP - noise SSSZ3
end

%get information for these waveforms
stations = [get(w_clean, 'station')];



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

%%
siteStruct = loadSiteTable('/raid/data/antelope/databases/PLUTONS/dbmerged');
minEl = min(siteStruct.elev);
siteSta = siteStruct.sta;
staStruct = struct();
SECS2DAY = 60 * 60 * 24;
for i=1:len
    for k = 1:numel(siteSta) 
        if strcmp(stations{i}, siteSta{k})
            staStruct(i).sta = get(w_clean(i),'station');
            staStruct(i).lat = siteStruct.lat(k);
            staStruct(i).lon = siteStruct.lon(k);
            staStruct(i).elev = siteStruct.elev(k); %elevation of station in km
            %staStruct(i).dist = taupTime([], eq(earthquake_number).depth, [], 'sta', [siteStruct.lat(k), siteStruct.lon(k)], 'evt', [eq(earthquake_number).lat, eq(earthquake_number).lon]);
            staStruct(i).dist = distance(eq(earthquake_number).lat, eq(earthquake_number).lon, siteStruct.lat(k), siteStruct.lon(k));
            [times, phasenames] = arrtimes(staStruct(i).dist, eq(earthquake_number).depth);
            %staStruct(i).timeDiff = times[];
            staStruct(i);
            w_clean(i) = addfield(w_clean(i), 'ELEV', staStruct(i).elev);
            w_clean(i) = addfield(w_clean(i), 'LAT', staStruct(i).lat);
            w_clean(i) = addfield(w_clean(i), 'LON', staStruct(i).lon);
            w_clean(i) = addfield(w_clean(i), 'DISTANCE', distance(eq(earthquake_number).lat, eq(earthquake_number).lon, siteStruct.lat(k), siteStruct.lon(k)) );
        end
    end
end

[Y,order] = sort([staStruct.dist]);


w_clean_sort = w_clean(order);

for cc=1:numel(w_clean_sort)
    get(w_clean_sort(cc),'station');
end

%Create time lag values from each station relative to the one closest to
%the earthquake


%%
% figure(3)
% for i=1:7
%     subplot_tight(7, 1, i)
%     plot(w_clean(i), 'xunit', 'date')
%     ylim([absMin absMax])
%     ylabel('nm/s')
% end

% grab a waveform
if earthquake_number == 1
    starttimes  = {'5/19/2011 15:02:02', '5/19/2011 15:02:02', '5/19/2011 15:02:03', '5/19/2011 15:02:03', '5/19/2011 15:02:03', '5/19/2011 15:02:04', '5/19/2011 15:02:04', '5/19/2011 15:02:04', '5/19/2011 15:02:04', '5/19/2011 15:02:06', '5/19/2011 15:02:06', '5/19/2011 15:02:06', '5/19/2011 15:02:07', ...
        '5/19/2011 15:02:07', '5/19/2011 15:02:08', '5/19/2011 15:02:08', '5/19/2011 15:02:08', '5/19/2011 15:02:09', '5/19/2011 15:02:09', '5/19/2011 15:02:09', '5/19/2011 15:02:10', '5/19/2011 15:02:11', '5/19/2011 15:02:11', '5/19/2011 15:02:12', '5/19/2011 15:02:13', '5/19/2011 15:02:13'};
    endtimes = {'5/19/2011 15:02:20', '5/19/2011 15:02:20', '5/19/2011 15:02:21', '5/19/2011 15:02:21', '5/19/2011 15:02:21', '5/19/2011 15:02:22', '5/19/2011 15:02:22', '5/19/2011 15:02:22', '5/19/2011 15:02:22', '5/19/2011 15:02:24', '5/19/2011 15:02:24', '5/19/2011 15:02:24', '5/19/2011 15:02:25', ...
        '5/19/2011 15:02:25', '5/19/2011 15:02:26', '5/19/2011 15:02:26', '5/19/2011 15:02:26', '5/19/2011 15:02:27', '5/19/2011 15:02:27', '5/19/2011 15:02:27', '5/19/2011 15:02:28', '5/19/2011 15:02:29', '5/19/2011 15:02:29', '5/19/2011 15:02:30', '5/19/2011 15:02:31', '5/19/2011 15:02:31'};
    w_clean_sort = subtime(w_clean_sort,starttimes,endtimes);
end

if strcmp(eq(earthquake_number).name, 'ESZ1') %change w_clean_sort for ESZ1 from 26x26 to 1x26
    for wavenum = 1:26
        w_clean_sort2(wavenum) = w_clean_sort(wavenum, wavenum);
    end
    w_clean_sort = w_clean_sort2;
end

close all;
% l = [10, 11, 13, 14, 15]; %Subset for PLLL, PLBR, PLSS, PLCM, and PLAR for KTSZ1 0.75-1.5
% for i = 1:numel(l)
%     w_clean_sort_trial(i) = w_clean_sort(l(i));
% end


%w_clean_test = w_clean_sort_trial(1);
days2secs = 60*60*24;
freq = get(w_clean_sort(1), 'freq');
data = get(w_clean_sort(1), 'data');
numel(data)
dnum(1) = datenum(get(w_clean_sort(1),'start'));
for l = 2:numel(data)
    dnum(l) = datenum(l/freq/days2secs+dnum(1));
end
if strcmp(eq(earthquake_number).name, 'KTSZ1')
    data = data(1:400);
    [ref_amp, ref_index] = nanmin(data);
    time_value_ref = dnum(ref_index); %reference time of minimum of first station
    start_time_ref = dnum(1); %start time of waveforms
    diff_time_ref = time_value_ref - start_time_ref %difference between start time of series and phase time
elseif strcmp(eq(earthquake_number).name, 'KTSZ3') || strcmp(eq(earthquake_number).name, 'KTSZ4') || strcmp(eq(earthquake_number).name, 'SSSZ2')
    if strcmp(eq(earthquake_number).name, 'KTSZ4') && fil(1) == 0.1875 && fil(2) == 3.000 
        data = data(300:500);
    elseif strcmp(eq(earthquake_number).name, 'KTSZ4') && fil(1) == 0.3750 && fil(2) == 1.500
        data = data(1:400);
    else
        data = data(1:500);
    end
    if strcmp(eq(earthquake_number).name, 'KTSZ3')
        [ref_amp, ref_index] = nanmax(data)
    else
        [ref_amp, ref_index] = nanmin(data)
        if strcmp(eq(earthquake_number).name, 'KTSZ4') && fil(1) == 0.1875 && fil(2) == 3.000
            ref_index = ref_index+300;
        end
    end
    time_value_ref = dnum(ref_index); %reference time of minimum of first station
    start_time_ref = dnum(1); %start time of waveforms
    diff_time_ref = time_value_ref - start_time_ref %difference between start time of series and phase time
elseif strcmp(eq(earthquake_number).name, 'KTSZ2')
    data = data(1:650);
    [ref_amp, ref_index] = nanmin(data);
    time_value_ref = dnum(ref_index) %reference time of minimum of first station
    start_time_ref = dnum(1) %start time of waveforms
    diff_time_ref = time_value_ref - start_time_ref %difference between start time of series and phase time
elseif strcmp(eq(earthquake_number).name, 'JSZ1') || strcmp(eq(earthquake_number).name, 'JSZ3')
    data = data(1:700);
    [ref_amp, ref_index] = nanmin(data);
    time_value_ref = dnum(ref_index); %reference time of minimum of first station
    start_time_ref = dnum(1); %start time of waveforms
    diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
elseif strcmp(eq(earthquake_number).name, 'JSZ2')
    data = data(1:600);
    [ref_amp, ref_index] = nanmin(data);
    time_value_ref = dnum(ref_index); %reference time of minimum of first station
    start_time_ref = dnum(1); %start time of waveforms
    diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
elseif strcmp(eq(earthquake_number).name, 'JSZ4')
    if fil(1) == 0.1875 && fil(2) == 3.000
        nums = ceil(numel(data));
        data = data((nums/3):600);
        [ref_amp, ref_index] = nanmax(data);
        ref_index = ref_index+(nums/3);
        time_value_ref = dnum(ref_index); %reference time of minimum of first station
        start_time_ref = dnum(1); %start time of waveforms
        diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
    else
        data = data(1:600);
        [ref_amp, ref_index] = nanmax(data);
        time_value_ref = dnum(ref_index); %reference time of minimum of first station
        start_time_ref = dnum(1); %start time of waveforms
        diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
    end
elseif strcmp(eq(earthquake_number).name, 'SSSZ1')
    if (fil(1)==0.7500 && fil(2)==1.500) || (fil(1)==0.75 && fil(2) == 3.000)
        data = data(1:900);
        [ref_amp, ref_index] = nanmax(data);
        time_value_ref = dnum(ref_index); %reference time of minimum of first station
        start_time_ref = dnum(1); %start time of waveforms
        diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
    else
        data = data(1:1000);
        [ref_amp, ref_index] = nanmax(data);
        time_value_ref = dnum(ref_index); %reference time of minimum of first station
        start_time_ref = dnum(1); %start time of waveforms
        diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
    end
elseif strcmp(eq(earthquake_number).name, 'SSSZ3')
    if fil(1)== 0.1875 && fil(2)==0.75
        data = data(100:400);
        [ref_amp, ref_index] = nanmax(data);
        time_value_ref = dnum(ref_index+100); %reference time of minimum of first station
        start_time_ref = dnum(1); %start time of waveforms
        diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
    else
        data = data(1:400);
        [ref_amp, ref_index] = nanmax(data);
        time_value_ref = dnum(ref_index); %reference time of minimum of first station
        start_time_ref = dnum(1); %start time of waveforms
        diff_time_ref = time_value_ref - start_time_ref; %difference between start time of series and phase time
    end
end

ref_index


tshift_days = [];
tshift_time_days = [];
for i=1:numel(tshift)
    tshift_days(i) = tshift(i)/days2secs;
    tshift_time_days(i) = tshift_days(i) + time_value_ref;
end


%%
if numel(tshift_time_days) == numel(w_clean_sort)
    [index_values, time_values, m_values] = edit_mulplt_eqSpecific(w_clean_sort, 0, absMax, absMin, eq(earthquake_number).name, fil, tshift_time_days);
end