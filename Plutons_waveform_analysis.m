%function Plutons_waveform_analysis(earthquake_number, fil);
%load the waveform data from Antelope database dbmerged
clc
clear all
close all
startup_seismiclab
addpath(genpath('/raid/apps/src/GEOTOOLS/matlab_util'))
addpath('data_func_matTaup/')
addpath('latlonutm/Codes/matlab/')
addpath('readhgt/')
NEWGISMODIR=fullfile('/raid/apps/matlab/toolbox/GISMO/startup_GISMO.m');
rmpath(genpath(NEWGISMODIR))
addpath('/raid/apps/src/gismotools/GISMO')
startup_GISMO

ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
earthquake_number = 9;
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
eq(8) = struct('name', 'JSZ3', 'snum', datenum(2011, 5, 10, 15, 44, 52), 'enum', datenum(2011, 5, 10, 15, 45, 2), 'lat', 43.29, 'lon', 130.94, 'depth', 544, 'mag', 5.4, 'evtime', datenum(2011, 5, 10, 15, 26, 5), 'freq', 1/1.3, 'az', 286, 'aoi', 4);
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
tshift = cross_corr(eq(earthquake_number), fil);

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
stations;



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
            w_clean(i) = addfield(w_clean(i), 'EX_ARR_TIME', (eq(earthquake_number).evtime + times(1)/SECS2DAY)+((((staStruct(i).elev-minEl)/cosd(eq(earthquake_number).aoi)))*(4.1/SECS2DAY))); %adding in the travel time from sea level to the station elevation using AOI
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
ex_arr_time = [get(w_clean_sort, 'EX_ARR_TIME')];

tlag = [];
for values = 1:numel(ex_arr_time)
    if values == 1
        seclag = 0;
    else
        seclag = (ex_arr_time(values)-ex_arr_time(1))*SECS2DAY;
    end
    tlag(values) = seclag;
end

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
    data = data(1:300);
    [ref_amp, ref_index] = nanmax(data);
    time_value_ref = dnum(ref_index); %reference time of minimum of first station
    start_time_ref = dnum(1); %start time of waveforms
    diff_time_ref = time_value_ref - start_time_ref %difference between start time of series and phase time
elseif strcmp(eq(earthquake_number).name, 'KTSZ3') || strcmp(eq(earthquake_number).name, 'KTSZ4') || strcmp(eq(earthquake_number).name, 'SSSZ2')
    if strcmp(eq(earthquake_number).name, 'KTSZ4') && fil(1) == 0.1875 && fil(2) == 3.000 
        data = data(300:500);
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
%time_values in dnum
%%

time_vals_ref = [];
if earthquake_number == 9
    ref_time = time_values(3);
else
    ref_time = time_values(1);
end
for values = 1:numel(time_values)
    timelag = (time_values(values)-ref_time)*SECS2DAY;
    time_vals_ref(values) = timelag;
end
slowness = tlag-time_vals_ref;
%time_values(17)-time_values(18)

%%
% w_clean_tp = taper(w_clean_sort_trial, 0.2);
% 
% s = spectralobject(200, [], 2, []);
% figure
% edit_specgram_iceweb(s, w_clean_tp, 0.75)
% data = [get(w_clean_tp,'data')];
% numel(data)



% %%
% for wavnum = 1:nwaveforms
%         data = get(w_clean_sort(wavnum),'data');
%         stn = get(w_clean_sort(wavnum),'station');
%         %[m,I] = corr_arrival_picks(data, stn, eq(earthquake_number).name, fil);
%         m = m_values(wavnum);
%         dnum=datenum(w_clean_sort(wavnum));
%         %time_value = dnum(I)
%         %data = data(index_values(wavnum)-150:index_values(wavnum)+150);
%         [maxValue,indexMax] = max(abs(fft(data)));
%         x = linspace(0,numel(data),numel(data));
%         frequency = indexMax * 150 / numel(data);
%         figure
%         r = abs(fft(data))/numel(data);
%         plot(x,r)
%         %plot(x, abs(fft(data))/numel(data))
%         yl(1) = min(abs(fft(data))/numel(data));
%         yl(2) = max(abs(fft(data))/numel(data));
%         hold on
%         line([frequency, frequency], [yl(1), yl(2)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
%         xlim([0 10])
%         xlabel('Frequency (Hz)')
%         ylabel('Strength')
%         title(sprintf('%s',stn))
%         w_clean_sort(wavnum) = addfield(w_clean_sort(wavnum), 'AMP_ABS', abs(m));
%         w_clean_sort(wavnum) = addfield(w_clean_sort(wavnum), 'sig_freq', frequency);
% end

%%
%[N,M] = size(w_clean_sort);
close all
nwaveforms = numel(w_clean_sort);
for i = 1:nwaveforms
    m = m_values(i);
    Fn = get(w_clean_sort(i),'NYQ'); %from https://code.google.com/p/gismotools/source/browse/trunk/GISMO/contributed/fft_tools/%2Bwf_fft/compute.m?r=321
    stn = get(w_clean_sort(i),'station');
    x = get(w_clean_sort(i),'data');
    if strcmp(eq(earthquake_number).name, 'JSZ1')
        if fil(1)==0.1875 && fil(2)==0.75
            if strcmp(stn, 'PLAR')
                x = x(index_values(i)-200:index_values(i)+197);
            end
        end
    else
        x = x(index_values(i)-200:index_values(i)+200);
    end
    NFFT=2.^(ceil(log(length(x))/log(2)));  % Next highest power of 2
    FFTX=fft(x,NFFT);                       % Take fft, padding with zeros.
    NumUniquePts = ceil((NFFT+1)/2);
    FFTX=FFTX(1:NumUniquePts);              % throw out neg frequencies
    MX=abs(FFTX);                           % Take magnitude of X
    MX=MX*2;                                % Multiply by 2 
    MX=MX/length(x);                        
    %PX=phase(FFTX);                           % Take magnitude of X
    f=(0:NumUniquePts-1)*2/NFFT;            
    f=f*Fn;
    figure
    plot(f,MX)
    hold on
    xlim([0 5])
    title(sprintf('%s',stn))
    w_clean_sort(i) = addfield(w_clean_sort(i),'FFT_FREQ',f');
    w_clean_sort(i) = addfield(w_clean_sort(i),'FFT_AMP',MX);
    w_clean_sort(i) = addfield(w_clean_sort(i), 'AMP_ABS', abs(m));
    %w(i) = addfield(w(i),'FFT_PHASE',PX);
    a = find(MX == max(MX));
    line([f(a), f(a)], [min(MX), max(MX)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
    w_clean_sort(i) = addfield(w_clean_sort(i),'FFT_DOM',f(a));
end;

sig_freq = [get(w_clean_sort, 'FFT_DOM')];
sta = [get(w_clean_sort, 'station')];

[stas,Q] = q_calc(w_clean_sort, 5.15, 20, sig_freq);

directory = sprintf('/home/a/akfarrell/Uturuncu/%s/text', eq(earthquake_number).name);
filename2 = sprintf('%s_output_diffFreq_%1.4f_%1.4f.txt',eq(earthquake_number).name,fil(1),fil(2));
dif=fopen(fullfile(directory,filename2), 'w');
for i=1:length(stas)
    fprintf(dif,'%s %10.5f %10.3f %10.4f\n',stas(i,:),Q(i),sig_freq(i), slowness(i));
    i
end
st = fclose('all');

%%
close all

topo_data = readhgt(-23:-22,-68:-67, 'merge','interp');
lat_sta = get(w_clean_sort, 'LAT');
lon_sta = get(w_clean_sort, 'LON');
elev = get(w_clean_sort, 'ELEV')*1000;
[lon, lat] = meshgrid(topo_data.lon, topo_data.lat);
%%

[dist, az_volc_eq] = distance(utu_lat, utu_lon, eq(earthquake_number).lat, eq(earthquake_number).lon);
[dist2, az_eq_volc] = distance(eq(earthquake_number).lat, eq(earthquake_number).lon, utu_lat, utu_lon);

% ----Calc line of intersection of planar wave with surface at closest station
if earthquake_number == 9
    lat_ref = lat_sta(1);
    lon_ref = lon_sta(1);
    elev_ref = elev(1);
elseif earthquake_number == 3
    lat_ref = lat_sta(1);
    lon_ref = lon_sta(1);
    elev_ref = elev(1);
elseif earthquake_number == 10 %%%%%%%%%CHANGE!
    lat_ref = lat_sta(1);
    lon_ref = lon_sta(1);
    elev_ref = elev(1);
end
%%
%x1 = lat_ref; y1 = lon_ref;
[x1,y1] = ll2utm(lat_ref, lon_ref);
if earthquake_number >=6 && earthquake_number<=9
    m = tand(az_volc_eq-270);
elseif earthquake_number >=2 && earthquake_number<=5
    m = tand(az_volc_eq-180);
end
C = -m*x1+y1; %intercept
A = m; %slope
B = -1;

%%

[x_of_sta, y_of_sta]= ll2utm(lat_sta, lon_sta); %calculate utm coordinates of station locations, in m
x_1 = [x_of_sta(1), y_of_sta(1), elev(1)];
x_2 = [x_of_sta(1)+100, A*(x_of_sta(1)+100)+C, elev(1)];


for j = 1:numel(sta)
    x_0 = [x_of_sta(j), y_of_sta(j), elev(j)];
    d2(j) = (abs(cross((x_2-x_1), (x_1-x_0)))/abs(x_2-x_1));
end
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
