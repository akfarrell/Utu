close all;

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
scnl = scnlobject('*', '*', 'PL');

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
%%
%disp(w_raw(1)) %check to make sure that the above worked

% %% Plot waveforms to see if the processing worked
% figure(1)
% len = numel(w_raw)
% for i=1:3
%     %subplot(len, 1, i)
%     figure
%     plot(w(i))
% end

%w_clean = waveform_clean(w_raw);
fil=[0.375 1.5];

w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));

stations = [get(w_clean, 'station')];
siteStruct = loadSiteTable('/raid/data/antelope/databases/PLUTONS/dbmerged');
minEl = min(siteStruct.elev);
siteSta = siteStruct.sta;
staStruct = struct();
SECS2DAY = 60 * 60 * 24;
len = numel(w_clean)
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
            %w_clean(i) = addfield(w_clean(i), 'EX_ARR_TIME', (eq(earthquake_number).evtime + times(1)/SECS2DAY)+((((staStruct(i).elev-minEl)/cosd(eq(earthquake_number).aoi)))*(4.1/SECS2DAY))); %adding in the travel time from sea level to the station elevation using AOI
        end
    end
end

[Y,order] = sort([staStruct.dist]);


w_clean_sort = w_clean(order);

%%
SECSPERDAY = 60 * 60 * 24;
%index_values = load('index_values.mat');

index_valuez = [571 549 534 605 566 542 591 611 611 613 574 549 592 579 603 639 615 662 618 623 600];
n = 3; 
index = reshape(repmat(index_valuez(:).',n,1),1,[]);

%chan = get(w_clean_sort, 'channel')
%chan(1:3:numel(chan)) %run to check that this is all HHE
%%
%[index_valuez(1)-150:1:index_valuez(1)+60]
num_vals = 1:3:numel(w_clean);
for wavnum = 1:3:numel(w_clean)
        dataE=get(w_clean_sort(wavnum),'data');
        dataN=get(w_clean_sort(wavnum+1),'data');
        dataZ=get(w_clean_sort(wavnum+2),'data');
        
        
        data_range = [index(wavnum)-150:1:index(wavnum)+60];
        N = numel(data_range);
        freqE = get(w_clean_sort(wavnum), 'freq');
        dnumE(1) = datenum(get(w_clean_sort(wavnum),'start'));
            for l = 2:numel(dataE)
                dnumE(l) = datenum((l/freqE)/SECSPERDAY+dnumE(1));
            end
            
        sta=get(w_clean_sort(wavnum),'station');
        chan=get(w_clean_sort(wavnum:wavnum+2),'channel');
        
        E_data = dataE(data_range);
        N_data = dataN(data_range);
        Z_data = dataZ(data_range);
        
        ZZ = (1/N)*sum(Z_data.*Z_data);
        ZN = (1/N)*sum(Z_data.*N_data);
        ZE = (1/N)*sum(Z_data.*E_data);
        
        NN = (1/N)*sum(N_data.*N_data);
        NE = (1/N)*sum(N_data.*E_data);
        
        EE = (1/N)*sum(E_data.*E_data);
        
        index_number = find(num_vals == wavnum);
        correlation_matrix{index_number} = [ ZZ ZN ZE ;...
                                             ZN NN NE;...
                                             ZE NE EE];
        
        %eig_val = eig(correlation_matrix)
        %[m, I] = nanmax(data);
        %data_subsetE = dataE()
        ZZ2(index_number) = ZZ;
        ZN2(index_number) = ZN;
        ZE2(index_number) = ZE;
        
        NN2(index_number) = NN;
        NE2(index_number) = NE;
        
        EE2(index_number) = EE;
end

%%
num_sta = numel(index_valuez)
ZZ_sum = sum(ZZ2)/num_sta;
ZN_sum = sum(ZN2)/num_sta;
ZE_sum = sum(ZE2)/num_sta;

NN_sum = sum(NN2)/num_sta;
NE_sum = sum(NE2)/num_sta;

EE_sum = sum(EE2)/num_sta;

sum_corr_matrix = [ ZZ_sum ZN_sum ZE_sum;...
                    ZN_sum NN_sum NE_sum;...
                    ZE_sum NE_sum EE_sum];
[eig_vec, eig_mat] = eig(sum_corr_matrix);



