function [P_az, P_inc] = directionality(eq_info, earthquake_number, index_valuez, fil, order)
close all;

%eq_info = eq(earthquake_number)
%function call: [P_az, P_inc] = directionality(eq(earthquake_number), earthquake_number, index_values, fil, order)

% startup_seismiclab
% addpath(genpath('/raid/apps/src/GEOTOOLS/matlab_util'))
% addpath('data_func_matTaup/')
% addpath('latlonutm/Codes/matlab/')
% addpath('readhgt/')
% NEWGISMODIR=fullfile('/raid/apps/matlab/toolbox/GISMO/startup_GISMO.m');
% rmpath(genpath(NEWGISMODIR))
% addpath('/raid/apps/src/gismotools/GISMO')
% startup_GISMO

ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
scnl = scnlobject('*', '*', 'PL')

utu_lat = -22.27;
utu_lon = -67.18;

w_raw = waveform(ds, scnl, eq_info.snum, eq_info.enum);

%chan1 = get(w_raw, 'channel')
%chan1(1:3:numel(chan1)) %run to check that this is all HHE
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

SECSPERDAY = 60 * 60 * 24;
%index_values = load('index_values.mat');
w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));

if earthquake_number == 4
    w_clean(12*3:12*3+2)=[]; %remove station PLLL - instrument noise KTSZ3
elseif earthquake_number == 5
    w_clean(1:7*3)=[]; %remove all LZ stations - KTSZ4
    w_clean(18*3:18*3+2)=[]; %remove PLSP - waveform unlike others
elseif earthquake_number == 7
    w_clean(6*3:6*3+2)=[]; %remove station PLQU - noise JSZ2
elseif earthquake_number == 8    
    w_clean(17*3:17*3+2)=[]; %remove station PLQU - noise JSZ3
    w_clean(20*3:20*3+2)=[]; %remove station PLSP - noise JSZ3
    w_clean(9*3:9*3+2)=[]; %remove station PLDK - noise JSZ3
elseif earthquake_number == 10
    w_clean(9*3:9*3+2)=[]; %remove station PLSM - noise SSSZ1
elseif earthquake_number == 11
    w_clean(6*3:6*3+2)=[]; %remove station PLQU - noise SSSZ2
elseif earthquake_number == 12
    w_clean(23*3:23*3+2)=[]; %remove station PLTP - noise SSSZ3
end

n = 3; 
index = reshape(repmat(index_valuez(:).',n,1),1,[]);

for counter = 1:numel(order)
    order(counter) = order(counter)*2+order(counter)-2;
end
order
numel(index)
orden = zeros(1,numel(index));

for counter = 1:numel(order)
    orden(counter*2+counter-2) = order(counter);
    orden(counter*2+counter-1) = orden(counter*2+counter-2)+1;
    orden(counter*2+counter) = orden(counter*2+counter-2)+2;
end

orden
%order = reshape(repmat(order(:).',n,1),1,[])

chan2 = get(w_clean, 'channel');
chan2(1:3:numel(chan2)); %run to check that this is all HHE
w_clean_sort = w_clean(orden);

%%


chan = get(w_clean_sort, 'channel')
chan(1:3:numel(chan)) %run to check that this is all HHE
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

%%
l1 = eig_mat(3,3);
l2 = eig_mat(2,2);
l3 = eig_mat(1,1);

u1 = eig_vec(:,3);
u2 = eig_vec(:,2);
u3 = eig_vec(:,1);

P_az = atand((u1(2)*sign(u1(1)))/(u1(3)*sign(u1(1))))
P_inc = acosd(abs(u1(1)))
end
