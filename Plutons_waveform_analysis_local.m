%load the waveform data from Antelope database dbmerged
ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
scnl = scnlobject('*', 'HHZ', 'PL');
snum = datenum(2011, 5, 19, 15, 1, 55);
enum = datenum(2011, 5, 19, 15, 3, 0);
w_raw = waveform(ds, scnl, snum, enum);
%disp(w(1)) %check to make sure that the above worked

%get information for these waveforms
stations = [get(w_raw, 'station')];
channels = [get(w_raw, 'channel')];
calibs = [get(w_raw, 'calib')];

% %% Plot waveforms to see if the processing worked
% figure(1)
% len = numel(w_raw)
% for i=1:3
%     %subplot(len, 1, i)
%     figure
%     plot(w(i))
% end

%%
w_clean = waveform_clean(w_raw);

%%

f = filterobject('b', [0.8 25], 2);
w_clean_filt = filtfilt(f, w_clean);

s = spectralobject(1024, 924, 5, []);
figure(2)
specgram2(s, w_clean_filt(1))

figure(3)
r = specgram_iceweb(s, w_clean_filt, 0.5);