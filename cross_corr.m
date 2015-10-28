function tshift = cross_corr(eq, fil)
%
% Abbreviated version of run_getwaveform.m and getwaveform_input.m
%

name = eq.name;
spdy = 86400;

% record section plotting parameters
isort = 2;      % =1 by azimuth, =2 by distance
iabs = 0;
T1 = fil(1);
T2 = fil(2);
trshift = 0;
tmark = [];
pmax = 50;
iintp = 0;
inorm = 1;
nfac = 1;
azcen = [];
iunit = 0;
tlims = [];     % time limits for plotting
imap = 1;
secs2days = 60*60*24;

%========================================
% extracted from getwaveform_input.m

idatabase = 6;

% source parameters (some can be empty)
originTime = eq.evtime;
elat = eq.lat;
elon = eq.lon;
edep_km = eq.depth;
eid = [];
mag = eq.mag;

chan = {'HHZ'};

duration_s = 10;
oshift = 1;
stasub = [];
cutoff = [];
samplerate = [];

% tshift (and trshift) are for plotting record sections only
tshift = oshift + trshift;

%startTime = originTime - max(oshift)/spdy;
%startTime = originTime - oshift/spdy;
%endTime   = originTime + duration_s/spdy;
startTime = eq.snum; %-(5/secs2days);
endTime = eq.enum; %+(5/secs2days);
dur_dy = endTime-startTime;
fprintf('origin time is %s\n',datestr(originTime,'yyyy-mm-dd HH:MM:SS.FFF'))
fprintf('startTime is %s\n',datestr(startTime,31));
fprintf('total length of time requested: %.2f s (= %.2f min = %.2f hours)\n',...
    dur_dy*spdy,dur_dy*3600,dur_dy*24);

%========================================

% additional user parameters
%sacdir = './';      % =[] to return waveform object only
sacdir = [];
iint = 0;            % integrate waveforms: =1 for displacement, =0 for velocity
iprocess = 1;        
%irs = input('Type 1 to plot record section or 0 to not, then enter: ');
irs = 1;

% get waveform object (optional: write sac files to a directory)
% I = [idatabase,startTime,endTime,chan,iint,iprocess,cutoff,samplerate,stasub,sacdir,originTime,elat,elon,edep_km,mag,eid];
% for i = 1:numel(I)
%     iscell(I(i))
% end

[w,s,site,sitechan] = getwaveform_test(idatabase,startTime,endTime,chan,iint,iprocess,cutoff,samplerate,stasub,sacdir,originTime,elat,elon,edep_km,mag,eid);
if strcmp(eq.name, 'ESZ1')
    starttimes  = {'5/19/2011 15:02:02', '5/19/2011 15:02:02', '5/19/2011 15:02:03', '5/19/2011 15:02:03', '5/19/2011 15:02:03', '5/19/2011 15:02:04', '5/19/2011 15:02:04', '5/19/2011 15:02:04', '5/19/2011 15:02:04', '5/19/2011 15:02:06', '5/19/2011 15:02:06', '5/19/2011 15:02:06', '5/19/2011 15:02:07', ...
        '5/19/2011 15:02:07', '5/19/2011 15:02:08', '5/19/2011 15:02:08', '5/19/2011 15:02:08', '5/19/2011 15:02:09', '5/19/2011 15:02:09', '5/19/2011 15:02:09', '5/19/2011 15:02:10', '5/19/2011 15:02:11', '5/19/2011 15:02:11', '5/19/2011 15:02:12', '5/19/2011 15:02:13', '5/19/2011 15:02:13'};
    endtimes = {'5/19/2011 15:02:20', '5/19/2011 15:02:20', '5/19/2011 15:02:21', '5/19/2011 15:02:21', '5/19/2011 15:02:21', '5/19/2011 15:02:22', '5/19/2011 15:02:22', '5/19/2011 15:02:22', '5/19/2011 15:02:22', '5/19/2011 15:02:24', '5/19/2011 15:02:24', '5/19/2011 15:02:24', '5/19/2011 15:02:25', ...
        '5/19/2011 15:02:25', '5/19/2011 15:02:26', '5/19/2011 15:02:26', '5/19/2011 15:02:26', '5/19/2011 15:02:27', '5/19/2011 15:02:27', '5/19/2011 15:02:27', '5/19/2011 15:02:28', '5/19/2011 15:02:29', '5/19/2011 15:02:29', '5/19/2011 15:02:30', '5/19/2011 15:02:31', '5/19/2011 15:02:31'};
    w = subtime(w,starttimes,endtimes);
end

if strcmp(eq.name, 'ESZ1') %change w_clean_sort for ESZ1 from 26x26 to 1x26
    for wavenum = 1:26
        w2(wavenum) = w(wavenum, wavenum);
    end
    w = w2;
end

if strcmp(eq.name, 'KTSZ3')
    w(12)=[]; %remove station PLLL - instrument noise KTSZ3
elseif strcmp(eq.name, 'KTSZ4')
    w(1:7)=[]; %remove all LZ stations - KTSZ4
    w(18)=[]; %remove PLSP - waveforms unlike other stations
elseif strcmp(eq.name, 'JSZ2')
    w(6)=[]; %remove station PLQU - noise JSZ2
elseif strcmp(eq.name, 'JSZ3')   
    w(17)=[]; %remove station PLQU - noise JSZ3
    w(20)=[]; %remove station PLSP - noise JSZ3
    w(9)=[]; %remove station PLDK - noise JSZ3
elseif strcmp(eq.name, 'SSSZ1')
    w(9)=[]; %remove station PLSM - noise SSSZ1
elseif strcmp(eq.name, 'SSSZ2')
    w(6)=[]; %remove station PLQU - noise SSSZ2
elseif strcmp(eq.name, 'SSSZ3')
    w(23)=[]; %remove station PLTP - noise SSSZ3
end

%------------------------------------------------------
% FROM HERE ON OUT, PLOTTING AND CHECKING ONLY

whos w s site sitechan

% check indexing and units
for ii=1:length(w)
   disp(sprintf('%3i %7s %3s %6s %6s %6s %10s',ii,get(w(ii),'channel'),get(w(ii),'KNETWK'),...
       get(w(ii),'station'),char(site(ii).sta),char(sitechan(ii).sta),get(w(ii),'units')));
end

% plot record section
if and(irs==1,~isempty(w))
    % assume all waveforms have the same event ID
    keid = get(w(1),'KEVNM');

    % plot record section
    T1
    T2
    plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,2,azcen,iunit,imap);
    % note: set imap = 0 in plotw_rs.m to NOT plot the map (and print the
    %       record section with the following command)
    iprint = 0;
    if iprint==1, orient tall; print(gcf,'-dpsc',sprintf('%s_%i_%i_isort%i_rs',get(w(1),'KEVNM'),round(T1),round(T2),isort)); end
    %orient tall; print(gcf,'-dpsc','sumatra2012hf_fullA');
end

imap = 0;

% now plot with time shifts
%tshift = 5*rand(length(w),1);
[svec,wfilt] = plotw_rs(w,isort,iabs,[],tmark,0.2,4,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);
wfilt = wfilt(svec);
%figure 4
plotw_rs(wfilt,0,iabs,[],tmark,[],[],pmax,iintp,inorm,tlims,2,azcen,iunit,imap);
% see nenana_trigger.m for example
get(wfilt,'station')
C = correlation(wfilt);
C = xcorr(C);
lagtimes = get(C,'lag');
tshift = -lagtimes(1,:);
if strcmp(name, 'KTSZ1')
    tshift(3) = tshift(3)+1.3000;
    tshift(4) = tshift(4)+1.1500;
    tshift(11) = tshift(11)+1.1000;
    tshift(13) = tshift(13)+1.1000;
    tshift(14) = tshift(14)+1.1500;
elseif strcmp(name, 'KTSZ4')
    tshift(1) = tshift(1)-3.500;
    tshift(2) = tshift(2)-3.500;
    tshift(3) = tshift(3)-5.500;
    %tshift(7) = tshift(7)-3.500;
    %tshift(7) = tshift(7)+10.000;
    tshift(11) = tshift(11)-3.500;%
    %tshift(12) = tshift(12)-3.500;
    %tshift(12) = tshift(12)+5.500;
    tshift(16) = tshift(16)-3.500;%
    tshift(17) = tshift(17)-3.500;
    %tshift(18) = tshift(18)-3.500;
    %tshift(19) = tshift(19)-3.500;
    %tshift(21) = tshift(21)-3.500;
    tshift(18) = tshift(18)-3.500;
    tshift(20) = tshift(20)-3.500;
    for i=1:numel(tshift)
        tshift(i) = tshift(i)+3.500;
    end
elseif strcmp(name, 'JSZ1')
    tshift(2) = tshift(2)+2.500;
    tshift(4) = tshift(4)+0.100;
    tshift(5) = tshift(5)-2.500;
    tshift(8) = tshift(8)-0.200;
    tshift(12) = tshift(12)+0.100;
elseif strcmp(name, 'JSZ2')
    tshift(3) = tshift(3)-1.7;
    tshift(4) = tshift(4)-1.7;
elseif strcmp(name, 'JSZ3')
    tshift(2)=tshift(2)-2.1;
    tshift(3)=tshift(3)-4.3;
    tshift(4)=tshift(4)+2.1;
    tshift(5)=tshift(5)+5.8;
    tshift(7)=tshift(7)-1;
    tshift(8)=tshift(8)-4.6;
    tshift(10)=tshift(10)+0.05;
    tshift(13)=tshift(13)+5.5;
    tshift(17)=tshift(17)+5.5;
    tshift(19)=tshift(19)+5.5;
    tshift(20)=tshift(20)+6;
    tshift(22)=tshift(22)-1.5;
elseif strcmp(name, 'JSZ4')
    tshift(2)=tshift(2)-1.2;
    tshift(6)=tshift(6)+4;
    tshift(10)=tshift(10)-4.3;
    tshift(12)=tshift(12)-1.2;
elseif strcmp(name, 'SSSZ1')
    tshift(2)=tshift(2)-0.4;
    tshift(3)=tshift(3)+2.1;
    tshift(7)=tshift(7)+4.0;
    tshift(10)=tshift(10)+4.1;
elseif strcmp(name, 'SSSZ2')
    tshift(3)=tshift(3)+4;
    tshift(5)=tshift(5)-3.5;
    tshift(6)=tshift(6)+4;
    tshift(7)=tshift(7)+2;
    tshift(8)=tshift(8)+5;
    tshift(9)=tshift(9)+4;
    tshift(10)=tshift(10)+6;
    tshift(11)=tshift(11)+1.5;
end

%plot(C);
%figure 5
plot(C,'lag'); set(gca,'xticklabel',get(wfilt,'station'),'yticklabel',get(wfilt,'station'))
%figure 6
plot(C,'corr'); set(gca,'xticklabel',get(wfilt,'station'),'yticklabel',get(wfilt,'station'))
%figure 7
plotw_rs(wfilt,isort,iabs,tshift,tmark,[],[],pmax,iintp,inorm,tlims,2,azcen,iunit,imap);

%figure 8
figure; hold on;
[rlat,rlon,sta] = getm(wfilt,'STLA','STLO','station');
scatter(rlon,rlat,20^2,tshift,'filled');
text(rlon,rlat,sta); colorbar
