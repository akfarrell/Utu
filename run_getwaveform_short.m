function cross_corr(eq, spectrogramFraction, mycolormap)
%
% Abbreviated version of run_getwaveform.m and getwaveform_input.m
%


spdy = 86400;

% record section plotting parameters
isort = 2;      % =1 by azimuth, =2 by distance
iabs = 0;
T1 = 0.1;
T2 = 2;
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

%========================================
% extracted from getwaveform_input.m

idatabase = 6;

% source parameters (some can be empty)
originTime = datenum('2011-08-21 12:46:50');
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
startTime = eq.snum;
endTime = eq.enum;
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
[w,s,site,sitechan] = getwaveform(idatabase,startTime,endTime,chan,iint,iprocess,cutoff,samplerate,stasub,sacdir,originTime,elat,elon,edep_km,mag,eid);

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
    plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);
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
plotw_rs(wfilt,0,iabs,[],tmark,[],[],pmax,iintp,inorm,tlims,2,azcen,iunit,imap);

% see nenana_trigger.m for example
get(wfilt,'station')
C = correlation(wfilt);
C = xcorr(C);
lagtimes = get(C,'lag');
tshift = -lagtimes(1,:);
%plot(C);
plot(C,'lag'); set(gca,'xticklabel',get(wfilt,'station'),'yticklabel',get(wfilt,'station'))
plot(C,'corr'); set(gca,'xticklabel',get(wfilt,'station'),'yticklabel',get(wfilt,'station'))

plotw_rs(wfilt,isort,iabs,tshift,tmark,[],[],pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);

figure; hold on;
[rlat,rlon,sta] = getm(wfilt,'STLA','STLO','station');
scatter(rlon,rlat,20^2,tshift,'filled');
text(rlon,rlat,sta); colorbar

