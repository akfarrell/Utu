function [ivec,w,wsyn] = plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap,wsyn)
%PLOTW_RS processes waveform object and plots as a record section
%
% See examples in run_plotw_rs.m
%
% INPUT
%   w       waveform object (WHAT ADDED FIELDS SHOULD WE REQUIRE?)
% INPUT (OPTIONAL) -- these can be omitted or set as [] to get default settings
%   rssort  =0 input sorting of records
%           =1 azimuthal sorting of records
%           =2 distance sorting of records
%           =3 alphabetical sorting of records (by station name or event ID)
%   iabs    =1 to plot by absolute distance or azimuth (this means unequal
%              vertical spacing between the waveforms)
%   tshift  time shift (in seconds) to apply to the waveforms ([] for none)
%   tmark   absolute time markers (Matlab serial date) ([] for none)
%           note: if tshift varies, then you can't use this option
%   Tfilt   =[Tmin Inf] for high-pass filter
%           =[0 Tmax] for low-pass filter
%           =[Tmin Tmax] for band-pass filter
%   T1      minimum period of filter: =[] for low-pass or no filter
%   T2      maximum period of bandpass: =[] for high-pass or no filter
%   pmax    maximum number of time series per record section
%   iintp   =1 to integrate, =-1 to differentiate, =0 for nothing
%   inorm   =0 for no amplitude scaling
%           =1 for scaling by max(abs(d(t))) per trace
%           =2 for [other options not yet implemented]
%   tlims   time limits for x-axis ([] to use default)
%           note that the 'tshift' field will shift each seismogram
%   nfac    controls spacing between seismograms
%   azstart   azimuthal angle for top record (applies with rssort=1 only)
%   iunit   units for distance record section (applies with rssort=2 only)
%           =1 for km sphere, =2 for deg sphere, =3 for km flat
%   imap    plot a simple map showing the station and source locations
%
% OUTPUT (OPTIONAL)
%   ivec  index in which w(i) are ordered in the plot;
%            w(ivec) gives the plotting order from top to bottom (if iabs=0)
%   w     processed waveform (typically data)
%   wsyn  processed waveform (typically synthetics)
%
% KEY: extract data using getwaveform.m OR use wset.m, prior to running this.
%
% DETAILS ABOUT THE TIME AXIS
%   The following variables are pertinent to the time axis:
%   tshift, tmark, and tlims. tmark is intended to mark absolute times with
%   vertical bars. If tshift varies from station to station, then the time
%   axis is relative time, and tmark cannot be used. If there is no tshift
%   specified, then t=0 will be the earliest time of all the seismograms,
%   and tlims will be specified w.r.t. this t=0.
%
% To print figures to files, set bprint=true below.
%
% FUTURE WORK:
%   - if padding zeros it should be done AFTER demean, taper, etc
%   - allow pmax to represent the first X sorted seismograms, while not
%     plotting the rest (which may have too-high SNR, for example)
%   - combine station and network code to allow for a station to have two
%     different network codes (e.g., MCK)
%   - CAUTION: THERE APPEARS TO BE SINGLE PRECISION VARIABLES ASSOCIATED
%              WITH THE NEW MATLAB RELEASE (2011b)
%
% See also
% /home/admin/share/matlab/PACKAGES/GISMO_OP/GISMO/@waveform/plot.m
%
% Carl Tape 11/21/2011
% Yun Wang 11/2011
%==========================================================================

disp('--> entering plotw_rs.m to plot a record section of waveforms');

spdy = 86400;       % seconds per day
synplot = 'r--';    % plotting a second set of seismograms
global maptime      % used for determining how long each step takes
                    % (see run_getwaveform_short.m)
%--------------------
% check input arguments

disp('input arguments:');
whos w rssort iabs tshift tmark T1 T2 pmax iintp inorm tlims nfac azstart iunit imap wsyn

narg0 = 16;
%if nargin < narg0-1, error('missing %i of the %i input arguments',narg0-nargin,narg0); end
if nargin == narg0
    disp('second waveform object detected -- waveforms will be superimposed');
    isyn = 1;
else
    isyn = 0;
    if nargin == narg0-1
        disp('all input arguments present');
    else
        disp(sprintf('missing %i of the %i input arguments',narg0-nargin,narg0));
        isyn = 0;
        if nargin==1
            rssort = []; iabs = []; tshift = []; tmark = [];
            T1 = []; T2 = []; pmax = []; iintp = []; inorm = [];
            tlims = []; nfac = []; azstart=[]; iunit=[]; imap=[];
        end
    end
end
% set default parameters
if isempty(rssort), rssort=2; end
if isempty(iabs), iabs=0; end
if isempty(iintp), iintp=0; end
if isempty(inorm), inorm=0; end
if isempty(nfac), nfac=1; end
if isempty(iunit), iunit=1; end
if isempty(imap), imap=1; end
% exit here if user enters impermissible values
if ~any(rssort==[0 1 2 3]), error('input rssort = %f must be 0, 1, 2, 3',rssort); end
if ~any(iabs==[0 1]), error('input iabs = %f must be 0 or 1',iabs); end
if ~any(inorm==[0 1]), error('input inorm = %f must be 0 or 1',inorm); end
if ~any(iintp==[-1 0 1]), error('input iintp = %f must be -1 or 0 or 1',iintp); end

% if vertical axis is absolute (distance or azimuth), then plot all
% waveforms on the same record section (pmax huge)
if iabs==1
    pmax = 1000;
    if ~any(rssort==[1 2])
        iabs, rssort
        error('if iabs=1, then rssort=1 (az) or rssort=2 (dist)');
    end
end

[nw,ncomp] = size(w);   % w could be nw x ncomp
if ~any(ncomp==[1:3])
    w = w(:);
    [nw,ncomp] = size(w);
end
w = w(:);               % convert w to vector
tshift = tshift(:);
nseis = length(w);
if isempty(pmax), pmax=nseis; end

% time shifts (allows for alignment of seismograms on a time that varies
% from one trace to the next, say, a P arrival)
if isempty(tshift)
    disp('no time shift applied (default)');
    tshift = zeros(nseis,1);
elseif length(tshift)==1
    disp(sprintf('time shift of %.2f s applied to all waveforms',tshift));
    tshift = tshift*ones(nseis,1);
elseif length(tshift)==nw
    disp(sprintf('input tshift has dimension %i x %i',size(tshift)));
    tshift = tshift(:);
    tshift = repmat(tshift,1,ncomp);
    disp(sprintf('output tshift has dimension %i x %i',size(tshift)));
else
    error('tshift (%i) must be same length as w (%i)',length(tshift),nseis);
end

% time markers (absolute times)
if isempty(tmark)
    nmark = 0;
    disp('no time markers');
else
    if length(unique(tshift)) > 1
        nmark = 0;
        disp('variable time shifts, so there can be no absolute time markers');
    else
        nmark = length(tmark);
        disp(sprintf('%i time markers to plot',nmark));
    end
end

% relative time or absolute time
if length(unique(tshift)) > 1
    itrel = 1;
else
    tshift0 = tshift(1);
    itrel = 0;
end

%--------------------

chans = get(w,'channel');

if nseis==1
    stchan = chans;
else
    % unique channels
    uchan = unique(get(w,'channel'));
    nuchan = length(uchan);
    stchan = [];
    for ii=1:nuchan, stchan = [stchan ' ' uchan{ii}]; end
end

% help waveform/get
% http://kiska.giseis.alaska.edu/Input/celso/matlabweb/waveform_suite/waveform.html
% (this saves the overhead with many calls to the get function)
[starttime,sta,rlat,rlon,elat,elon,edep,eid,mag,netwk] = ...
    getm(w,'start','station','STLA','STLO','EVLA','EVLO','EVDP','KEVNM','mag','KNETWK');

% FUTURE WORK
% NEED TO EXIT IF ANY OF THE ABOVE FIELDS ARE EMPTY (isempty does not work)

% we assume that there is either one event or one station in the set of waveforms
if nseis==1
    nsta = 1;
    neve = 1;
    irs = 1;
    eid = cellstr(eid);
    slabs = cellstr([sta '.' netwk]);
    sta = cellstr(sta);
else
    %nsta = length(sta);     % temporary (MCK with two networks)
    nsta = length(unique(sta));
    neve = length(unique(eid));
    %[nsta,~] = size(unique([rlon(:) rlat(:)],'rows'));
    %[neve,~] = size(unique([elon(:) elat(:)],'rows'));
    if and(neve==1,nsta>0)
        irs=1;  % 1 event, multiple stations
        %nsta = nseis;
        %slabs = sta;
        slabs = repmat(cellstr(''),nseis,1);
        for ii = 1:nseis
            if nuchan > 1
                slabs{ii} = [sta{ii} '.' netwk{ii} '.' chans{ii}];
            else
                slabs{ii} = [sta{ii} '.' netwk{ii}];
            end
        end
    elseif and(nsta==1,neve>0)
        irs=0;  % 1 station, multiple events
        slabs = eid;
        neve = nseis;
    else
        disp(sprintf('nsta = %i, neve = %i -- error with waveform data',nsta,neve));
        disp('check headers KEVNM and station:');
        get(w,'KEVNM')
        get(w,'station')
        disp('check headers KEVNM and station');
        error('must have a single event or station common to all input waveforms');
    end
end
disp(sprintf('%i event, %i stations',neve,nsta));

% reference start time
[tstartmin,imin] = min(starttime);
disp(sprintf('minimum start time of all waveforms is %s (%s)',...
    datestr(tstartmin,31),sta{imin}));
tref = tstartmin;
if and(itrel==1,length(tmark)==1)
    tref = tmark;
end
stref = sprintf('reference time is %s',datestr(tref,31));
disp(stref);
disp('--> this will be subtracted from all time vectors')

% compute distances and azimuths to stations (or events)
if irs==1
    lat1 = elat(1)*ones(nseis,1);
    lon1 = elon(1)*ones(nseis,1);
    lat2 = rlat(:);
    lon2 = rlon(:);
else
    lat1 = rlat(1)*ones(nseis,1);
    lon1 = rlon(1)*ones(nseis,1);
    lat2 = elat(:);
    lon2 = elon(:);
end
%whos lat1 lon1 lat2 lon2
[dist,azi] = distance(lat1,lon1,lat2,lon2, 'degrees');

% range of distances and azimuths
if ~any(iunit==[1 2 3])
    disp('WARNING: iunit not 1,2,3 -- setting iunit = 1');
    iunit = 1;
end
switch iunit
    case 1
        sunit = 'km';
        dist = deg2km(dist);
    case 2
        sunit = 'deg';
    case 3
        % input lon2 and lon1 are assumed to be utmx and utmy,
        % so dist (from above) is over-written
        sunit = 'km';
        dist = 1e-3 * sqrt( (lon2-lon1).^2 + (lat2-lat1).^2 );
end
dran = max(dist) - min(dist);
aran = max(azi) - min(azi);
disp(sprintf('distance range is %.1f - %.1f = %.1f %s',max(dist),min(dist),dran,sunit));
disp(sprintf(' azimuth range is %.1f - %.1f = %.1f',max(azi),min(azi),aran));

% sort
sortlab = {'input','azimuth','distance','label'};
switch rssort
    case 0      % default is no sorting
        ivec = [1:nseis]';
    case 1      % azimuth
        if isempty(azstart)
            azstart = 0;
            disp('WARNING: azstart not specified -- setting azstart = 0');
        end
        [~,ivec] = sort(wrapTo360(azi-azstart)); 
    case 2      % distance
        [~,ivec] = sort(dist);
    case 3      % alphabetical
        [~,ivec] = sort(slabs);
end
ivec = ivec(:);

% labels for record section -- UNSORTED
rlabels = repmat(cellstr(''),nseis,1);
if iabs==0
    for jj=1:nseis
        % variable time shift
        if length(unique(tshift)) > 1
            stshift = sprintf('DT %.1f',tshift(jj));
            %stshift = sprintf('DT %.1f',tshift(jj)-min(tshift));
        else
            stshift = '';
        end
        if irs==1   % 1 event, multiple stations
            rlabels{jj} = sprintf('%s (%.0f, %.0f %s) %s',...
                slabs{jj},azi(jj),dist(jj),sunit,stshift);
        else        % 1 station, multiple events
            rlabels{jj} = sprintf('%s (%.0f, %.0f %s) %.0f km) %s',...
                slabs{jj},azi(jj),dist(jj),sunit,edep(jj),stshift);
        end
    end
else
    rlabels = slabs;
end

%----------------------FILTER AND PLOT----------------------
   
RTAPER = 0.05;
T1 = 0.75;
T2 = 1.5;
% filter
if and(isempty(T1),isempty(T2))
    ifilter = 0;
    disp('NO FILTER WILL BE APPLIED');
    stfilt = '--';
else
    %w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
    ifilter = 1;
    npoles = 2;
    
    % fill gaps with mean value
    w = fillgaps(w,'meanAll');
      
    % these operations might depend on whether the input is displacements
    % (which could have static offsets) or velocities
    disp('pre-processing: detrend, demean, taper');
    w = detrend(w);
    w = demean(w); 
    w = taper(w,RTAPER);
    
    Tmax_for_mHz = 100;     % list mHz, not Hz, for T2 >= Tmax_for_mHz
    if and(~isempty(T1),isempty(T2))
        disp(sprintf('%i-pole low-pass T > %.1f s (f < %.2f Hz)',npoles,T1,1/T1));
        f = filterobject('L',1/T1,npoles);
        stfilt = sprintf('T > %.1f s (f < %.2f Hz)',T1,1/T1);
        if T1 >= Tmax_for_mHz, stfilt = sprintf('T > %.1f s (f < %.1f mHz)',T1,1/T1*1e3); end
        
    elseif and(isempty(T1),~isempty(T2))
        disp(sprintf('%i-pole high-pass T < %.1f s (f > %.2f Hz)',npoles,T2,1/T2));
        f = filterobject('H',1/T2,npoles);        
        stfilt = sprintf('T < %.1f s (f > %.2f Hz)',T2,1/T2);
        if T2 >= Tmax_for_mHz, stfilt = sprintf('T < %.1f s (f > %.1f mHz)',T2,1/T2*1e3); end
        
    elseif and(~isempty(T1),~isempty(T2))
        disp(sprintf('%i-pole band-pass filter between T = %.1f - %.1f s',npoles,T1,T2));
        f = filterobject('B',[1/T2 1/T1],npoles);
        stfilt = sprintf('T = %.1f-%.1f s (%.2f-%.2f Hz)',T1,T2,1/T2,1/T1);
        if T2 >= Tmax_for_mHz, stfilt = sprintf('T = %.1f-%.1f s (%.1f-%.1f mHz)',T1,T2,1/T2*1e3,1/T1*1e3); end
    end
    
    %w = filtfilt(f,w);   % apply filter
    
    try
       w = filtfilt(f,w);   % apply filter
    catch
       for ii=1:nw
           disp(sprintf('%s, %.6f',get(w(ii),'station'),get(w(ii),'freq')));
           w(ii) = filtfilt(f,w(ii));
       end
    end
    
    % apply identical filtering to synthetics, if present
    if isyn==1
        wsyn = detrend(wsyn);
        wsyn = demean(wsyn); 
        wsyn = taper(wsyn,RTAPER);
        wsyn = filtfilt(f,wsyn);
    end
end

% temporary (if no filter specified)
%w = demean(w); 

% integrate or differentiate
% note: units of w will automatically change
if iintp==1
    if ifilter==0, w = detrend(w); end
    w = integrate(w);
    if isyn==1
        if ifilter==0, wsyn = detrend(wsyn); end
        wsyn = integrate(wsyn);
    end
end
if iintp==-1
    w = diff(w);        % help waveform/diff
    if isyn==1, wsyn = diff(wsyn); end
end

% get units, which may have changed
units = get(w,'units');
if nseis==1
    stunit = units;
else
    % unique units
    uunit = unique(get(w,'units'));
    nuunit = length(uunit);
    if nuunit > 1, disp('WARNING: multiple units are present on waveforms'); end
    stunit = [];
    for ii=1:nuunit, stunit = [stunit ' ' uunit{ii}]; end
end

%---------Plot record sections---------

nfig = ceil(nseis/pmax);
fsize = 10;
kk = 0;
az1 = azstart;
azinc = 45;
azbin = wrapTo360(azstart : azinc : azstart+360);
azbin_fwid = 0.1;

% compute amplitude scaling for plots
% NOTE: This will normalize by the full time series, not simply the time
% specified within the window denoted by tlims.
% MORE OPTIONS HERE (correct for geometrical spreading, etc)
nlabs = {'none','max(abs(d_i))'};
nlab = ['norm --> ' nlabs{inorm+1}];
nvec = ones(nseis,1);       % normalization vector for seismograms
wmax = max(max(abs(fillgaps(w,0))));
wsep = 1.5*median(max(abs(fillgaps(w,0)))); % will handle empty records
yshift = wsep*nfac;
if inorm==1
   nvec = max(abs(w)); 
   yshift = nfac;
end
if iabs==1
    if inorm==0
        nvec = nfac*max(max(abs(w)))*ones(nseis,1);
    else
        if rssort==1
            yran=aran;
        else
            yran=dran;
            if iunit==1, yran = km2deg(dran); end
        end
        nvec = nfac * nseis/yran*max(abs(w));
    end
    % sign flip is needed, since y-axis direction is flipped
    nvec = -nvec;
end
disp(sprintf('wmax = %.3e, wsep = %.3e, yshift = %.3e',wmax,wsep,yshift));

% debugging for variable nvec:
%disp(sprintf('summary of nvec (%i): min/mean/median/max = %.3e / %.3e / %.3e / %.3e',...
%    length(nvec),min(abs(nvec)),mean(abs(nvec)),median(abs(nvec)),max(abs(nvec))));

% get time limits for record section
% NOTE: THIS SHOULD BE AVOIDED SINCE IT INVOLVES READING ALL THE WAVEFORMS
if isempty(tlims)
    tic
    tlim1 = zeros(nw,1);
    tlim2 = zeros(nw,1);
    for ii=1:nw
        ti3 = get(w(ii),'timevector');
        % KEY: time vector for plotting
        tplot = (ti3 - tref)*spdy - tshift(ii);
        tlim1(ii) = min(tplot);
        tlim2(ii) = max(tplot);
    end
    disp(sprintf('it took %.1f s to establish the time limits for plotting',toc));
    tlims = [min(tlim1) max(tlim2)];
end

for pp = 1:nfig
    disp(sprintf('record section page %i/%i (max %i per page)',pp,nfig,nseis));
    
    %scrsz = get(0,'ScreenSize');
    %figure('Position',[1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2]);
    %set(gcf,'Color',[1,1,1]); % Set backgroud color to white
    %figure('Position',[1 1 800 1000],'color',[1,1,1]);
    figure; hold on;
    
    if nfig==1, jmax = nseis; else
        if pp==nfig, jmax = mod(nseis,pmax); if jmax==0, jmax=pmax; end
        else jmax = pmax; end
    end
    %disp(sprintf('TESTING: pp=%i, nfig=%i, nseis=%i, pmax=%i, jmax=%i',pp,nfig,nseis,pmax,jmax));
    
    % initialize arrays
    dy = zeros(jmax,1);
    rlabs = repmat(cellstr(''),jmax,1);
    dplotmax = 0;
    for jj=1:jmax       % loop over seismograms
        kk = kk+1;
        ii = ivec(kk);  % key sorting
        rlabs(jj) = rlabels(ii);

        % get seismogram
        [ti3,di3] = getm(w(ii),'timevector','data');
        % KEY: time vector for plotting
        tplot = (ti3 - tref)*spdy - tshift(ii);
        
        if iabs==0
            % plot from top to bottom
            dy(jj) = (jmax + 1 - jj)*yshift;
        else
            % use absolute scaling (e.g., plot records at their actual distance)
            if rssort==1, dy(jj) = azi(ii); else dy(jj) = dist(ii); end
        end
        dplot = di3/nvec(ii);
        dplotshift = dplot + dy(jj);
        % get the max value of the seismogram within the plotted time interval
        btplot = and(tplot > tlims(1),tplot < tlims(2));
        dplotmax = max([dplotmax max(abs(dplot(btplot)))]);
        
        plot(tplot,dplotshift,'b');
             
        % label for amplitude of one record
        if jj==1
            [~,imx] = max(abs(di3));
            stmx = sprintf('%s max %.2e %s at t = %.1f s',...
                get(w(ii),'station'),di3(imx),get(w(ii),'units'),tplot(imx));
        end
        
        if isyn==1
            % be careful about the reference time
            [tisyn,disyn,tstartsyn] = getm(wsyn(ii),'timevector','data','start');
            tplot = (tisyn - tstartsyn)*spdy - tshift(ii);
            dplotshift = disyn/nvec(ii) + dy(jj);
            plot(tplot,dplotshift,synplot);
        end
        
        % partition if plotting by azimuth (see azbin above)
        if and(rssort==1,iabs==0)
            % azimuth of previous (az0) and current (az1) stations in the sorted list
            az0 = az1;
            az1 = azi(ii);
            % (this boolean clause could probably be simplified)
            if or( and(az1 > az0, any(and(az1 > azbin, azbin > az0))), and(az1 < az0,  any(or(az1 > azbin, azbin > az0))) ) 
                % note: it would nice if these bars extended to the LEFT,
                % outside the plotting axes
                tp1 = tlims(1);
                tp2 = tlims(1) + azbin_fwid*(tlims(2)-tlims(1));
                %disp(sprintf('%i %i %.2f %.2f',pp,jj,tp1,tp2));  % testing
                if ~exist('wsyn','var'), pc = 'r'; else pc = 'k'; end
                plot([tp1 tp2],[1 1]*dy(jj)+yshift/2,pc,'linewidth',2);
            end
        end
    end  % jj (loop over seismograms)
    
    %if isempty(tlims), tlims = [min(tplot) max(tplot)]; end
    %if isempty(tlims), tlims = [min(tlim1) max(tlim2)]; end
    if iabs==0
        % this may chop off a large-amplitude trace near the boundary, but
        % at least the gap will be the same for a sequence of plots
        ylims = yshift*[0 jmax+1];
        
        % this will provide more space if one of the traces has a large
        % amplitude, but the gap may vary for a sequence of plots
        %ylims = yshift*[1 jmax] + dplotmax*[-1 1];
        
    else
        % using actual values of distance or azimuth
        if rssort==1
            ylims = [min(azi)-5 max(azi)+5];
            if max(azi)-min(azi) > 300
               ylims = [-10 370]; 
            end
        else
            ylims = [min(dist)-0.1*dran max(dist)+0.1*dran];
        end
    end
    disp(sprintf('tlims = %.2f to %.2f',tlims(1),tlims(2)));
    disp(sprintf('ylims = %.3e to %.3e (yshift = %.3e, jmax=%i)',ylims(1),ylims(2),yshift,jmax));
    ax0 = [tlims ylims];
    % plot tshift marker
    plot([0 0],ax0(3:4),'r','linewidth',2);
    % plot absolute-time marker
    for mm=1:nmark
        tm = (tmark(mm) - tstartmin)*spdy - tshift0;
        plot(tm*[1 1],ax0(3:4),'r','linewidth',2);
    end
    axis(ax0);
    % PLOT LABELS FOR EACH WAVEFORM
    %set(gca,'ytick',dy,'yticklabel',rlabels(ivec),'fontsize',fsize-2);
    if iabs==0
        set(gca,'ytick',flipud(dy(:)),'yticklabel',flipud(rlabs(:)),'fontsize',fsize-2);
    else
        t0 = tlims(2) + 0.02*diff(tlims);
        text(t0*ones(jmax,1),dy,rlabs,'fontsize',fsize-2,'interpreter','none');
        if rssort==1, ylabel('Azimuth, degrees'); else ylabel(['Distance, ' sunit]); end
        set(gca,'ydir','reverse');  % note: also need to flip the sign of the waveforms
    end
    xlabel('Time (s)','FontWeight','Bold');

    % title string explaining the time axis
    if itrel==0
        t1a = tstartmin + (tshift0+tlims(1))/spdy;
        t2a = t1a + diff(tlims);
        stdur = sprintf('%s + %.2f s',datestr(double(t1a),31),t2a-t1a);
    else
        % relative time shifts
        stdur = sprintf('variable time shifts: %s',stref);
%         if length(tmark)==1
%             stdur = sprintf('variable time shifts: %s',stref);
%         else
%             % just pick the iith record as an example -- this should be the
%             % bottom time series on the record section
%             ix = ii;
%             t1a = tstartmin + (tshift(ix)+tlims(1))/spdy;
%             t2a = t1a + diff(tlims);
%             stdur = sprintf('variable time shifts: w(%i) (%s) is %s + %.2f s',...
%                 ix,sta{ix},datestr(double(t1a),31),t2a-t1a);
%         end
    end
    
    % plot title
    % note: might want to label the channel, too
    %if ifilter==1, stfilt = sprintf('T = %.1f-%.1f s (%.2f-%.2f Hz)',T1,T2,1/T2,1/T1); else stfilt = '--'; end
    st0 = sprintf('%s [%s, %s]',stchan,stunit,stfilt);
    if irs==1
       stline1 = sprintf('event %s (%s, M%.1f, %.1f, %.1f, z = %.1f km)',...
            eid{1},datestr(starttime(1),29),mag(1),elon(1),elat(1),edep(1));
       stline2 = sprintf('%i / %i seismograms (%i stations) ordered by %s, %s',...
           jmax,nseis,nsta,sortlab{rssort+1},nlab);
    else
       stline1 = sprintf('station %s (%.1f, %.1f)',sta{1},rlon(1),rlat(1));
       stline2 = sprintf('%i / %i seismograms (%i events) ordered by %s, %s',...
           jmax,nseis,neve,sortlab{rssort+1},nlab);   
    end
    title({[stdur ';  ' stmx],sprintf('%s %s',st0,stline1),stline2},'fontsize',fsize,'interpreter','none');

    % print figure
    % To make composite PDF:
    % for file in `ls rs_*.ps` ; do ps2pdf $file ; done
    % pdcat -r rs_*.pdf all_rs.pdf
    bprint = false;
    if bprint, orient tall;
        ofile = sprintf('rs_%s',eid{1});
        if nfig > 1, ofile = sprintf('rs_%s_p%2.2i',eid{1},pp); end
        print(gcf,'-dpsc',ofile);
    end

end

%--------------------------------------------------------------------------
% plot map
% future work would be to add some more options for alaska plots (alaska_basemap.m)
startmaptime = tic;
if imap==1
    if iunit==3, fac = 1e-3; else fac = 1; end
    iplotsrc = 1;
    fsize = 10;
    figure; hold on;
    plot(rlon*fac,rlat*fac,'bV');
    if iplotsrc==1, plot(elon*fac,elat*fac,'kp','markersize',20,'markerfacecolor','r'); end
    % plot station labels (or eid labels)
    if irs==1
        text(rlon*fac,rlat*fac,sta,'fontsize',fsize,'interpreter','none');
    else
        text(elon*fac,elat*fac,eid,'fontsize',fsize,'interpreter','none');
    end
    title([stline1 ' ' stchan],'fontsize',10,'interpreter','none');
    %print(gcf,'-dpsc',sprintf('plotw_rs_map_%s',eid{1}));
end
maptime = toc(startmaptime);

disp('--> leaving plotw_rs.m');

%==========================================================================
