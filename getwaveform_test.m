function [w,s,site,sitechan] = getwaveform(idatabase,startTime,endTime,chan,...
    iint,iprocess,cutoff,samplerate,stasub0,sacdir,originTime,elat0,elon0,edep,emag,keid)
%GETWAVEFORM extract waveforms from several antelope databases
% 
% Example function to get waveforms and write sac files with appropriate
% headers. Uses GISMO functions and sac utilities {rsac.m,wsac.m,lh.m,ch.m}.
%
% 1. Extract waveforms from the specified databases
% 3. Return [w,s,site,sitechan] for plotting purpose and other possible usage 
% 3. optional: save waveforms as sac files with headers
%
% OUTPUT
%   w           waveform object
%   s           station list array
%   site        station struct array with fields: sta,ondate,offdate,lat,lon,elev
%   sitechan    station channel struct array with fields: sta,chan,ondate,offdate,hang,vang
%
% INPUT
%   idatabase   database indices (see db_index.m)         
%   startTime   waveform window start time in Matlab time format
%   endTime     waveform window end time in Matlab time format
%   chan        instrument channel cell array, e.g. chan={'BHZ','BHN','BHE'};
%   iint        =1 to integrate, =-1 to differentiate, =0 for nothing
%   iprocess    iprocess = 0 to retrieve raw data in DIGITAL COUNT
%               iprocess = 1 to retrieve data with calibration factor applied,
%                   but without normalized instrument response
%               iprocess = 2 to retrieve data with both calibration factor and
%                   normalized instrument response applied
%                   When iprocess = 2, data processing precedure includes
%                   removing mean using demean.m, removing linear trend
%                   using "detrend.m" and deconvolving instrument 
%               response using "response_apply.m". When iprocess = 0, one can leave "cutoff"
%               and "samplerate" blank []
%   cutoff      cutoff frequency array [low high] for bandpass filtering
%               during instrument deconvolution (iprocess=2).
%               Make sure low >= 4/(length of wavefrom) and high <=min(Nyquist)/4
%               =[] to use default value [4/(length of wavefrom) (min(Nyquist)/4)]
%   samplerate  rate for resampling waveform object
%               If [], keep original sampling rate as in the datebase.
%   stasub      restrict stations into a subset (length 0,2,4,6)
%                 =[] to ignore
%                 =[dmin_km dmax_km] for stations within a range of distances from (elon,elat)
%                 =[lonmin lonmax latmin latmax] for lon-lat box
%                 =[clon clat dmin_deg dmax_deg azmin azmax] for stations
%                  within the specified distance range of (clon,clat) and
%                  within the specified azimuth range
%   sacdir      =[] to NOT write any sac files
%               directory where subfolders with sac files will be written
%                 iprocess = 0, subfolder = 'RawData_COUNT';
%                 iprocess = 1, subfolder = 'ProcessedData1/';
%                 iprocess = 2, subfolder = 'ProcessedData2/';
%   elat        event latitude (or reference point latitude)
%   elon        event lontitude (or reference point longitude)
% OPTIONAL HEADERS
%   originTime  event origin time in Matlab time format (for sac header 'o')
%   edep        event depth (km)
%   emag        emag
%   keid        event ID (string); =[] to use default date-time
%
% Output sac file name looks like: 20000504_042133.FIB.BHZ.AK.sac
% 20000504_042133: wavefrom start time
% FIB: station name
% BHZ: channel name
% AK:  network name
% Network name explanation: (assigned by the FDSN)
% Network name begin with X Y Z shows a temporary experiment 
% Description on permanent networks -- http://www.iris.edu/stations/networks.txt
% Description on temporary networks -- http://www.iris.edu/SeismiQuery/tempNets.phtml
%                                      or http://www.iris.edu/mda#tnetlist
% AK: Alaska Regional Network
% AT: Alaska Tsunami Warning Seismic System
% AU: Geoscience Australia
% AV: Alaska Volcano Observatory (AVO)
% BK: Berkeley Digital Seismograph Network
% CI: Caltech Regional Seismic Network
% CN: Canadian National Seismic Network
% CU: Caribbean Network (USGS)
% GT: Global Telemetered Seismograph Network (USAF/USGS)
% IC: New China Digital Seismograph Network
% II: IRIS/IDA
% IM: International Miscellaneous Stations
% IU: IRIS/USGS
% MS: IRIS/Singapore,Singapore National Network
% TA: USArray Transportable Array (Earthscope/IRIS)chan
% TS: TERRAscope (Southern California Seismic Network)
% US: US National Seismic Network
% XE: BEAAR -- Broadband Experiment Across Alaskan Range
% XM: Alaska Peninsul BB
% XR: ARCTIC -- Structure and Rotation of the Inner Core
% YM: Denali fault aftershocks RAMP
% YV: MOOS -- Multidisiplinary Observations of Subduction
% XF: YAHTSE 
% XP: PLUTONS -- Investigating the pluton growth and volcanism - Central Andes 
%
% Yun Wang & Carl Tape 11/2011
%
% FUTURE WORK:
% 1. allow for a simple scnl input -- for example someone wants to
%    deconvolve the instrument response for a single station 
% 2. make elon and elat optional
%
% -------------------------------------NOTES-------------------------------
% 1. CALIBRATION_APPLIED: 
%    In "uaf_continuous" database, some records have calibration missing or
%    wrong number which doesn't match the calibration table. This reflects
%    a bug in Antelope/Matlab routines that interacts with the database. 
%    The solution is applying calibration correction based on calib in
%    wfdisc and calibration table.
% 2. UNITS:
%    For some waveforms "units" are "null", which means unit is missing
%    from "wfdisc" (not necessarily missing calibration factor).
%    In this case, call function "segtype2units" to pull unit out of 
%    "calibration" table (segtype field)
% -------------------------------------------------------------------------

disp('==================== entering getwaveform.m ====================');
    % start the clock to time how long this takes

%----------------------------------------
% ADDITIONAL USER PARAMETERS

% search 'global' below
global FCUT1_PAR FCUT2_PAR ell0 spdy npoles elon elat stasub ideconspeed,


% we could instead use nested functions
stasub = stasub0;
elon = elon0;
elat = elat0;

% parameters associated with deconvolution
npoles = 2; 
FCUT1_PAR = 4.00;   % cutoff(1) >= FCUT1_PAR / duration
FCUT2_PAR = 0.50;   % cutoff(2) <= FCUT2_PAR * fNyq (original was 0.25)
%FCUT2_PAR = 0.85;   % temporary: to allow higher frequencies in deconvolution
NFRAC = 0.005;      % fraction of record to cut from both ends due to
                    % spurious signals at the endpoints
                    
% for computing azimuths and distances in Matlab we choose an ellipsoid (see discussion below)
% NOTE: units of ellipsoid determine the units of arc-distance in 'distance'
etype = 'ellipsoid';
%etype = 'wgs84';        % almanac('earth')
%etype = 'sphere'; 
ell0 = almanac('earth',etype,'kilometers');

spdy = 86400;   % seconds per day

NDB = 8;    % number of databases available at UAF (db_index.m)

%----------------------------------------

if ~isempty(samplerate)
    if ~(round(samplerate) == samplerate) 
       error('samplerate needs to be an integer');
    end
end

% event headers
% note 1: if [] is input, then sac headers edep and emag will be -nan
% note 2: we should consider making elon and elat OPTIONAL also
if isempty(elat), error('elat must not be empty'); end
if isempty(elon), error('elon must not be empty'); end
if isempty(edep), edep = NaN; end
if isempty(emag), emag = NaN; end
 % etag: eid 'tag' that is used for the KEVNM header added to w object
if isempty(keid)
    if isempty(originTime)
        etag = datestr(startTime,'yyyymmddHHMMSSFFF');
    else
        etag = datestr(originTime,'yyyymmddHHMMSSFFF');
    end
else
    if isnumeric(keid)
        etag = num2str(keid);
    else
        etag = keid;
    end
end   

% save sac files or not
if isempty(sacdir)
    isavesac = 0;
else
    isavesac = 1;
    if ~isempty(originTime)
        botime = true;
        disp('output sac files will have o header set to input origin time');
    else
        botime = false;
        disp('output sac files will have NOT have o header set');
    end
    
    %------------Directory to store sac files-----------%
    if iprocess == 0
        sacdirsub = strcat(sacdir,'RawData_COUNT/');
    elseif iprocess == 1
        sacdirsub = strcat(sacdir,'ProcessedData1/');
    elseif iprocess == 2
        sacdirsub = strcat(sacdir,'ProcessedData2/');
    else
        error('iprocess should be 0, 1 or 2');   
    end
    disp(['OUTPUT BASE DIRECTORY FOR SAC FILES: ' sacdirsub]);
end

% excluding stations
if isempty(stasub)
   disp('no stations will be excluded');
else
   if length(stasub)==2
       disp(sprintf('keep stations that are between %.1f km and %.1f km of (%.2f, %.2f)',...
           stasub(1),stasub(2),elon,elat));
   elseif length(stasub)==4
       disp(sprintf('keep stations inside the lon-lat box [%.2f %.2f %.2f %.2f]',stasub));
   elseif length(stasub)==6
       disp(sprintf('keep stations that are %.2f and %.2f deg of (%.2f, %.2f)',stasub([3 4 1 2])));
       disp(sprintf('   and have azimuths between %.2f and %.2f deg',stasub(5:6)));
   else
       error('stasub must be [] or length 2, 4, or 6');
   end
end

% scnlobject(station, channel, network, location);
% network and location are not actually used in "waveform.m"
% e.g., BHZ becomes BHZ*
nchan = length(chan);
for ii=1:nchan, chanex(ii) = strcat(chan(ii),'*'); end
scnl = scnlobject('*',chanex,'','');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Extract waveforms from database and save unfiltered waveforms as sac files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------READ WAVEFORMS INTO MATLAB-----------------%
% waveforms have calibration applied already

if isavesac==1
    evtdir = char(strcat(sacdirsub,etag,'/'));
    if exist(evtdir,'dir')
        fprintf('we.re here, will try to not rmdir\n');
%        rmdir(evtdir,'s');
    end
    mkdir(evtdir);
    disp(['OUTPUT DIRECTORY FOR SAC FILES: ' evtdir]);
else
    evtdir = './';
end

% keep log file
logfile = [evtdir 'log.txt'];
fid = fopen(logfile,'w');
iflag = 0;

if length(unique(idatabase)) ~= length(idatabase)
    idatabase
    error('repeated index in idatabase');
end
if idatabase==0
    idatabase = [1:NDB];
end
ndatabase = length(idatabase);

w = [];
s = [];
cha = [];
site = [];
sitechan = [];
disp(sprintf('looping over %i/%i databases',ndatabase,NDB));
for ii=1:ndatabase
    idb = idatabase(ii);
    disp(sprintf('========== database %i (%i/%i) ==========',idb,ii,NDB));
    % KEY COMMAND TO GET WAVEFORMS
    [w1,s1,cha1,site1,sitechan1,iflag] = getw(idb,scnl,iflag,fid,...
        startTime,endTime,iprocess,cutoff,samplerate);
    w = [w w1];
    s = [s s1];
    cha = [cha cha1];
    site = [site site1];
    sitechan = [sitechan sitechan1];
end

% close log file
fclose(fid);

% POST-PROCESSING OF WAVEFORMS AND HEADERS

if ~isempty(w)
    % number of waveforms
    nseis = length(w); 
    
    % resample waveform in order to speed up deconvolution
    % WARNING: THIS IS ONLY ACCURATE IF YOU'RE FAR FROM THE NYQUIST;
    %   NEAR THE NYQUIST, THE INSTRUMENT RESPONSE TAKES INTO ACCOUNT
    %   THE ORIGINAL SAMPLING RATE.
    if and(~isempty(samplerate),ideconspeed==0)
        for ii=1:nseis
            origrate = get(w(ii),'freq');
            if samplerate ~= round(origrate)
                np = origrate/samplerate;
                % when original sampling rate is not an integer (numerical issue)
                if ~(round(origrate) == origrate)
                    if (abs(round(np)-np) < 1e-3) np = round(np); end
                end
                if mod(np,1)==0     % check if integer
                    % see waveform/resample
                    w(ii) = resample(w(ii),'builtin',np);
                else
                   disp('WARNING -- RESAMPLING AT NON-INTEGER INCREMENTS OF ORIGINAL SAMPLES:');
                   %error('resampling must be an integer');
                   dtold = 1/origrate/spdy;
                   dtnew = 1/samplerate/spdy;
                   txi = get(w(ii),'start');
                   txf = get(w(ii),'end') - dtold;
                   tnew = [txi:dtnew:txf]';
                   dxi = get(w(ii),'data');
                   xmethod = 'linear';
                   dnew = interp1(get(w(ii),'timevector'),get(w(ii),'data'),tnew,xmethod);

                   if samplerate > origrate
                       % so what appears to be the Nyquist is actually
                       % based on an interpolated time series; the actual
                       % Nyquist is lower than implied by the sample rate
                       warning('you are upsampling your waveforms from %i sps to %i sps',origrate,samplerate); 
                       %error('you probably do not want to upsample your waveforms from %i sps to %i sps',origrate,samplerate); 
                   end
                   
                   w(ii) = set(w(ii),'data',dnew);
                   w(ii) = addhistory(w(ii),sprintf('Resampled as %s with interp1: %.3f',xmethod,np));
                   w(ii) = set(w(ii),'freq',samplerate);
                end
                disp(sprintf('  %s %s orig rate %i sps, new rate %i sps, np = %.2f',...
                       get(w(ii),'station'),get(w(ii),'channel'),origrate,samplerate,np));
            end
        end
    end  
  
    % WARNING: In some cases of deconvolution, tapering ('taper') is needed
    %          prior to integration; iint here is mainly useful for when
    %          writing sac files; in other cases you can do the post-processing
    %          outside of getwaveform.m
    if and(FCUT2_PAR > 0.5,iprocess==2)
      %w = taper(w);
      disp(sprintf('CUTTING THE FIRST %.4f AND LAST %.4f PAR OF WAVEFORMS',NFRAC,NFRAC));
      for ii=1:nseis
          nl = get(w(ii),'data_length');
          i1 = ceil(NFRAC*nl);
          i2 = floor((1-NFRAC)*nl);
          w(ii) = extract(w(ii),'index',i1,i2);
      end
    end

    % integrate, differentiate, or do nothing
    % default waveforms are velocity so iint=1 will give displacement,
    % and iint=-1 will give acceleration
    % NOTE: iint here is mainly useful for when writing sac files;
    %       in other cases you can do the post-processing outside of getwaveform.m
    if iint == 1  
      w = detrend(w);
      w = integrate(w);
    elseif iint == -1
      w = diff(w);
    end

    % add fields to waveform object
    % note: safer to get elon and elat from the site table
    w = addfield(w,'KEVNM',char(etag));
    %w = addfield(w,'NEVID',eid);
    w = addfield(w,'EVLA',elat);
    w = addfield(w,'EVLO',elon);
    w = addfield(w,'EVDP',edep);
    w = addfield(w,'MAG',emag);

    nseis = length(w);
    for ii = 1:nseis
      % site headers
      %w(ii) = addfield(w(ii),'KNETWK',get(w(ii),'network'));
      w(ii) = addfield(w(ii),'STLA',site(ii).lat);
      w(ii) = addfield(w(ii),'STLO',site(ii).lon);
      w(ii) = addfield(w(ii),'STEL',site(ii).elev);
      w(ii) = addfield(w(ii),'CMPAZ',sitechan(ii).hang);
      w(ii) = addfield(w(ii),'CMPINC',sitechan(ii).vang);

      % compute distances and azimuths
      % note 1: this uses the ellipsoid specified by ell0
      % note 2: probably these should be left done in sac (see note below)
      [dist,az] = distance(elat,elon,site(ii).lat,site(ii).lon,ell0);
      baz = azimuth(site(ii).lat,site(ii).lon,elat,elon,ell0);
      gcarc = km2deg(dist);     % assumes distance is on SPHERE
      %disp(sprintf('DIST %.1f AZ %.1f BAZ %.1f GCARC %.1f',dist,az,baz,gcarc));

      w(ii) = addfield(w(ii),'DIST',dist);
      w(ii) = addfield(w(ii),'AZ',az);
      w(ii) = addfield(w(ii),'BAZ',baz);
      w(ii) = addfield(w(ii),'GCARC',gcarc);
    end  

    % apply Hanna and Long (2012) station correction
    %w = apply_HLazcorr(w);

    %-------Save waveform as sac files-------%
    if isavesac==1
      knetwk = get(w,'network');
      if ~iscell(knetwk)    % convert knetwk into a cell array
         knetwk = {knetwk};
      end

      % savesac will apply all sac headers that are attached to w as added fields
      % note: if fname is not specified, a default file name will be applied
      for ii = 1:length(w)
          scn{ii} = strcat(get(w(ii),'station'),'.',get(w(ii),'channel'),'.',get(w(ii),'network'));
          fname = strcat(etag,'.',scn{ii},'.sac');
          savesac(w(ii),evtdir,fname);
          %savesac(w(ii),evtdir);
      end      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % 2. Add origin time header info into sac files
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % NOTE: If LCALDA = TRUE, then the four sac headers
      % DIST, AZ, BAZ, GCARC will be automatically updated as soon as STLA,
      % STLO, EVLA, EVLO are added. BUT THIS ONLY OCCURS IF YOU ARE IN SAC,
      % NOT WITH THE MATLAB SCRIPTS.
      % Here we compute the four fields in Matlab, but ideally this would be
      % left to sac. Without adding the fields, one can go into sac and type
      % 'w over' -- this will automatically fill the four header fields.
      % It would be best to leave these calculations to sac, since they
      % depend on the assumed ellipsoid for the Earth (see comparison below).
      % Question: What ellipsoid does sac assume?

      % add a marker for the origin time
      % (see also @waveform/private/waveform2sacheader.m)
      if botime
          [nzyr,nzmo,nzdy,nzhr,nzmn,nsec] = datevec(originTime);
          nzjdy = floor(datenum(nzyr,nzmo,nzdy,nzhr,nzmn,nsec) - datenum(nzyr-1,12,31,0,0,0));
          nzsc = fix(nsec);
          nzms = (nsec - fix(nsec)) * 1000;

          % start and end times relative to origin time, in seconds
          tstart = (startTime - originTime)*spdy;
          tend = (endTime - originTime)*spdy;

          for ii = 1:nseis
              % read the file that has the default file name (savesac.m above)
              fstring = strcat(evtdir,'*.',scn{ii},'*');
              %fstring = strcat(evtdir,datestr(startTime,'yyyymmdd_'),'*.',stag);
              sacevt  = dir(fstring);
              sacfile = strcat(evtdir,sacevt.name);
              sacvals = rsac(sacfile);

              sacvals = ch(sacvals,'NZYEAR',nzyr,'NZJDAY',nzjdy,'NZHOUR',nzhr,'NZMIN',nzmn,'NZSEC',nzsc,'NZMSEC',nzms);
              sacvals = ch(sacvals,'B',tstart,'E',tend,'O',0);

              % write sac file with new headers
              % note: this will write over the same file
              ofile = char(strcat(evtdir,etag,'.',scn{ii},'.sac'));
              wsac(ofile,sacvals);
              % delete original sac file (from savesac.m)
              %delete(sacfile);
          end
      end

      % saving station name (needed as input another code)--added by vipul
      % note: what about stations like MCK with two different networks?

%      fid1 = fopen(strcat(evtdir,'st_rec.dat'),'w');    % original
%      st = unique(get(w,'station'));                    % original

      % sort by distance
      fid1 = fopen(strcat(evtdir,'st_rec_dist.dat'),'w');
      [~,iindex] = unique(get(w,'dist'));   % order by 'dist' or 'AZ'
      st = get(w(iindex(:)),'station');  % list stations (in order already)
      st = cellstr(st);
      for jj=1:length(st)
         fprintf(fid1,'%s\n',st{jj});
      end
      fclose(fid1);

      % repeat but sort by azimuth
      fid1 = fopen(strcat(evtdir,'st_rec_az.dat'),'w');
      [~,iindex] = unique(get(w,'AZ'));   % order by 'dist' or 'AZ'
      st = get(w(iindex(:)),'station');  % list stations (in order already)
      st = cellstr(st);
      for jj=1:length(st)
         fprintf(fid1,'%s\n',st{jj});
      end

      fclose(fid1);

    end   % if isavesac

else
    disp('No waveform retrieved');  % w is empty
    return
end

% 20090407_201215.KDAK.BHZ.IR.sac
% sac computations:    
%        DIST = 4.388422e+02
%       GCARC = 3.951996e+00
%          AZ = 2.026800e+02
%         BAZ = 2.022986e+01
%       GCARC = 3.951996e+00
% matlab computations:
% saclst dist gcarc az baz f 20090407_201215.KDAK.BHZ.IR.sac
%                sac: 438.842    3.952       202.68      20.2299
%   ellipsoid, wgs84: 438.841    3.94659     202.682     20.2278
%             sphere: 437.91     3.93822     202.649     20.1943

% check the channels
disp('input channels requested:'); chan
wchan = get(w,'channel');
if ~iscell(wchan), wchan={wchan}; end
uwchan = unique(wchan);
disp('checking how many waveforms for each channel are returned:');
for ii=1:length(uwchan)
    disp(sprintf('channel %6s : %3i',uwchan{ii},sum(strcmp(uwchan{ii},wchan))));
end

%disp(sprintf('%.1f s to execute getwaveform.m',toc));
disp('==================== leaving getwaveform.m ====================');

%==========================================================================

function index2 = soutside(site1,s1,fid)

global ell0 elon elat stasub

% indices of stations that are outside the designated bounding box
index2 = [];

if ~isempty(stasub)
    
    if length(stasub)==2
        % stations within a range of distances from the source
        dmin_km = stasub(1);
        dmax_km = stasub(2);
        
        for ii = 1:length(s1)
            % note 1: uses ell0 ellipsoid specified above
            % note 2: we might get the distance from an added DIST field
            dtemp = distance(site1(ii).lat,site1(ii).lon,elat,elon,ell0);
            if or(dtemp < dmin_km,dtemp > dmax_km)
                if dtemp < dmin_km
                    stext = sprintf('Station %4s is %.2f km (< %.2f km) from the source\n',s1{ii},dtemp,dmin_km);
                else
                    stext = sprintf('Station %4s is %.2f km (> %.2f km) from the source\n',s1{ii},dtemp,dmax_km);
                end
                index2 = [index2 ii];
                fprintf(stext); 
                fprintf(fid,stext);
            end
        end
        
    elseif length(stasub)==4
        % bounding box for stations
        lonmin = wrapTo360(stasub(1));
        lonmax = wrapTo360(stasub(2));
        latmin = stasub(3);
        latmax = stasub(4);

        % regions are labeled from 1 to 8, starting clockwise from NORTH
        stregions = {'north','northeast','east','southeast','south','southwest','west','northwest'};

        for ii = 1:length(s1)
            rlon = wrapTo360(site1(ii).lon);
            if site1(ii).lat < latmin || site1(ii).lat > latmax || rlon < lonmin || rlon > lonmax
                isouth = site1(ii).lat < latmin;
                inorth = site1(ii).lat > latmax;
                iwest  = rlon < lonmin;
                ieast  = rlon > lonmax;
                if inorth==1
                    if iwest==1, ir=8; elseif ieast==1, ir=2; else ir = 1; end
                elseif isouth==1
                    if iwest==1, ir=6; elseif ieast==1, ir=4; else ir = 5; end
                elseif iwest==1
                    ir=7;
                else
                    ir=3;
                end
                stext = sprintf('Station %4s is %s of study region\n',s1{ii},stregions{ir});
                fprintf(stext); 
                fprintf(fid,stext);
                index2 = [index2 ii];
            end
        end
        
    elseif length(stasub)==6
        clon = stasub(1);
        clat = stasub(2);
        dmin_deg = stasub(3);
        dmax_deg = stasub(4);
        azmin = stasub(5);
        azmax = stasub(6);
        stcen = sprintf('%.2f, %.2f',clon,clat);
        
        for ii = 1:length(s1)
            % note 1: uses ell0 ellipsoid specified above
            % note 2: we might get the distance from an added DIST field
            [dtemp,az] = distance(clat,clon,site1(ii).lat,site1(ii).lon,ell0);
            gcarc = km2deg(dtemp);     % assumes distance is on SPHERE
            
            if or(gcarc < dmin_deg,gcarc > dmax_deg)
                if gcarc < dmin_deg
                    stext = sprintf('Station %4s is %.2f deg (< %.2f deg) from %s\n',s1{ii},gcarc,dmin_deg,stcen);
                else
                    stext = sprintf('Station %4s is %.2f deg (> %.2f deg) from %s\n',s1{ii},gcarc,dmax_deg,stcen);
                end
                index2 = [index2 ii];
                fprintf(stext); 
                fprintf(fid,stext);
            elseif or(az < azmin,az > azmax)
                if az < azmin
                    stext = sprintf('Station %4s azimuth is %.2f deg (< %.2f deg) from %s\n',s1{ii},az,azmin,stcen);
                else
                    stext = sprintf('Station %4s azimuth is %.2f deg (> %.2f deg) from %s\n',s1{ii},az,azmax,stcen);
                end
                index2 = [index2 ii];
                fprintf(stext); 
                fprintf(fid,stext);
            end
        end
        
    end
end

%==========================================================================

function [w,s,cha,site,sitechan,iflag] = getw(idatabase,scnl,iflag,fid,...
    startTime,endTime,iprocess,cutoff,samplerate)

global FCUT1_PAR FCUT2_PAR spdy npoles ideconspeed

w = []; % initialize
s = {}; 
cha = {};
site = [];
sitechan = [];
index1 = [];  % catch stations without site info
index4 = [];  % catch stations without sitechan info  
index5 = [];  % catch stations without a calib record in the calibration table during target time range 
index6 = [];  % catch records with error in applying instrument response

[dsta,db,ds,Net,basename] = db_index_test(idatabase);

disp(sprintf('getw: entering %s database',basename));

disp('calling waveform(ds,scnl,startTime,endTime)');
%ds
get(scnl,'station'), get(scnl,'channel'), get(scnl,'network'), get(scnl,'location')
datestr(startTime), datestr(endTime)

if 0==1
    % if there were no database errors, this would work fine
    w = waveform(ds,scnl,startTime,endTime);
    
else
    % the workaround work solution will only work for single channels,
    % so here we loop over the input channels
    chanex = get(scnl,'channel');
    if ~iscell(chanex), chanex = {chanex}; end
    w = [];
    for ii=1:length(chanex)
        scnlx = scnlobject('*',chanex(ii),'','')
        ds
        startTime
        endTime
        try
            w0 = waveform(ds,scnlx,startTime,endTime)
        %--temporary solution to go around "trload_css" failure for uaf_continuous database--
            if isempty(w) && strcmp(basename,'uaf_continuous')
               disp('*****"trload_css" failed, do the workaround*****');
               w0 = waveform(ds,scnlx,startTime,endTime,true); 
            end
        %------------------------------------------------------------------------------------
        catch
        end
        w = [w w0];
    end
end

%-------Instrument response frequency range--------%
% used when iprocess = 2
% used in "response_get_from_db", not in "response_apply"
if iprocess == 2
    frequencies = logspace(-3,3,1000);  % 0.001 - 1000Hz     
end

% The specified startTime could be much different from the start time of
% the waveforms, since the waveforms are only present over a particular
% time interval. To query the database some absolute time is needed --
% this is represented by tpick (search below).

if ~isempty(w)
    w = fillgaps(w,'meanAll');
    s = get(w,'station'); % read station list
    cha = get(w,'channel');
    if ~iscell(s) % When there is only one station, convert s1 into a cell array
       s = {s};
       cha = {cha};
    end
     if ~(strcmp(basename,'uaf_continuous'))
       % note: network info used when calling savesac.m
       if Net ~= '--'
          w = set(w,'network',Net);
          w = addfield(w,'KNETWK',Net);
       else
          for ii = 1:length(s)
              % retrieve network info
              %try
                 Net = db_get_net_info(s(ii),dsta);
              %catch
              %   error('error retreiving network for %s',s{ii}); 
              %end
              Net = char(Net);
              w(ii) = set(w(ii),'network',Net);
              w(ii) = addfield(w(ii),'KNETWK',Net);
          end
       end
     end
    for ii = 1:length(s)
        %tpick = startTime;
        tpick = get(w(ii),'start');
        if (iflag == 0 && iprocess ~= 2)
           disp('station         database    date            site_info');
           fprintf(fid,'%s\n','station         database    date            site_info');
           iflag = 1; 
        end    
        try 
           % read site table from database 
           site = [site db_get_site_info(s{ii},tpick,dsta)];
           if iprocess ~=2
              fprintf('%7s %16s %s %13s\n',s{ii},basename,datestr(tpick,29),'YES'); 
              fprintf(fid,'%7s %16s %s %13s\n',s{ii},basename,datestr(tpick,29),'YES');
           end
        catch
             if iprocess ~=2
                fprintf('%7s %16s %s %13s\n',s{ii},basename,datestr(tpick,29),'NO'); 
                fprintf(fid,'%7s %16s %s %13s\n',s{ii},basename,datestr(tpick,29),'NO');
             end
             index1 = [index1 ii];
        end
    end
    % delete records without site info
    w(index1) = [];
    s(index1) = [];
    cha(index1) = [];
end

% delete records outside the target region
if ~isempty(s)
    index2 = soutside(site,s,fid);
    w(index2) = [];
    s(index2) = [];
    site(index2) = [];
    cha(index2) = [];
end

if (strcmp(basename,'uaf_continuous'))
   if ~isempty(s)
       for ii = 1:length(s)
%            % Search for network name for uaf_continous stations in uafnet 
%            ind = find(strcmp(s(ii),stnm));
%            if isempty(ind)
%               % Catch stations within study region but not listed in uafnet file for some reason
%               index3 = [index3,i]; 
%               fprintf('%s %s %s\n','Station ',s{ii},'is not listed in uafnet');
%               fprintf(fid,'%s %s %s\n','Station ',s{ii},'is not listed uafnet');
%            elseif length(ind) == 1
%               Net = char(netwk(ind)); 
% %          disp(s{ii})
% %          disp(netwk(ind));
%               w(ii) = set(w(ii),'network',Net);
%               w(ii) = addfield(w(ii),'KNETWK',Net);
%            else
%              fprintf('%s %s\n','ERROR: more than one record found for station',s{ii});
%              fprintf(fid,'%s %s\n','ERROR: more than one record found for station',s{ii});
%            end
            %try
                Net = db_get_net_info(s(ii),dsta); % retrieve network info
            %catch
            %    error('error retreiving network for %s',s{ii});
            %end
            Net = char(Net);
            w(ii) = set(w(ii),'network',Net);
            w(ii) = addfield(w(ii),'KNETWK',Net);
        end
%        w(index3) = []; % Delete records not listed in uafnet file
%        s(index3) = [];
%        cha(index3) = [];
%        site(index3) = [];
        
   end
end  

if ~isempty(s)
    for ii = 1:length(s)
        %tpick = startTime;
        tpick = get(w(ii),'start');
        try 
           % read sitechan table from database
           sitechan = [sitechan db_get_sitechan_info(s{ii},cha{ii},datestr(tpick,29),dsta)];
            
        catch
           if iprocess ~= 2
              fprintf('db_get_sitechan_info(''%s'',''%s'',''%s'',''%s'')\n',...
                  s{ii},cha{ii},datestr(tpick,29),dsta);
              fprintf('%7s %16s %s %13s\n',s{ii},basename,datestr(tpick,29),'No sitechan'); 
              fprintf(fid,'%7s %16s %s %13s\n',s{ii},basename,datestr(tpick,29),'No sitechan');
             
           end 
           index4 = [index4 ii];
           
        end
    end
    % delete records without sitechan info
    w(index4) = [];
    s(index4) = [];
    cha(index4) = [];
    site(index4) = [];
end

if and(~isempty(s),iprocess == 0)
    w = remove_calib(w);
end

%----------------CALIBRATION CORRECTION------------------------------------
% Waveforms retrieved by "GISMO/waveform" use calibration numbers stored
% in "wfdisc" for speed purpose; However, sometimes the calibration number
% in "wfdisc" is not found or correct. So waveforms need to be corrected
% according the "calibration" numbers in "calibration" table.
if and(~isempty(s),iprocess ~=0) 
    if ~(strcmp(basename,'uaf_continuous'))
        dbp = dbopen(db,'r'); 
    end
    for ii = 1:length(s)
          %tpick = startTime;
          tpick = get(w(ii),'start');
          if (strcmp(basename,'uaf_continuous'))
              %if (strcmp(get(w(ii),'KNETWK'),'AV'))
              %    dbp = dbopen(dbav,'r');
              %else
                  dbp = dbopen(dsta,'r');
              %end
          end
          dbp = dblookup_table(dbp,'calibration');
          dbp = dbsubset(dbp, sprintf('sta=="%s" && chan=="%s" && calibration.time<=%f &&  calibration.endtime>=%f',...
              s{ii},cha{ii},datenum2epoch(tpick),datenum2epoch(tpick)));
          try 
             calib1 = dbgetv(dbp,'calib');      % 'official' calib that resides in response tables
             calib2 = get(w(ii),'calib');       % 'unofficial' calib that resides in wfdisc
             if (abs(calib1-calib2) > 1e-03)    % calib1 is not equal to calib2
             %if abs(log(calib1/calib2)) > 1e-3
                disp('------------------------------------------------------');
                if (calib2 == 0)
                   fprintf('WARNING: No calib found for %s in wfdisc (%s)\n',s{ii},datestr(tpick));
                   fprintf(fid,'WARNING: No calib found for %s in wfdisc (%s)\n',s{ii},datestr(tpick));
                   calib = calib1;
                else 
                   fprintf('WARNING: Wrong calib for %s in wfdisc (%s)\n',s{ii},datestr(tpick));
                   fprintf(fid,'WARNING: Wrong calib for %s in wfdisc (%s)\n',s{ii},datestr(tpick)); 
                   calib = calib1/calib2;
                end        
                fprintf('Apply calibration correction\n');  
                fprintf(fid,'Apply calibration correction\n');  
                newData = get(w(ii),'DATA') * calib;
                w(ii) = set(w(ii),'DATA',newData);
                w(ii) = addhistory(w(ii),'Calibration correction applied');
                w(ii) = set(w(ii),'calib',calib1);
                w(ii) = set(w(ii),'calibration_applied','YES');
             end
             % segtype2units: pull unit out of "calibration" table (segtype field)
             % when unit is missing from "wfdisc"
             if strcmp(get(w(ii),'UNITS'),'null')
                w(ii) = segtype2units(w(ii),dbp);
             end
          catch
             fprintf('WARNING: No calib record for station %s (%s) in calib table (%s)\n',s{ii},cha{ii},datestr(tpick));
             fprintf('Waveform removed \n');  
             fprintf(fid,'WARNING: No calib record for station %s (%s) in calib table (%s)\n',s{ii},cha{ii},datestr(tpick));
             fprintf(fid,'Waveform removed \n');  
             index5 = [index5 ii];
             
%              %% testing
%              tpick = datenum('2012/04/11 09:18:37.00');
%              sta = 'COLA';
%              cha = 'BH1';  % BH2
%              dsta = '/aerun/sum/params/Stations/master_stations';
%              dbp = dbopen(dsta,'r');
%              dbp = dblookup_table(dbp,'calibration');
%              dbp = dbsubset(dbp, sprintf('sta=="%s" && chan=="%s" && calibration.time<=%f &&  calibration.endtime>=%f',...
%                 sta,cha,datenum2epoch(tpick),datenum2epoch(tpick)));
%              calib1 = dbgetv(dbp,'calib')
%              %%
          end
    end
    w(index5) = [];  % Delete records without calib
    s(index5) = [];
    cha(index5) = [];
    site(index5) = [];
    sitechan(index5) = [];
    dbclose(dbp);
end    
%--------------------------------CALIBRATION CORRECTION END----------------
%    
%-----------------------------INSTRUMENT RESPONSE DECONVOLUTION------------

%--------Filter applied during instrument response deconvolution-------%
% used when iprocess = 2
% make sure cutoff frequncies be within (0,1) after normalized by Nyquist frequency
% Nyquist frequency = sampling rate / 2
% Nyquist = get(w,'NYQ')
% freq    = get(w,'FREQ')
% For example, if you are looking at broadband data, BEAAR MOOS ARCTIC have
% sampling rate at 50 HZ and AEIC has varying sampling rate at 20HZ, 40HZ,
% 50HZ, 100HZ, then the low pass frequency f2 in filterobject would be less
% than 10 Hz.
% 'B' & 'npoles = 2' in "filterobject" leads to 4th order Butterworth
% bandpass filter within "response_apply".
% When a very broad freqeuency range ("cutoff") is chosen (e.g. [0.001 4.5]),
% npoles = 2 and 4 make nearly no difference to output waveforms except a
% small segment at the end; however a narrow frequency band is chosen
% (e.g. [0.01 0.05]), npoles = 2 is preferred; npole = 4 makes filter
% unstable, which leads to 0 amplitude in output waveforms.

if ~isempty(w)
   if iprocess == 2
      if and(~isempty(samplerate),ideconspeed==1)
         Nyquist = samplerate/2; 
      else
         % THIS MAY CAUSE PROBLEMS SINCE IT TIES cutoff TO THE RECORD(S)
         % WITH THE MINIMUM SAMPLE RATE, BUT NOTE THAT THE NYQUIST IS ONLY
         % USED TO DECIDE ON THE CUTOFF FREQUENCIES.
         %get(w,'NYQ')
         [Nyquist,imin] = min(get(w,'NYQ'));
         fprintf('\nNyquist is set to the overall minimum: %.2f Hz at %s (and others)\n',Nyquist,get(w(imin),'station'));
      end
      duration = (endTime - startTime)*spdy;
      
      if ~isempty(cutoff)
         %if (1/cutoff(1) > duration/FCUT1_PAR)
         if (cutoff(1) < FCUT1_PAR/duration)
            fprintf('\n*********************************************');
            fprintf('\nERROR: High pass frequency in cutoff needs to be greater than %.2f/(length of waveform) = %f Hz',FCUT1_PAR,duration);
            fprintf('\n length of waveform = %.2f s, f1 = %.4f Hz, 1/f1 = %.2f s',duration,cutoff(1),1/cutoff(1));
            fprintf('\n So within the time window there are at least %.2f wavelengths of the longest period waveform\n',FCUT1_PAR);
            error('exit here');
         else
             fprintf('High pass frequency in cutoff = %.3f Hz\n',cutoff(1));
         end
         if (cutoff(2) > Nyquist*FCUT2_PAR)
            fprintf('\n*********************************************');
            fprintf('\nERROR: Low pass frequency in cutoff needs to be <= than %.2f of the Nyquist frequency %f Hz',FCUT2_PAR,Nyquist);
            fprintf('\n fNyq = %.3f Hz, fNyq * %.2f = %.3f Hz which is < f2 = %.4f Hz\n',Nyquist,FCUT2_PAR,Nyquist*FCUT2_PAR,cutoff(2));
            error('exit here');
         else
            fprintf('Low pass frequency in cutoff = %.3f Hz\n',cutoff(2));
         end
         filterObj = filterobject('B',[cutoff(1) cutoff(2)],npoles);  
      else
         % set default cutoff frequencies for deconvolution
         fprintf('\ncutoff frequencies not specified for deconvolution');
         fprintf('\n  FCUT1_PAR = %.3e Hz',FCUT1_PAR);
         fprintf('\n  duration = %.3e Hz',duration);
         fprintf('\n  Nyquist = %.3e Hz',Nyquist);
         fprintf('\n  FCUT2_PAR = %.3e Hz',FCUT2_PAR);
         fprintf('\n  --> default %i poles and [FCUT1_PAR/duration Nyquist*FCUT2_PAR] = [%.3e %.3e] Hz\n',...
             npoles,FCUT1_PAR/duration,Nyquist*FCUT2_PAR);
         filterObj = filterobject('B',[FCUT1_PAR/duration Nyquist*FCUT2_PAR],npoles);
      end
   end
end

if ~isempty(s) && iprocess == 2
    if iflag == 0
       disp('station  channel       database       response_file');
       fprintf(fid,'%s\n','station  channel       database       response_file');
       iflag = 1;
    end
    
    % demean and detrend before filtering
    w = demean(w);
    w = detrend(w);
    
    for ii = 1:length(s)
        %tpick = startTime;
        tpick = get(w(ii),'start');
        
%         if (strcmp(basename,'uaf_continuous'))
%            %if (strcmp(get(w(ii),'KNETWK'),'AV'))
%            %   dsta = dbav;
%            %else
%               dsta = '/aerun/sum/params/Stations/master_stations'; 
%            %end
%         end
    
        %------------------------------------------------------------------ 
        % resample waveform in order to speed up deconvolution
        % WARNING: THIS IS ONLY ACCURATE IF YOU'RE FAR FROM THE NYQUIST; NEAR
        % THE NYQUIST, THE INSTRUMENT RESPONSE TAKES INTO ACCOUNT THE ORIGINAL
        % SAMPLING RATE.
        if and(~isempty(samplerate),ideconspeed==1)
            disp('RESAMPLING TO SPEED UP DECONVOLUTION');
            origrate = get(w(ii),'freq');
            np = origrate/samplerate;
            % when original sampling rate is not an integer (numerical issue)
            if ~(round(origrate) == origrate)
                if (abs(round(np)-np) < 1e-3), np = round(np); end
            end
            if mod(np,1)==0     % check if integer
                % see waveform/resample
                w(ii) = resample(w(ii),'builtin',np);
            else
               error('resampling interval must be an integer');
            end
        end   
        %------------------------------------------------------------------ 
        try
           response(ii) = response_get_from_db(s{ii},cha{ii},tpick,frequencies,dsta);
           [~,resfile,ext] = fileparts(response(ii).respFile);
           resfile = [resfile,ext];
           fprintf('%7s %7s %16s %s   %s\n',s{ii},cha{ii},basename,datestr(tpick,29),resfile);
           fprintf(fid,'%7s %7s %16s %s   %s\n',s{ii},cha{ii},basename,datestr(tpick,29),resfile);
           
           % add name of response file
           w(ii) = addfield(w(ii),'RESPFILE',resfile);
           
           try
               w(ii) = response_apply(w(ii),filterObj,'antelope',dsta); 
           catch 
              %dsta = '/aerun/sum/params/Stations/master_stations';
              %w = response_apply(w,filterobject('B',[0.01 5],2),'antelope',dsta);
              fprintf('%25s %15s\n','','Error in response_apply');
              fprintf(fid,'%25s %15s\n','','Error in response_apply');
              index6 = [index6 ii];
           end 
        catch
           fprintf('%7s %7s %16s %s %s\n',s{ii},cha{ii},basename,'  ','Not found');
           fprintf(fid,'%7s %7s %16s %s %s\n',s{ii},cha{ii},basename,'  ','Not found');
           index6 = [index6 ii];
        end
    end
    
    w(index6) = [];  % Delete records with error in applying instrument response
    s(index6) = [];
    cha(index6) = [];
    site(index6) = [];
    sitechan(index6) = [];
end    

%----------------INSTRUMENT RESPONSE DECONVOLUTION END---------------------

function woutput = segtype2units(winput,dbp)
% PULL UNIT OUT OF "CALIBRATION" TABLE (SEGTYPE FIELD) WHEN IT IS MISSING FROM "WFDISC"
%'segtype' in antelope datasets indicate the natural units of the detector
segtype = dbgetv(dbp,'segtype');
segTypes = 'ABDHIJKMPRSTVWabcdfhimnoprstuvw-';
segUnits = {'A','nm / sec / sec','acceleration';
  'B', '25 mw / m / m','UV (sunburn) index(NOAA)';
  'D', 'nm', 'displacement';
  'H','Pa','hydroacoustic';
  'I','Pa','infrasound';
  'J','watts','power (Joulses/sec) (UCSD)';
  'K','kPa','generic pressure (UCSB)';
  'M','mm','Wood-Anderson drum recorder';
  'P','mb','barometric pressure';
  'R','mm','rain fall (UCSD)';
  'S','nm / m','strain';
  'T','sec','time';
  'V','nm / sec','velocity';
  'W','watts / m / m', 'insolation';
  'a','deg', 'azimuth'
  'b','bits/ sec', 'bit rate';
  'c','counts', 'dimensionless integer';
  'd','m', 'depth or height (e.g., water)';
  'f','micromoles / sec / m /m', 'photoactive radiation flux';
  'h','pH','hydrogen ion concentration';
  'i','amp','electric current'
  'm','bitmap','dimensionless bitmap';
  'n','nanoradians','angle (tilt)';
  'o','mg/l','diliution of oxygen (Mark VanScoy)';
  'p','percent','percentage';
  'r','in','rainfall (UCSD)';
  's','m / sec', 'speed (e.g., wind)';
  't','C','temperature';
  'u','microsiemens/cm','conductivity';
  'v','volts','electric potential';
  'w','rad / sec', 'rotation rate';
%  '-','null','null'};
  ' ','null','null'};
if isempty(segtype)
  segtype=  '-';
end
if ~ismember(segtype,segTypes)
  segtype=  '-';
end
thisseg = find(segtype==segTypes);
units = segUnits{thisseg,2};
%type_of_data = segUnits{thisseg,3};   
woutput = set(winput,'UNITS', units);

%--PULL UNIT FROM CALIBRATION TABLE WHEN IT IS MISSING FROM WFDISC END-----

%==========================================================================
