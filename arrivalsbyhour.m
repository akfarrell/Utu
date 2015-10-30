function arrivalsbyhour(db, epochtimestart, epochtimeend, stepsize)
%ARRIVALSBYHOUR
%  arrivalsbyhour(db, epochtimestart, epochendtime, stepsize)
% arrivalsbyhour('mydb', 1431907200, 1443744000, 60*60)
db
    fout = fopen('cassfile','w');
    for thisepoch = epochtimestart:stepsize:epochtimeend
        subset_expr = sprintf('time >= %f && time < %f', thisepoch, thisepoch+stepsize)
        arrivals = load_arrivals(db, subset_expr);
        if ~isempty(arrivals.arid)
            %w = arrivals2waveforms('/home/c/cmsmith10/Dissertation/dbs/dbSAK', arrivals, 5, 10, 2, 13000);
            w = arrivals2waveforms(db, arrivals, 5, 10, 2, 13000);
            m = max(abs(detrend(w)));
            for count=1:numel(m)
                fprintf(fout,'%d %e\n',arrivals.arid(count), m(count));
            end
        end
    end
    fclose(fout);

end


function arrival = load_arrivals(db, subset_expr)
    % LOAD_ARRIVALS Load arrivals from a CSS3.0 database


    %% Load the earthquake catalogue

    fprintf('Loading arrivals from %s',db);

    % Open database
    db = dbopen(db,'r');

    % Apply subset expression
    db = dblookup_table(db,'arrival');
    if exist('subset_expr','var')
        db = dbsubset(db,subset_expr);
    end
    
    if dbnrecs(db) == 0
        arrival = struct('arid',[]);
        return
    end

    % Get the values
    [sta,chan,time,arid,iphase] = dbgetv(db,'sta','chan','time','arid','iphase');
    if strcmp(class(sta),'char')
       sta = {sta};
       chan = {chan};
       iphase = {iphase};
    end
    %time = epoch2datenum(time);

    % Close database link
    dbclose(db);

    % Display counts
    fprintf('\n%d arrivals\n',numel(time));

    arrival.sta = sta;
    arrival.chan = chan;
    arrival.time = time;
    arrival.arid = arid;
    arrival.iphase = iphase;
end

function w = arrivals2waveforms(dbpath, arrivals, pretrig, posttrig, taper_seconds, nwaveforms)
% ARRIVALS2WAVEFORMS Load waveform objects corresponding to an arrivals structure
%	w = arrivals2waveforms(dbpath, arrivals, pretrig, posttrig)
	w = waveform(); % Initialize output

	% Check input variables. 
	if nargin<4
		disp('Incorrect number of input variables')
		return
	end

	% Create datasource
	ds = datasource('antelope', dbpath);
	% Loop over arrivals structure
	for i=1:min([nwaveforms numel(arrivals.arid)])
        scnl = scnlobject(arrivals.sta{i}, arrivals.chan{i});
		snum = epoch2datenum(arrivals.time(i) - pretrig - taper_seconds);
		enum = epoch2datenum(arrivals.time(i) + posttrig + taper_seconds);
        try
%             arrivals.time(i)
%             ds
%             scnl
%             datestr(snum)
%             datestr(enum)
            thisw = waveform(ds, scnl, snum, enum);
%             plot(w)
%             choice = input(' any key to continue');
        end
        if numel(thisw)==1
            w(i) = thisw;
        else
            w(i)=waveform();
        end
	end
end

