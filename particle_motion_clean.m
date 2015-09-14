%Cassandra Smith Feb 12 2015

% create a datasource object
    dbpath = '/home/c/cmsmith10/dbpavlof2007_224_273';
    ds = datasource('antelope', dbpath);

    chan = 'EHZ';
    network = 'AV';
    location = '--';
    sta = 'PV6';
        
    epoch_start = 1189058675.820;
    epoch_end   = 1189058676;

    %create our filterobject
    f = filterobject('B', [10 15], 2)
    
    % create a scnl object
    scnl = scnlobject(sta, chan, network, location);

    % get our waveform object
    snum = epoch2datenum(epoch_start);
    enum = epoch2datenum(epoch_end);
    w = waveform(ds, scnl, snum, enum);
    
    fwave = filtfilt(f,w)

    dataEHZ = get(fwave, 'data')

    
    chan2 = 'EHN';
   
    % create a scnl object
    scnl2 = scnlobject(sta, chan2, network, location);

    % get our waveform object

    w2 = waveform(ds, scnl2, snum, enum);
    
    fwave2 = filtfilt(f,w2)
    
    dataEHN = get(fwave2, 'data')

    
    chan3 = 'EHE';
   
    % create a scnl object
    scnl3 = scnlobject(sta, chan3, network, location);

    % get our waveform object

    w3 = waveform(ds, scnl3, snum, enum);
    
    fwave3 = filtfilt(f,w3)
    
    dataEHE = get(fwave3, 'data')
    
    plot3(dataEHE, dataEHN, dataEHZ)
    grid on
    box on
    xlabel 'EHE'
    ylabel 'EHN'
    zlabel 'EHZ'