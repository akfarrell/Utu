function w = load_waveforms(dbpath, start, posttrig)
% arrivals2waveforms  create waveform data for a time period around each
% arrival
%   w = arrivals2waveforms(dbpath, arrivals, pretrig, posttrig) creates
%   waveform objects for each arrival
%
% w = arrivals2waveforms(dbpath, arrivals, pretrig, posttrig)
%
%   Inputs:
%       dbpath - the path to the database containing your wfdisc table
%       arrivals - the output from loadArrivalTable.
%       start - start time of waveform to load
%       posttrig - number of seconds after start of waveform to load
%   Output:
%       w - waveform object of specified data
%
%   Author: Alexandra Farrell 2015/3/27 from work by Glenn Thompson
    ds = datasource('antelope', dbpath); %create datasource containing the path to the database and the data format
    scnl = scnlobject('*', 'BHZ'); %screate scnlobject containing the station and component
    snum = datenum(start); %starting time
    enum = datenum(start) + posttrig; %ending time
    w = waveform(ds, scnl, snum, enum); %create a waveform object for each arribal with the above defined specification
end