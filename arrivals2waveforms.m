function w = arrivals2waveforms(dbpath, arrivals, pretrig, posttrig)
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
%       pretrig - number of seconds before arrival to load
%       posttrig - number of seconds after arrival to load
%   Output:
%       w - 
%
%   Author: Alexandra Farrell 2014/11/12 from work by Glenn Thompson
    ds = datasource('antelope', dbpath); %create datasource containing the path to the database and the data format
    %scnl = scnlobject('*', 'HHZ'); %screate scnlobject containing the station and component
    secs2days = 60*60*24;
    pretrig = pretrig/secs2days;
    posttrig = posttrig/secs2days;
    numel(arrivals)
    for index = 1:numel(arrivals)
        index
        arrivals(index).orid
        snum = epoch2datenum(arrivals(index).time)-pretrig; %define start time as arrival time minus the number of seconds in pretrig
        enum = epoch2datenum(arrivals(index).time)+posttrig; %define end time as arrival time plus the number of seconds in posttrig
        scnl = scnlobject(arrivals(index).sta, 'HHZ', '', '');
        w(index) = waveform(ds, scnl, snum, enum); %create a waveform object for each arrival with the above defined specifications
    end
end