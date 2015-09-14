function [dsta,db,ds,net,bname] = db_index(ind)
%DB_INDEX given index, return database path for databases at UAF/GI
%
% INPUT
%   ind     index of database
%
% OUTPUT
%   dsta    path to stations
%   db      path to waveforms
%   ds      datasource class
%   net     network
%   bname   base name for database
%
% EXAMPLE: [dsta,db,ds,net,bname] = db_index(2)
%

NDB = 8;

% check input
if ~any(ind==[1:NDB]), error('database index (%i) must be 1-%i',ind,NDB); end

% waveform database
% uaf_continuous --> '/aerun/sum/db/archive/archive_[YEAR]/archive_[YEAR]_[MONTH]_[DAY]'
db1 = 'uaf_continuous'; 
db2 = '/home/admin/databases/BEAAR/wf/beaar';
db3 = '/home/admin/databases/ARCTIC/wf/arctic';
db4 = '/home/admin/databases/MOOS/wf/moos';
db5 = '/home/admin/databases/YAHTSE/wf/yahtse';
%db6 = '/home/admin/databases/PLUTONS/wf/plutons';
db6 = '/raid/data/antelope/databases/PLUTONS/dbmerged';
db7 = '/home/admin/databases/SUMATRA/data/sumatra';
db8 = '/home/admin/databases/1995_MFSZ/wf/1995_msfz';  % note typo mfsz

ds1 = datasource(db1);
ds2 = datasource('antelope',db2);
ds3 = datasource('antelope',db3);
ds4 = datasource('antelope',db4);
ds5 = datasource('antelope',db5);
ds6 = datasource('antelope',db6);
ds7 = datasource('antelope',db7);
ds8 = datasource('antelope',db8);

% networks
% net='--'; means that network will be retrieved by db_get_snetsta_info.m 
nets = {'--','XE','XR','YV','XF','XP','--','--'};
net = nets{ind};

% base names for databases
bnames = {'uaf_continuous','BEAAR','ARCTIC','MOOS','YAHTSE','PLUTONS',...
    'SUMATRA','1995_MFSZ'};
bname = bnames{ind};

switch ind
    case 1, db = db1; ds = ds1;
        % Station database for 'uaf_continous' (for general purpose) 
        % note: to check date of last update: ls -ltr /aerun/sum/params/Stations/
        dsta = '/aerun/sum/params/Stations/master_stations';
    case 2, db = db2; ds = ds2; dsta = db;
    case 3, db = db3; ds = ds3; dsta = db;
    case 4, db = db4; ds = ds4; dsta = db;
    case 5, db = db5; ds = ds5; dsta = db;
    case 6, db = db6; ds = ds6; dsta = db;
    case 7, db = db7; ds = ds7; dsta = db;
    case 8, db = db8; ds = ds8; dsta = db;
end

%==========================================================================
