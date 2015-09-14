function [ParrivalStruct, SarrivalStruct] = loadArrival(databasePath, expr)
% Write a function to load the P arrival times from an arrival table into a structure.
% (Hint: base this on the loadOriginTable function above. Also you don't need to open any other tables - like assoc or event).
% 
%   loadOriginTable  Load an origin table
%   originStruct = loadOriginTable(databasePath) loads an origin table from
%   the database indicated by databasePath.
%
%   originStruct = loadOriginTable(databasePath, subsetExpression) applies
%   dbsubset(subsetExpression) before reading in the origin table records.
%
%   Inputs:
%       databasePath - the path to the database descriptor
%       subsetExpression - (optional) a dbsubset expression to apply
%   Output:
%       originStruct - a structure containing the fields lat, lon, depth,
%       time, ml, orid and evid. Each of these is a vector of length equal
%       to the number of records
%
%   Author: Glenn Thompson 2014/10/01; Modified by: Alexandra Farrell
%   2014/10/15-2015/10/16 and 2014/11/12

    %arrivalStruct = struct('time', []); %Create arrivalStruct structure to be filled in later   
    
    if ~exist(databasePath, 'file') %Check to make sure that the database exists
        warning(sprintf('Database %s does not exist',databasePath));
        return
    end  
    
    db = dbopen(databasePath, 'r'); %Open to read database at specified databasePath
    db = dblookup_table(db, 'arrival'); %Look up arrival table, assign that to variable db   
    dbe = dblookup_table(db, 'origin');
    dbe2 = dblookup_table(db, 'assoc');
    dbe = dbjoin(dbe, dbe2); %join origin and assoc tables
    db = dbjoin(dbe, db); %join ^ and arrival tables
    db2 = dbsubset(db, 'phase==''S''');  %subset the database in variable db to preferred orid's, assign to db2
    db = dbsubset(db, 'phase==''P''');   %subset the database in variable db to the P phase, reassign to db
     
    if nargin>1 %check to see if a string subsetExpression was given; if so, db is assigned as the given subset
        db = dbsubset(db, expr);
        db2 = dbsubset(db2, expr);
    end    
    
    db = dbsort(db, 'time');   %sort database by time
    db2 = dbsort(db2, 'time'); %sort database by time
    [lat, lon, depth, time, ml, orid, evid, sta] = ...
    dbgetv(db, 'lat', 'lon', 'depth', 'time', 'ml', 'orid', 'evid', 'sta'); %assign variables to each column specified    
    for n = 1:numel(lat)    
        ParrivalStruct(n) = struct('lat', lat(n), 'lon', lon(n), 'depth', depth(n), 'time', time(n), 'ml', ml(n), 'orid', orid(n), 'evid', evid(n), 'sta', sta(n));
    end
    
    [lat, lon, depth, time, ml, orid, evid, sta] = ...
    dbgetv(db2, 'lat', 'lon', 'depth', 'time', 'ml', 'orid', 'evid', 'sta');
    for i = 1:numel(lat)
        SarrivalStruct(i) = struct('lat', lat(i), 'lon', lon(i), 'depth', depth(i), 'time', time(i), 'ml', ml(i), 'orid', orid(i), 'evid', evid(i), 'sta', sta(i));
    end
    dbclose(db)