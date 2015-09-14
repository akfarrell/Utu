function siteStruct = loadSiteTable(databasePath)
% loadSiteTable  Load a site table
%   siteStruct = loadSiteTable(databasePath) loads an origin table from
%   the database indicated by databasePath.
%
%   siteStruct = loadSiteTable(databasePath, subsetExpression) applies
%   dbsubset(subsetExpression) before reading in the site table records.
%
%   Inputs:
%       databasePath - the path to the database descriptor
%       subsetExpression - (optional) a dbsubset expression to apply
%   Output:
%       siteStruct - a structure containing the fields sta, ondate,
%       offdate, lat, lon, elev. Each of these is a vector of length equal
%       to the number of records
%
%   Author: Ophelia George 2014/10/21
%% Initialize the structure with pertinent fields/ ensure that the file exist
    siteStruct = struct('sta', [], 'ondate', [], 'offdate', [], 'lat', [], 'lon', [], 'elev', []);
    if ~exist(databasePath, 'file')
        warning(sprintf('Database %s does not exist',databasePath));
        return
    end

%% If the database exists, then open and begin searching for necessary fields.

    db = dbopen(databasePath, 'r'); %open the database
    db = dblookup_table(db, 'site'); %find the origin table
    
    %Subset the new table based on any subset expression supplied by the
    %user
%     if exist(subsetExpression, 'var') & isstr(subsetExpression)
%         db = dbsubset(db, subset_Expression);
%     end
    
    %sort the database from the most southern to the most northern latitudes
    db = dbsort(db, 'lat'); 
    
    %grab the variables needed in the siteStruct from the new database

    [siteStruct.sta, siteStruct.ondate, siteStruct.offdate, siteStruct.lat, siteStruct.lon, siteStruct.elev] = ...
        dbgetv(db, 'sta', 'ondate', 'offdate', 'lat', 'lon', 'elev'); 
    
    %close out the database once all the needed info has been extracted
    dbclose(db); 
    
end



