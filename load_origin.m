function origin = load_origin(dbpath, subset_expr)
 %LOAD ORIGIN
   if ~exist(dbpath, 'file')
       warning(sprintf('dbpath %s does not exist', dbpath));
   end
   db = dbopen(dbpath, 'r');
   db = dblookup_table(db, 'origin');
   %db = dbsubset(db, subset_expr); %to subset data - i.e. by magnitude, etc.
   %db = dbjoin(db, dbe);
   %db = dbsubset(db, 'orid==prefor'); %to load in only the preferred origins for each event
   if nargin>1
       db = dbsubset(db, subset_expr);
   end
   db = dbsort(db, 'time');
   [lat, lon, depth, time, ml, orid, evid] = dbgetv(db, 'lat', 'lon', 'depth', 'time', 'ml', 'orid', 'evid');
   dbclose(db);
   origin = struct('lat', lat, 'lon', lon, 'depth', depth, 'time', time, 'ml', ml, 'orid', orid, 'evid', evid);
end