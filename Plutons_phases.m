dbpath = '/raid/data/antelope/databases/PLUTONS/dbmerged';

%Load in arrival table joined to origin/assoc tables, with both S
[Parrival, Sarrival] = loadArrival(dbpath, 'lat>=-22.6698 && lat<=-22.1415 && lon<=-66.7865 && lon>=-67.6381 && depth<=30 && time>=''4/15/2010 (105) 0:0:0.00000'' && orid <= 27000');
lat = Parrival.lat; lon = Parrival.lon; depth = Parrival.depth; time = Parrival.time; ml = Parrival.ml; orid = Parrival.orid; evid = Parrival.evid;


%load in origin table for the hell of it
origin = load_origin(dbpath, 'lat>=-22.6698 && lat<=-22.1415 && lon<=-66.7865 && lon>=-67.6381 && depth<=30 && time>=''4/15/2010 (105) 0:0:0.00000'' && orid <= 27000');

w = arrivals2waveforms(dbpath, Parrival, 5, 20);