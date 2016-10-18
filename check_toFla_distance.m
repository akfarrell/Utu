close all;
cell_array = {'058A', -81.80, 27.06, '059A', -81.14, 26.97, '059Z', -81.44, 26.34, '060A', -80.36, 27.04, ...
'060Z', -80.56, 26.41, '449A', -87.22, 30.76, '450A', -86.59, 30.80, '451A', -85.75, 30.62, '452A', ....
-85.18, 30.85, '453A', -84.32, 30.85, '454A', -83.63, 30.71, '455A', -83.03, 30.74, '456A', -82.02, 30.72, ...
'457A', -81.56, 30.62, '552A', -85.29, 30.13, '553A', -84.43, 30.19, '554A', -83.68, 30.08, '555A', ...
-82.97, 30.12, '556A', -82.41, 30.00, '557A', -81.73, 30.02, '655A', -83.26, 29.51, '656A', -82.53, 29.37, ...
'657A', -81.87, 29.59, '658A', -81.26, 29.42, '757A', -82.07, 28.94, '758A', -81.20, 28.96, '857A', ...
-82.23, 28.27, '858A', -81.36, 28.21, '859A', -80.90, 28.06, '957A', -82.24, 27.67, '958A', -81.75, 27.59, ...
'959A', -80.88, 27.52};

names = 1:3:numel(cell_array);
long = 2:3:numel(cell_array);
lat = 3:3:numel(cell_array);

cell_names = cell_array(names);
cell_lat = cell_array(lat);
cell_long = cell_array(long);

eq_lat = -69.6;
eq_lon = -17.7;

distances = zeros(1,numel(cell_names))
azimuths = zeros(1,numel(cell_names))

for i = 1:numel(cell_names)
    [distances(i),azimuths(i)] = distance(cell_lat{i}, cell_long{i}, eq_lat, eq_lon);
end

deg_dist = distances;
distances = distances*111.12;

[sort_distances,I] = sort(distances)
names_sort = cell_names(I)
azimuths_sort = azimuths(I)
deg_dist_sort = deg_dist(I)
lat_s = cell_lat(I)
lon_s = cell_long(I)

h = figure;
latlim = [26 31.25]; %[southern_limit northern_limit] 
lonlim = [-88 -79.5]; %[western_limit eastern_limit]
set(h, 'Position', [1000 1000 1000 1000])
worldmap(latlim, lonlim);
setm(gca, 'fontsize',11);
for i = 1:numel(cell_names)
    scatterm(lat_s{i}, lon_s{i}, 'k')
    order = sprintf('%d',i);
    dist = sprintf('%3.2f', deg_dist_sort(i));
    textm(lat_s{i}, lon_s{i}+0.1, order)
    textm(lat_s{i}-0.15, lon_s{i}-0.1, names_sort{i})
    textm(lat_s{i}+0.15, lon_s{i}-0.1, dist)
    clear order
    clear dist
end