close all 
j = figure
 worldmap('World')
 load coast
 [latcells, loncells] = polysplit(lat, long);
 numel(latcells)
 plotm(lat, long, 'k')
 hold on
 scatterm(eq(9).lat, eq(9).lon, 'r', 'filled')
 scatterm(eq(10).lat, eq(10).lon, 'b', 'filled')
 scatterm(eq(3).lat, eq(3).lon, 'g', 'filled')
 