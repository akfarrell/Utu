close all
%r=openfig('lvztop_depths_surfaceplot');
r=openfig('Alex_Fig');
hold on
ch=get(gca,'ch');
x=get(ch,'xd')';
y=get(ch,'yd')';

%lat_sta
%lon_sta


z=get(ch,'zd');
zsize = size(z);
z_mean = nanmean(nanmean(z));
z_max = max(max(z));
z_min = min(min(z));

zmin = z; zmax = z; zmean = z;
for i = 1:zsize(1)
    warning('on')
    for k = 1:zsize(2)
        if isnan(z(i,k))
            zmin(i,k) = z_min; %replace NaN's with min(min(z))
            zmax(i,k) = z_max;
            zmean(i,k) = z_mean;
        end
    end
end
warning('Replacing NaNs with min(z)!!!!')
%warning('Replacing NaNs with zeros!!!!')
warning('off')

for i = 1:numel(lat_sta)
    depths_at_stas_mean(i) = interp2(x,y,zmean,lon_sta(i), lat_sta(i));
    depths_at_stas_min(i) = interp2(x,y,zmin,lon_sta(i), lat_sta(i));
    depths_at_stas_max(i) = interp2(x,y,zmax,lon_sta(i), lat_sta(i));
end

% vecs = [21 19 13]
% depths_at_stas_mean(vecs)
% depths_at_stas_min(vecs)
% depths_at_stas_max(vecs)
%%
depths_at_stas = depths_at_stas_mean;
[lowest_val,index_low] = min(depths_at_stas)
[highest_val,index_high] = max(depths_at_stas)
sta(index_low)
sta(index_high)
scatter3(lon_sta, lat_sta, depths_at_stas, 'k', 'filled')
text(lon_sta-0.01, lat_sta, depths_at_stas, sta)
hold off

directory = sprintf('/home/a/akfarrell/Uturuncu/Synth_data');
filename = sprintf('elevs_APMB_stas.fig');
filename_wPath = fullfile(directory,filename);
savefig(r, filename_wPath);