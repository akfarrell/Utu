function [m,I] = corr_arrival_picks(data, sta, name, fil)
if strcmp(name, 'KTSZ1')
    [m,I] = KTSZ1_corr_arrival_picks(data, sta, fil);
else
    [m,I] = nanmax(abs(data));
end