function [m,I] = KTSZ1_corr_arrival_picks(data, sta, fil)
    if fil(1)==0.375 && fil(2)==1.5
        if strcmp(sta, 'PLRV') || strcmp(sta, 'PLSS')
            data = data(1:500);
            [m,I] = nanmin(data);
        else
            [m,I] = nanmin(data);
        end
    elseif fil(1)==0.375 && fil(2)==0.75
        if strcmp(sta, 'PLRV') || strcmp(sta, 'PLSS') || strcmp(sta, 'PLBR') 
            data = data(1:500);
            [m,I] = nanmin(data);
        else
            [m,I] = nanmin(data);
        end
    elseif fil(1)==0.75 && fil(2)==1.5
        if strcmp(sta, 'PLMD') || strcmp(sta, 'PLLL') || strcmp(sta, 'PLSS') || strcmp(sta, 'PLAR')
            data = data(1:500);
            [m,I] = nanmin(data);
        elseif strcmp(sta, 'PLCM') || strcmp(sta, 'PLRV')
            data = data(1:400);
            [m,I] = nanmin(data);
        elseif strcmp(sta, 'PLSQ')
            data = data(1:300);
            [m,I] = nanmin(data);
        else
            [m,I] = nanmin(data);
        end
    elseif (fil(1)==0.1875 && fil(2)==0.75) || (fil(1)==0.375 && fil(2)==1.5)
        if strcmp(sta, 'PLRV')
            data = data(1:450);
            [m,I] = nanmin(data);
        else
            [m,I] = nanmin(data);
        end
    elseif fil(1)==0.75 && fil(2)==3
        if strcmp(sta, 'PLMD') || strcmp(sta, 'PLLL') || strcmp(sta, 'PLQU')
            data = data(1:500);
            [m,I] = nanmin(data);
        elseif strcmp(sta, 'PLSQ') || strcmp(sta, 'PLRV')
            data = data(1:350);
            [m,I] = nanmin(data);
        else
            [m,I] = nanmin(data);
        end
    elseif fil(1)==0.1875 && fil(2)==3
        if strcmp(sta, 'PLRV') || strcmp(sta, 'PLQU')
            data = data(1:500);
            [m,I] = nanmin(data);
        else
            [m,I] = nanmin(data);
        end
    else
        [m,I] = nanmin(data);
    end
