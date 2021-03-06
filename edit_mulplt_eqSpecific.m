function [index_values, time_values, m_values] = edit_mulplt_eqSpecific(w, alignWaveforms, absMax, absMin, name, fil,tshift, id, delay2)
%added last input for eq_noSwaves to tell difference between split-up
%waveform objects
%MULPLT Plot multiple waveform objects in a figure. is inspired by the 
%Seisan program of the same name
%   mulplt(w, alignWaveforms) 
%   where:
%       w = a vector of waveform objects
%       alignWaveforms is either true or false (default)
%   mulplt(w) will plot a record section, i.e. each waveform is plotted
%   against absolute time.
%   mulplt(w, true) will align the waveforms on their start times.

% Glenn Thompson 2014/11/05, generalized after a function I wrote in 2000
% to operate on Seisan files only
% Modified 4/2015 by Alexandra Farrell

    %w = waveform_nonempty(w); % get rid of empty waveform objects
    if numel(w)==0
        warning('no waveforms to plot')
        return
    end
    
    if ~exist('alignWaveforms', 'var')
            alignWaveforms = false;
    end
    
    if ~exist('id', 'var')
        id = 'eq';
    end
    index_values = [];
    time_values = [];
    m_values = [];
    % get the first start time and last end time
    [starttimes endtimes]=gettimerange(w);
    snum = nanmin(starttimes);
    enum = nanmax(endtimes);
    
    
    % get the longest duration - in mode=='align' 
    durations = endtimes - starttimes;
    maxduration = nanmax(durations); 
    SECSPERDAY = 60 * 60 * 24;
    
    nwaveforms = numel(w);
    
	h=figure;
	trace_height=0.95/nwaveforms;
	left=0.125;
	width=0.8;
    %set(h, 'Position', [1000 1000 width+525 trace_height*nwaveforms+1000])
    set(h, 'Position', [1000 1000 width+1000 trace_height*nwaveforms+1250])
    for wavnum = 1:nwaveforms
        data=get(w(wavnum),'data');
        freq = get(w(wavnum), 'freq');
        dnum(1) = datenum(get(w(wavnum),'start'));
            for l = 2:numel(data)
                dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
            end
        sta=get(w(wavnum),'station');
        chan=get(w(wavnum),'channel');
        ax(wavnum)=axes('Position',[left+0.025 0.95-wavnum*trace_height+0.05 width trace_height*.75]); 
        
        if alignWaveforms
            plot((dnum-min(dnum))*SECSPERDAY, data,'-k');
            set(gca, 'XLim', [0 maxduration*SECSPERDAY]);
        else
            if strcmp(id, 'selected')
                if strcmp(chan,'HHZ') == 1 || strcmp(chan, 'HHE') == 1
                    plot(dnum, data,'-k');
                end
            else
              plot(dnum, data,'-k');  
            end
            set(gca, 'XLim', [snum enum]);
        end
        ylim([absMin absMax])
        if strcmp(id,'selected')
            if strcmp(chan,'HHZ') == 1 || strcmp(chan, 'HHE') == 1
                if strcmp(sta,'PL07')==1 || strcmp(sta,'PLSE') == 1 || strcmp(sta,'PLMK') ==1
                   ylabel(sprintf('%s\n%s ',sta,chan),'FontSize',9,'Rotation',90, 'FontWeight', 'bold');
                else
                   ylabel(sprintf('%s\n%s ',sta,chan),'FontSize',9,'Rotation',90); 
                end
            end
        else
            ylabel(sprintf('%s\n%s ',sta,chan),'FontSize',9,'Rotation',90);
        end
       
 %%       
        %set(gca,'YTick',[],'YTickLabel',['']);
        datetick('x', 'keeplimits');
        if wavnum<nwaveforms
           set(gca,'XTickLabel',['']);
        end
        time_value = tshift(wavnum);
        wavnum;
        tolerance = 0.0000001;
        for i=1:numel(dnum)
            t = find(dnum(i)>(time_value-tolerance) & dnum(i)<(time_value+tolerance));
            if t~=0
                index = i;
            end
        end
        
        if fil(2)/fil(1)==2
            range_val = ceil(15/((fil(2)-fil(1))/1.2));
        elseif fil(2)/fil(1)==4
            range_val = ceil(15/((fil(2)-fil(1))/4));
        elseif fil(2)/fil(1)==16
            range_val = 50;
        end

        if strcmp(name, 'KTSZ4')
            if fil(1) == 0.375
                range_val = range_val+50;
%             elseif fil(1) == 0.1875 && fil(2) == 3.000
%                 range_val = range_val + 100
            %elseif fil(1) == 0.1875
             %   range_val = range_val+2;
            end
        end
        index+range_val;
        %index-range_val;

        if strcmp(name, 'JSZ3') && fil(1) == 0.1875 && fil(2) == 0.7500 && (wavnum == 18 || wavnum == 22) %exception for this filter PLLL & PLWB
            range_val = range_val-50;
        elseif strcmp(name, 'SSSZ1') && fil(1) == 0.1875 && fil(2) == 3.000
            range_val = 10;
        end
        index
        range_val
        if strcmp(name, 'JSZ1') || strcmp(name, 'JSZ2') || strcmp(name, 'JSZ3') || strcmp(name, 'JSZ4') || strcmp(name, 'SSSZ1') || strcmp(name, 'SSSZ2') || strcmp(name, 'SSSZ3')
            if strcmp(name, 'JSZ3') && fil(1) == 0.75 && fil(2) == 1.5000 && (wavnum == 10 || wavnum == 15)
                data = data(index-range_val:index+range_val+30);
            elseif strcmp(name, 'SSSZ2') && fil(1) == 0.1875 && fil(2) == 3.000 %&& (wavnum == 10 || wavnum == 15)
                range_val = range_val+15;
                if wavnum == 1
                    data = data(index+20:index+range_val+50);
                else
                    data = data(index-range_val:index+range_val);
                end
            elseif strcmp(name, 'SSSZ3') && (fil(1) == 0.1875 && fil(2)== 3.000) && wavnum~=24
                range_val = range_val + 5;
                data=data(index-range_val:index+range_val);
            elseif strcmp(name, 'SSSZ3') && (fil(1)==0.7500 && fil(2)==1.500) && wavnum~=24
                range_val = range_val + 10;
                data=data(index-range_val:index+range_val);
            elseif strcmp(name, 'SSSZ3') && ((fil(1) == 0.375 && fil(2) == 1.500) || (fil(1) == 0.1875 && fil(2) == 3.000) || (fil(1)==0.7500 && fil(2)==1.500) || (fil(1)==0.1875 && fil(2)==0.7500) || (fil(1)==0.375 && fil(2)==0.7500) || (fil(1)==0.7500 && fil(2)==3.000)) && (wavnum==24 ||(fil(1)==0.1875 && fil(2)==0.75 && wavnum==2))
                range_val = 90;
                if fil(1)==0.1875 && fil(2)==0.75 && wavnum==2
                    data = data(index+10:index+range_val);
                else
                    data = data(index+50:index+range_val);
                end
            elseif strcmp(name, 'SSSZ3') && fil(1)==0.1875 && fil(2)==0.7500 && wavnum == 19
                data = data(index-range_val:index+range_val+10);
            elseif strcmp(name, 'SSSZ3') && fil(1)==0.1875 && fil(2)==0.7500 &&  wavnum==10
                data = data(index-range_val:index);
            elseif strcmp(name, 'SSSZ3') && fil(1)==0.7500 && fil(2)==3.000 && (wavnum == 4 || wavnum == 13 || wavnum == 15 || wavnum == 20)
                range_val = range_val+30;
                data = data(index-range_val: index);
            elseif strcmp(name, 'SSSZ2') && fil(1)==0.3750 && fil(2)==0.7500 && wavnum == 8
                data = data(index-range_val:index);
            elseif strcmp(name, 'SSSZ3') && ((fil(1)==0.3750 && fil(2)==0.7500 && wavnum == 14) || ((fil(1)==0.3750 && fil(2)==1.500) && (wavnum == 14 || wavnum == 5)))
                if fil(1)==0.375 && fil(2)==1.500 && wavnum == 5
                    data = data(index:index+90);
                else
                    data = data(index-range_val:index+90);
                end
            elseif strcmp(name,'KTSZ4')
                if fil(1) == 0.375 && fil(2) == 1.5
                    data = data(index-100:index+50);
                elseif fil(1) == 0.1875 && fil(2) == 3.000
                    data = data(index-20:index+range_val+100);
                end
            else
                data = data(index-range_val:index+range_val);
                numel(data)
                pasta = 'dessert'
            end
        elseif strcmp(name, 'SSSZ5')
            if wavnum == 11
                range_val = -10;
                data = data(index-range_val:index+range_val+70);
            elseif wavnum == 3
                range_val = 70;
                data = data(index-range_val:index+range_val);
            elseif wavnum == 9
                range_val = 20;
                data = data(index-range_val:index+50);
            elseif wavnum == 20
                range_val = 0;
                data = data(index-range_val:index+130);
            elseif wavnum == 18
                range_val = -50;
                data = data(index-range_val:index+150);
            elseif wavnum == 4
                range_val = 100;
                data = data(index-range_val:index);
            else
                hu = 'bu'
                range_val = 30;
                data = data(index-range_val:index+range_val+70);
            end
        else
            data = data(index-range_val:index+range_val);
            turkey = 'dinner'
        end
        
        
        numel(data);
        if strcmp(name, 'KTSZ3') || strcmp(name, 'JSZ4') || strcmp(name, 'SSSZ3') || strcmp(name, 'SSSZ1') || strcmp(name, 'SSSZ4') || strcmp(name, 'SSSZ5')
            [m, I] = nanmax(data);
        else
            [m, I] = nanmin(data);
        end
        
        if strcmp(name, 'KTSZ4')
            if fil(1) == 0.375 && fil(2) == 1.5
                time_value = dnum(I+(index-100)-1);
                I_fullData = I+(index-100)-1;
            elseif fil(1) == 0.1875 && fil(2) == 3.000
                time_value = dnum(I+(index-20));
                I_fullData = I+(index-20)-1;
            else
                time_value = dnum(I+(index-range_val));
                I_fullData = I+(index-ranve_val)-1;
            end
        elseif strcmp(name, 'SSSZ2') && fil(1) == 0.1875 && fil(2) == 3.000 && wavnum==1
            time_value = dnum(I+(index+20));
            I_fullData = I+(index+20)-1;
        elseif strcmp(name, 'SSSZ3') && ((fil(1) == 0.375 && fil(2) == 1.500) || (fil(1) == 0.1875 && fil(2) == 3.000) || (fil(1)==0.7500 && fil(2)==1.500) || (fil(1)==0.375 && fil(2)==0.7500) || (fil(1)==0.1875 && fil(2)==0.7500) || (fil(1)==0.7500 && fil(2)==3.000)) && (wavnum==24 ||(fil(1)==0.1875 && fil(2)==0.75 && wavnum==2))
            if fil(1)==0.1875 && fil(2)==0.75 && wavnum==2
                time_value = dnum(I+(index+10));
                I_fullData = I+(index+10)-1;
            else
                time_value = dnum(I+(index+50));
                I_fullData = I+(index+50)-1;
            end
            
        elseif strcmp(name, 'SSSZ3') && fil(1)==0.375 && fil(2)==1.5 && wavnum==5
            time_value = dnum(I+(index));
            I_fullData = I+(index)-1;
        else
            wavnum
            time_value = dnum(I+(index-range_val));
            I_fullData = I+(index-range_val)-1;
            cat = 'dog'
        end
time_value
%         derp = dnum(index);  
%         time_value1 = dnum(I+(index-range_val)-200);
%         time_value2 = dnum(I+(index-range_val)+200);
%         time_value3 = dnum(I+(index-range_val)-100);
%         time_value4 = dnum(I+(index-range_val)+60);
%         time_value5 = dnum(I+(index-range_val)-150);
%         time_value6 = dnum(I+(index-range_val)+150);
        %index_values(wavnum) = index+I;
        index_values(wavnum) = I_fullData;
        time_values(wavnum) = time_value;
        m_values(wavnum) = m;
        time_value1 = dnum(index_values(wavnum)-110);
        time_value2 = dnum(index_values(wavnum)+90);
        time_value3 = dnum(index_values(wavnum)-150);
        time_value4 = dnum(index_values(wavnum)+60);
        hold on
        yl=ylim;
        %line([get(w(wavnum), 'EX_ARR_TIME'), get(w(wavnum), 'EX_ARR_TIME')], [yl(1), yl(2)], 'Color', 'k');
        %hold on
        line([time_value, time_value], [yl(1), yl(2)], 'Color', 'k', 'LineWidth', 2);
        if exist('delay2', 'var')
            line([time_value-(delay2(wavnum)/SECSPERDAY), time_value-(delay2(wavnum)/SECSPERDAY)], [yl(1), yl(2)], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 4);
        end
        %line([time_value1, time_value1], [yl(1), yl(2)], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 4);
        %line([time_value2, time_value2], [yl(1), yl(2)], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 4);
         %line([time_value3, time_value3], [yl(1), yl(2)], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 4);
         %line([time_value4, time_value4], [yl(1), yl(2)], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 4);
%         line([time_value5, time_value5], [yl(1), yl(2)], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 4);
%         line([time_value6, time_value6], [yl(1), yl(2)], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 4);
        %data_construct = [dnum; data]'
%         if wavnum==1
%            title('','FontSize',10);
%         end
        %axis tight;
        
        % display mean on left, max on right
         a=axis;
%         tpos=a(1)+(a(2)-a(1))*.02; %horizontal position
         dpos=a(3)+(a(4)-a(3))*.7; %vertical position
%         text(tpos,dpos,sprintf('%5.0f',nanmean(data)),'FontSize',10,'Color','k');
% %         tpos=a(1)+(a(2)-a(1))*.4;
% %         text(tpos,dpos,sprintf(' %s',datestr(starttimes(wavnum),30)),'FontSize',10,'Color','k');
         tpos=a(1)+(a(2)-a(1))*.9;
         %text(tpos,dpos,sprintf('%5.0f',nanmax(abs(data))),'FontSize',10,'Color','k');

        
    end
%     if exist('ax','var')
%         linkaxes(ax,'x');
%         %samexaxis();
%         hlink = linkprop(ax,'XLim');
%         if ~alignWaveforms
%             %datetick('x', 'keeplimits');
%         end
%     end
xlabel('Time');
if exist('delay2', 'var')
    filename = sprintf('%s_%s_waveforms_%1.4f_%1.4f_delay.pdf',name,id,fil(1),fil(2));
else
    filename = sprintf('%s_%s_waveforms_%1.4f_%1.4f.pdf',name,id,fil(1),fil(2));
end
directory = sprintf('/home/a/akfarrell/Uturuncu/%s/figures', name);
filename_wPath = fullfile(directory,filename);
fp = fillPage(h, 'margins', [0 0 0 0], 'papersize', [8.5 11]);
print(h, '-dpdf', '-r400',filename_wPath);
%hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'pdf');