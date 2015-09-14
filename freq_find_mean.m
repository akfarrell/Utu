close all; clc
sig_len = (eq(earthquake_number).enum - eq(earthquake_number).snum)*SECS2DAY;
nwaveforms = numel(w_clean_sort);

for wavnum = 1:nwaveforms
        data = get(w_clean_sort(wavnum),'data');
        stn = get(w_clean_sort(wavnum),'station');
        %[m,I] = corr_arrival_picks(data, stn, eq(earthquake_number).name, fil);
        m = m_values(wavnum);
        dnum=datenum(w_clean_sort(wavnum));
        %time_value = dnum(I)
        data = data(index_values(wavnum)-150:index_values(wavnum)+150);
        %[maxValue,indexMax] = max(abs(fft(data)))
        
        r=abs(fft(data));
        maxValue = mean(r(1:10));
        %tolerance = 1647; %for +-100
        tolerance = 998;
        for i=1:numel(data)
            t = find(data(i)>(maxValue-tolerance) & data(i)<(maxValue+tolerance));
            if t~=0
                indexMax = i
            end
        end

        r=abs(fft(data));
        frequency = indexMax * 100 / numel(data);
%         frequency_center = indexMax * 100 / numel(data);
%         frequency = int(abs(fft(data)))
        figure
        x = linspace(0,numel(data), numel(data));
        plot(x, abs(fft(data)))
        yl(1) = min(abs(fft(data)));
        yl(2) = max(abs(fft(data)));
        hold on
        line([frequency, frequency], [yl(1), yl(2)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
        xlim([0 10])
        xlabel('Frequency (Hz)')
        ylabel('Strength')
        title(sprintf('%s',stn))
        w_clean_sort(wavnum) = addfield(w_clean_sort(wavnum), 'AMP_ABS', abs(m));
        w_clean_sort(wavnum) = addfield(w_clean_sort(wavnum), 'sig_freq', frequency);
end