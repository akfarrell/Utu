for wavnum = 1:nwaveforms
        data = get(w_clean_sort(wavnum),'data');
        stn = get(w_clean_sort(wavnum),'station');
        %[m,I] = corr_arrival_picks(data, stn, eq(earthquake_number).name, fil);
        m = m_values(wavnum);
        dnum=datenum(w_clean_sort(wavnum));
        %time_value = dnum(I)
        data = data(index_values(wavnum)-100:index_values(wavnum)+100);
        [maxValue,indexMax] = max(abs(fft(data)));
        x = linspace(0,100,numel(data));
        frequency = x(indexMax)
        %frequency = indexMax * 100 / numel(data);
        figure
        r = abs(fft(data))/numel(data)
        plot(x, r)
        %plot(x, abs(fft(data))/numel(data))
        yl(1) = min(abs(fft(data))/numel(data));
        yl(2) = max(abs(fft(data))/numel(data));
        hold on
        line([frequency, frequency], [yl(1), yl(2)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
        xlim([0 10])
        xlabel('Frequency (Hz)')
        ylabel('Strength')
        title(sprintf('%s',stn))
        w_clean_sort(wavnum) = addfield(w_clean_sort(wavnum), 'AMP_ABS', abs(m));
        w_clean_sort(wavnum) = addfield(w_clean_sort(wavnum), 'sig_freq', frequency);
end