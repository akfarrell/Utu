function [index_values, time_values, m_values] = amp_calc(w, index_values, range_values)
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
% Modified 4/2015 and 10/2015 by Alexandra Farrell

    %w = waveform_nonempty(w); % get rid of empty waveform objects
    if numel(w)==0
        warning('no waveforms to plot')
        return
    end
    
    w_raw = waveform(ds, scnl, index)
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

    for wavnum = 1:nwaveforms
        data=get(w(wavnum),'data');
        freq = get(w(wavnum), 'freq');
        dnum(1) = datenum(get(w(wavnum),'start'));
            for l = 2:numel(data)
                dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
            end
        sta=get(w(wavnum),'station');
        chan=get(w(wavnum),'channel');

       
 %%       
        data = data(index-range_val:index+range_val);
        
        numel(data);
        if %add case where you'd be looking for maximum, or comment out the loop and just put nanmax or nanmin
            [m, I] = nanmax(data);
        else
            [m, I] = nanmin(data);
        end
        time_value = dnum(I+(index-range_val));

        index_values(wavnum) = index+I;
        time_values(wavnum) = time_value;
        m_values(wavnum) = m;
        clear(w)
    end
