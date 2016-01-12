function all_eq_all_filt_struct = text_read()
% text_read  Read in text files for each eq, all filters
%   all_eq_all_filt_struct = text_read()  creates output struct of all eq's
%   all filters
%
% all_eq_all_filt_struct = text_read()
%
%   Inputs:
%       none
%   Output:
%       all_eq_all_filt_struct - struct with fields:
%           struct.EQ_name.TextFileName.Station
%           struct.EQ_name.TextFileName.Q
%           struct.EQ_name.TextFileName.Frequency
%
%       ------- ex: -------
%       all_eq_all_filt_struct.ESZ1.ESZ1_output_diffFreq_0pt1875_0pt3750.Q
%
%   Author: Alexandra Farrell 2015/6/16
%
clear
%directories = dir('*SZ*');
directories = dir('JSZ*');
for k = length(directories):-1:1
    if ~directories(k).isdir
        directories(k) = [];
        continue
    end
end

for t = 1:length(directories)
    dirName1 = [directories(t).name];
    dirName = fullfile(dirName1,'text')
    cd(dirName)
    files = dir('*.txt');
    
    for i=1:length(files)
        fileName = [files(i).name];
        corrFileName = strrep(fileName, '.txt','');
        corrFileName = strrep(corrFileName, '.','pt');
        fid = fopen(files(i).name);
        C = textscan(fid, '%s %10.5f %10.3f %10.4f');
        fclose(fid);
        f = {'Station','Q','Freq', 'Time_offset'};
        structure.(dirName1).(corrFileName) = cell2struct(C,f,2);
    end
    cd ../..
end
all_eq_all_filt_struct = structure;