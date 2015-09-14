function all_eq_all_filt_struct = text_read()
%read in each eq, all filters
clear
directories = dir('*SZ*');
for k = length(directories):-1:1
    if ~directories(k).isdir
        directories(k) = [];
        continue
    end
end

for t = 1:length(directories)
    dirName = [directories(t).name];
    cd(dirName)
    files = dir('*.txt');
    
    for i=1:length(files)
        fileName = [files(i).name];
        corrFileName = strrep(fileName, '.txt','');
        corrFileName = strrep(corrFileName, '.','pt');
        fid = fopen(files(i).name);
        C = textscan(fid, '%s %10.5f %10.3f');
        fclose(fid);
        f = {'Station','Q','Freq'};
        struct.(dirName).(corrFileName) = cell2struct(C,f,2);
    end
    cd ..
end

all_eq_all_filt_struct = struct;