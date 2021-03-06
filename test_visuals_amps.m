%Creates cell arrays with visualization of data from text files. Bam.
clc
siteStruct = loadSiteTable('/raid/data/antelope/databases/PLUTONS/dbmerged'); %load in site table to compare
siteSta = siteStruct.sta;
siteSta = siteSta(9:numel(siteSta)); %remove Lazufre stations
siteSta = sort(siteSta);
%%
%eqs = {'SSSZ1', 'SSSZ2', 'SSSZ3', 'KTSZ1', 'KTSZ2', 'KTSZ3', 'JSZ1',
%'JSZ3', 'JSZ4'}; % for norm 9
%eqs = {'KTSZ2', 'JSZ4', 'SSSZ1'}; % for norm 3
eqs = {'SSSZ1', 'SSSZ2', 'SSSZ3', 'KTSZ1', 'KTSZ2', 'KTSZ3', 'JSZ1', 'JSZ3', 'JSZ4'};
for i = 1:numel(eqs)
    dirName = fullfile(eqs{i},'text')
    cd(dirName)
    files = dir('*_output_diffFreq_0.3750_1.5000.txt');
    fileName = files.name
    fid = fopen(fileName)
    C = textscan(fid, '%s %10.5f %10.3f %10.4f %10.4f %10.4f');
    fclose(fid);
    f = {'Station','Q','Freq', 'Time_offset', 'Normalized_ems', 'Reverse_norm_ems'};
    structure.(eqs{i}) = cell2struct(C,f,2);
    cd ../..
end
t = structure;
%%

cellyN = cell(1+numel(siteSta), 4); %make cell array
cellyN(2:numel(siteSta)+1,1) = siteSta;
mynames = fieldnames(t);
cellyN(1,2:numel(mynames)+1) = mynames;
cellyR = cellyN;
%%
scale_factor = 150;
%for number, need to use celly{index, index} = number;
for k = 1:numel(siteSta)
    for j = 1:numel(mynames)
        mynames{j};
        Normalized_ems = t.(mynames{j}).Normalized_ems;
        Stationzz = t.(mynames{j}).Station;
        Reverse_norm_ems = t.(mynames{j}).Reverse_norm_ems;
        Reverse_norm_ems = (Reverse_norm_ems+0.01);
        Reverse_norm_ems = Reverse_norm_ems/max(Reverse_norm_ems);
        len = numel(Stationzz);
        for i=1:len
            if strcmp(Stationzz{i}, siteSta{k})
                cellyN{k+1,j+1} = Normalized_ems(i)*scale_factor;
                cellyR{k+1,j+1} = Reverse_norm_ems(i)*scale_factor;
                matN(k,j) = Normalized_ems(i)*scale_factor;
                matR(k,j) = Reverse_norm_ems(i)*scale_factor;
            end
            if j == 2
                
            end
        end
    end
end


cellyN
cellyR
matN
matR
%%
% ---------- Method for making average ---------%

derpy = size(matR);
numRows = derpy(1);
for i=1:numRows
    c(i)=nnz(matR(i,:));
end

sumR = sum(matR,2);
sumN = sum(matN,2);
meanR = sumR'./c
meanN = sumN'./c

%%
% ------------ Plotting ------------%
num_eqs = derpy(2);
%ID - 'norm' or 'revnorm'
%plotting_allEq_amps(siteStruct, siteSta, meanN, num_eqs, 'norm', c,scale_factor)
%plotting_allEq_amps(siteStruct, siteSta, meanR, num_eqs, 'revnorm', c)
