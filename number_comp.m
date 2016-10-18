%Creates cell arrays with visualization of data from text files. Bam.
clc
siteStruct = loadSiteTable('/raid/data/antelope/databases/PLUTONS/dbmerged'); %load in site table to compare
siteSta = siteStruct.sta;
siteSta = siteSta(9:numel(siteSta)); %remove Lazufre stations
siteSta = sort(siteSta);

cellyQ = cell(1+numel(siteSta), 5); %make cell array
cellyQ(2:numel(siteSta)+1,1) = siteSta;

t = text_read();
mynames = fieldnames(t);
cellyQ(1,2:numel(mynames)+1) = mynames;
cellyT = cellyQ;

%for number, need to use celly{index, index} = number;
for k = 1:numel(siteSta)
    for j = 1:numel(mynames)
        %mynames{j}
        mynames2 = fieldnames(t.(mynames{j}));
        filename = sprintf('%s_output_diffFreq_0pt3750_1pt5000',mynames{j});
        index = find(strcmp(mynames2, filename));
        Q = t.(mynames{j}).(mynames2{index}).Q;
        Stationzz = t.(mynames{j}).(mynames2{index}).Station;
        Time_offset = t.(mynames{j}).(mynames2{index}).Time_offset;
        len = numel(Stationzz);
        for i=1:len
            if strcmp(Stationzz{i}, siteSta{k})
                cellyQ{k+1,j+1} = Q(i);
                cellyT{k+1,j+1} = Time_offset(i);
            end
        end
    end
end
cellyQ
cellyT