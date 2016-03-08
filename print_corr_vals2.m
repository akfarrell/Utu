directory = sprintf('/home/a/akfarrell/Uturuncu/%s/text', eq(earthquake_number).name);
filename2 = sprintf('%s_corrVals2_%1.4f_%1.4f.txt',eq(earthquake_number).name,fil(1),fil(2));
dif=fopen(fullfile(directory,filename2), 'w');
fprintf(dif, '    ');
for i =1:length(stas)
    fprintf(dif, '     %s ', stas(i,:));
end
fprintf(dif, '\n');
for i=1:length(stas)
    fprintf(dif,'%s' ,stas(i,:));
    fprintf(dif, '%10.3f',corr_vals(i,:));
    fprintf(dif, '\n');
end
fprintf(dif, '\n\n');
fprintf(dif, 'min ');
fprintf(dif, '%10.3f', min(corr_vals));
fprintf(dif, '\nmax ');
sorted_vals = sort(corr_vals, 'descend');
fprintf(dif, '%10.3f', sorted_vals(2,:));
fprintf(dif, '\nmean');
fprintf(dif, '%10.3f', mean(corr_vals));
st = fclose('all');
