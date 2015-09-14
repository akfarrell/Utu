% Frequency vs Period visualizations
x = linspace(0.1,3,1000);
y_series = [];
for i=1:numel(x)
    y_series(i) = 1/x(i);
end

figure
plot(x,x,'k',x,y_series,'r');
title('Frequency (k) vs T (r)')

diff_series = [];
for i=2:numel(y_series)
    diff_series(i-1) = y_series(i-1)-y_series(i);
end

diff_series(1000) = 0;
figure
plot(x,diff_series);
title('Difference between values of T')

s = [x;y_series]'

% dif=fopen('output_FvsT.txt','w');
% for i=1:length(x)
% fprintf(dif,'%10.5f %10.5f\n',x(i),y_series(i));
% end