t = linspace(0,2*pi,1000);
for i = 1:numel(t)
    x(i) = 16*sin(t(i))^3;
    y(i) = 13* cos(t(i)) - 5*cos(2*t(i)) - 2*cos(3*t(i)) - cos(4*t(i));
end

close all;
h = figure;
plot(x,y, 'Color', 'r', 'LineWidth', 2)
fill(x,y,'r')
hgexport(h, 'heart.png', hgexport('factorystyle'), 'Format', 'png');