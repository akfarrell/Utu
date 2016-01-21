function y = funcshun(b)
t = (0:0.1:5);
for i = 1:numel(t)
    y(i) = t*exp(-b*t);
end
plot(t, y)