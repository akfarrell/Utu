%% Q4
x = linspace(0,5,100);
for i = 1:numel(x)
    fx(i) = cos(2*pi*x(i));
end
plot(x,fx)

%4.1
legend('fx')
title('Function Plot')
grid()

%4.2
hold on
for i = 1:numel(x)
    gx(i) = x(i)*exp(-x(i));
end
plot(x,gx, 'r')

%4.3
figure
subplot(2,1,1), plot(x,fx)
subplot(2,1,2), plot(x,gx, 'r')

%% Q5
close all
z1 = 3+4i;
z2 = 5-4i;

%5.1
real(z1)
imag(z1)

%5.2
abs(z2)
angle(z2)

%5.3
z1+z2

%5.4
plot(z1, 'o')
xlabel('Real')
ylabel('Imaginary')

%% Q6
%doing this in the same scrip to keep things consolidated

numbs = (1:1:10);
for i = 1:numel(numbs)
    squares(i) = numbs(i)^2
end
sum(squares)
%or:
sum(numbs.^2)

%% Q7
%in funcshun.m
b = 2
y = funcshun(b)

%% Q8
syms t
xt = cos(t);
diff(xt)
diff(xt,2)

%% Q9
int(t^3, t, -1,1)
int(exp(-t^2), t)
