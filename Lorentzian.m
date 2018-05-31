function [ya] = Lorentzian()
ak = 5;
fk = 654;
N = 2249; 
f = 250:N;
dk = 10; 
y1 = ak * 1./(1.+((fk-f)./(dk/2)).^2);
ak = 10;
fk = 1002;
N = 2249 
f = 250:N;
dk = 20;
y2 = ak * 1./(1.+((fk-f)./(dk/2)).^2);
ak = 15;
fk = 1455;
N = 2249; 
f = 250:N;
dk =20;
y3 =  ak * 1./(1.+((fk-f)./(dk/2)).^2);
x = y1+y2+y3;
h = awgn(x,5,10)
y=2500./f+4
b=250:1:2249
ya=x
%plot(ya);
%xlswrite('C:\Users\fangzheng\Desktop\Lorentzian1.xlsx',x.')
end