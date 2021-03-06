b=[1, -2, 1];
a=[1, -1.95998, 0.9615];
d=[2, -2, 0];
c=[1, -0.96148, 0];
g=[0.0013094, 367.1433];
% g =[1, 1];
% g=[0.980369946384960777585604319028789177537 / 2, 0.980740725797664247842533313814783468843/2];

b0=g(1) * g(2) * conv(b, d);
a0=conv(a, c);

sos1=[g(1) * b,a];
sos2=[g(2) * d,c];
sos=[sos1;sos2];
[b3, a3]=sos2tf(sos);

W=linspace(-pi, pi, 256);
[h3,w3]=freqz(b3, a3, W);
[h2,w2]=freqz(sos,W);
[h1,w1]=freqz(b0, a0, W);
% h2
plot(w1/pi,20*log10(abs(h1)));
hold on;
plot(w2/pi,20*log10(abs(h2)), '*');
hold on;
plot(w3/pi, 20*log10(abs(h3)), 'g');
