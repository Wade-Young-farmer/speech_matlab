close all
N = 6;
x1 = [1,1,1,1,0,0,0,0];
x2 = [0,0,1,1,1,1,0,0];
figure
subplot(2,1,1),stem(0:7,x1)
xlim([-2,8]),ylim([0,2])
xlabel('m'),title('x_1(m)')

subplot(2,1,2),stem(0:7,x2)
xlim([-2,8]),ylim([0,2])
xlabel('m'),title('x_2(m)')

x3 = [1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1];
figure,
subplot(2,1,1),
stem(-8:8,x3)
xlim([-8,8]),ylim([0,2])
xlabel('m'),title('x_2(m)以N = 6为周期延拓((x_2(m)))_6')

subplot(2,1,2),
stem(-8:8,fliplr(x3))
xlim([-8,8]),ylim([0,2])
xlabel('m'),title('x_2(m)以N = 6为周期延拓 反褶((x_2(-m)))_6')


x_m = fliplr(x3);
x4 = [1,0,0,1,1,1,x_m];
figure
for i = 0:N-1
    subplot(7,1,i+1),
    stem(0:6,x4(9-i:15-i))
    xlim([-2,8]),ylim([0,1.5])
    
    ylabelString = '(x_2(%d-m))_6';
    str = sprintf(ylabelString,i);

    ylabel(str,'rotation',90,'FontSize',8);
    %xlabel('m'),title('x_2(m)以N = 6为周期延拓 反褶((x_2(-m)))_6')
end
subplot(7,1,7)
    stem(0:6,x1(1:7))
    xlim([-2,8]),ylim([0,1.5])
    
    ylabelString = '(x_2(%d-m))_6';
    str = sprintf(ylabelString,i);

    ylabel('x_1(m)','rotation',90,'FontSize',8);
    
    
x = cconv(x1,x2,6)
figure,stem(0:5,x)
xlim([-2,8]),ylim([0,5])




