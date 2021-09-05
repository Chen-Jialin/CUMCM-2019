clc;clear;close all;
load r_theta.mat;
fo = fitoptions('Method','nonlinearLeastSquares',...
               'Lower',[-Inf,-Inf,-Inf],...
               'Upper',[Inf,Inf,Inf],...
               'StartPoint',[1 1 1]);
ft = fittype('a*cos(b*x)+c','options',fo);
[curve1,gof1] = fit(theta,r,ft)
plot(theta,r,'k.','markersize',8)
hold on
grid on
plot(curve1,'k')
legend('数据点','拟合线','fontsize',14)
xlim([0,2 * pi])
xlabel('极角\theta / rad','interpreter','tex','fontsize',18)
ylabel('极径r / mm','interpreter','tex','fontsize',18)
