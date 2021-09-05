%调参结果  两秒  t_01 = 0.900  t_02 = 0.752  s = 1.75
%调参结果  五秒  t_01 = 0.600  t_02 = 0.752  s = 1.36 
%调参结果  十秒  t_01 = 0.700  t_02 = 0.752  s = 8.62
clc;clear;close all;
C=0.85;
A=0.7^2*pi;
V=500*5^2*pi;
a1=-0.000000656979419;
a2=0.000522920735556;
a3=0.804277706982135;

s = 1.75;
tspan1 = 0:0.02:s*1000;
tspan2 = s*1000:0.02:15*1000;
p0=100; % 初始压强
ph=160; % 高压端压力
rho=@(p)(a1*p^2+a2*p+a3);
rho_h=rho(ph);

t_01 = 0.900;
p=@(t,p)((rho_h*C*A.*sqrt(2*(160-p)/rho_h).*uQin(t,t_01)-Qout(t).*(a1*p.^2+a2*p+a3))./(V*(2*a1*p+a2))); % dp/dt=f(p,t)
Sol1=RK(p,tspan1,p0); % 四阶龙格库塔计算微分方程
x1=tspan1;
y1=Sol1;

t_02 = 0.752;
p1 = y1(1,length(Sol1));
p=@(t,p)((rho_h*C*A.*sqrt(2*(160-p)/rho_h).*uQin(t,t_02)-Qout(t).*(a1*p.^2+a2*p+a3))./(V*(2*a1*p+a2)));
Sol2 = RK(p,tspan2,p1);
x2 = tspan2;
y2 = Sol2;

figure(1)
plot(x1,y1,'k-','linewidth',1)
hold on
plot(x2,y2,'k-','linewidth',1)
grid on
xlabel('时间t / ms','interpreter','tex','fontsize',18)
ylabel('管内油压P / MPa','interpreter','tex','fontsize',18)

function Sol=RK(f,xspan,y0) % 四阶龙格库塔
n=length(xspan)-1;
x0=xspan(1);
xr=xspan(n);
h=(xr-x0)/(n-1);
xn=x0;
yn=y0;
Sol=xspan;
Sol(1)=y0;
for ii = 1:n
    K1=f(xn,yn);
    K2=f(xn+h/2,yn+h*K1/2);
    K3=f(xn+h/2,yn+h*K2/2);
    K4=f(xn+h,yn+h*K3);
    yn=yn+h/6*(K1+2*K2+2*K3+K4);
    xn=xn+h;
    Sol(ii+1)=yn;
end
end

function Sol=O1(f,xspan,y0) % 显示欧拉
n=length(xspan)-1;
x0=xspan(1);
xr=xspan(n);
h=(xr-x0)/(n-1);
xn=x0;
yn=y0;
Sol=xspan;
Sol(1)=y0;
for ii = 1:n
    yn=yn+h*f(xn,yn);
    xn=xn+h;
    Sol(ii+1)=yn;
end
end

function y = Qout(t) % 流出流量与时间的关系
y=t;
for i = 1:length(t)
if mod(t(i),100)<=0.2
    y(i)=100*mod(t(i),100);
else 
    if mod(t(i),100)<=2.2
        y(i)=20;
    else
        if mod(t(i),100)<=2.4
            y(i)=100*(2.4-mod(t(i),100));
        else
            y(i)=0;
        end
    end
end
end
end

function y=uQin(t,t0) % 单位阶跃函数
y=t;
for i = 1:length(t)
    if mod(t(i),t0+10)<=t0
        y(i)=1;
    else
        y(i)=0;
    end
end
end