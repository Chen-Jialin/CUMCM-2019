clc;clear;close all;
C=0.85;
A=0.7^2*pi;
V=500*5^2*pi;
a1=-0.000000656979419;
a2=0.000522920735556;
a3=0.804277706982135;
tspan=0:0.1:20*1000; % 计算时长：t_left:h:t_right
p0=100; % 初始压强
ph=160; % 高压端压力
rho=@(p)(a1*p^2+a2*p+a3);
rho_h=rho(ph);
p_tar=p0;

t0=0.286:0.0002:0.29;
x=t0;
y=t0;
h=tspan(2)-tspan(1);
T=max(tspan)-min(tspan);
for ii = 1:length(t0)
p=@(t,p)((rho_h*C*A.*sqrt(2*(160-p)/rho_h).*uQin(t,t0(ii))-Qout(t).*(a1*p.^2+a2*p+a3))./(V*(2*a1*p+a2))); % dp/dt=f(p,t)
Sol=RK(p,tspan,p0); % 四阶龙格库塔计算微分方程
y(ii)=sum((Sol-p_tar).^2)*h/T;
end
plot(x,y,'k-','linewidth',3)
grid on
xlabel('单向阀每次开启时长t_0 / ms','interpreter','tex','fontsize',18)
ylabel('目标函数Var(\Delta P) / (MPa)^2','interpreter','tex','fontsize',18)

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