clc;clear;close all;
C=0.85;
A=0.7^2*pi;
V=500*5^2*pi;
a1=-0.000000656979419;
a2=0.000522920735556;
a3=0.804277706982135;
tspan=0:0.02:15*1000; % ����ʱ����t_left:h:t_right
p0=100; % ��ʼѹǿ
ph=160; % ��ѹ��ѹ��
rho=@(p)(a1*p^2+a2*p+a3);
rho_h=rho(ph);
epsilon = 0.25;

for t0=0.745:0.001:0.765
    p=@(t,p)((rho_h*C*A.*sqrt(2*(160-p)/rho_h).*uQin(t,t0)-Qout(t).*(a1*p.^2+a2*p+a3))./(V*(2*a1*p+a2))); % dp/dt=f(p,t)
    Sol=RK(p,tspan,p0); % �Ľ������������΢�ַ���
    x=tspan;
    y=Sol;
%     figure(1)
%     plot(x,y)
%     hold on
    maxi = mean(y(745000:750000));
    for j = 1:150 % ����ﵽ���ֵ����Ҫ��ʱ��
        hf = mean(y(1,(j-1) * 5000 + 1:(j-1) * 5000 + 5000));
        jf = 0;
        if hf >= (150 - epsilon)
            jf = ((j-1)*5000)/50;
            break
        end
    end
    figure(2)
    plot(t0,maxi,'k.','markersize',8)
    xlabel('����ÿ�ο���ʱ��t_0 / ms','interpreter','tex','fontsize',18)
    ylabel('������ѹ�ȶ�ֵP_{�ȶ�} / MPa','interpreter','tex','fontsize',18)
    hold on
    grid on
    xlim([0.744,0.766])
    figure(3)
    plot(t0,jf,'k.')
    xlabel('����ÿ�ο���ʱ��t_0 / ms','interpreter','tex','fontsize',18)
    ylabel('������ѹ�ﵽ�ȶ�ֵ����ʱ��T_{�ȶ�} / ms','interpreter','tex','fontsize',14)
    hold on
    grid on
    xlim([0.744,0.766])
end

function Sol=RK(f,xspan,y0) % �Ľ��������
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

function y = Qout(t) % ����������ʱ��Ĺ�ϵ
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

function y=uQin(t,t0) % ��λ��Ծ����
y=t;
for i = 1:length(t)
    if mod(t(i),t0+10)<=t0
        y(i)=1;
    else
        y(i)=0;
    end
end
end