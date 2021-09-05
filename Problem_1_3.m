clc;clear;close all;
C=0.85;
A=0.7^2*pi;
V=500*5^2*pi;
a=-9.957*10^(-5);
b=0.00246;
c=0.0007509;
d=-0.127922301898197;
tspan=0:0.02:15*1000; % ����ʱ����t_left:h:t_right
p0=100; % ��ʼѹǿ
ph=160; % ��ѹ��ѹ��
rho=@(p)(exp(a/b*exp(b*p)-c/b*exp(-b*p)-d));
rho_h=rho(ph);

for t0 = 0.030:0.001:0.071 % ����ʱ��/ms
    p=@(t,p)(rho_h*C*A.*sqrt(2*(160-p)/rho_h).*uQin(t,t0)-Qout(t).*(exp(a/b*exp(b*p)-c/b*exp(-b*p)-d))./(V.*(exp(a/b*exp(b*p)-c/b*exp(-b*p)-d))*(a*exp(b*p)+c*exp(-b*p)))); % dp/dt=f(p,t)
    Sol=RK(p,tspan,p0); % �Ľ������������΢�ַ���
    x=tspan; % ʱ��t
    y=Sol; % ѹǿp
    s = 0;
    h = 0;
    jf = 0;
    for i = 745000:750000  %ѹǿ���Դﵽ�����ֵ
        s = s + y(1,i);
    end
    maxi = s/5000;
    
    for j = 1:6600  %�ﵽ���ֵ����Ҫ��֮��
        for k = 1:50
            h = h + y(1,(j-1)*50+k);
            hf = h/5000;
        end
        if hf >= maxi
            jf = j;
        break
        else
            continue
        end
    end   
    
    figure(1)
    plot(t0,maxi,'k.','markersize',8)% �ȶ�ѹǿ
    hold on
    figure(2)
    plot(t0,jf,'k.','markersize',8)% ѹǿ�ﵽ�ȶ������ʱ�� 
    hold on
end
figure(1)
xlabel('����ÿ�ο�����ʱ��t_0/ms','interpreter','tex')
ylabel('������ѹ���ս����ȶ�ֵP_{�ȶ�}/MPa')
grid on
figure(2)
xlabel('����ÿ�ο���ʱ��t_0/ms','interpreter','tex')
ylabel('������ѹ�ﵽ�����ȶ�����ʱ��T_{�ȶ�}/ms')
grid on

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