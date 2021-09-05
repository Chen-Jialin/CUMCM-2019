C=0.85;
A=0.7^2*pi;
V=500*5^2*pi;
S=2.5^2*pi;
a1=-0.000000656979419;
a2=0.000522920735556;
a3=0.804277706982135;
a=[a1,a2,a3];
b1=2.413;
b3=4.826;
tspan=0:1:10*1000; 
h=tspan(2)-tspan(1);
TT=max(tspan)-min(tspan);
ph=0.5;
rho_h=(a1*0.5^2+a2*0.5+a3);
p0_2=100; 
Omiga=0.01:0.002:0.05;
Error=Omiga;
for jj = 1:length(Omiga)
omiga=Omiga(jj);
To=2*pi/omiga;
tt=0:0.15:min(To,max(tspan));
yy=tt;
Q=0;
for ii = 1:length(tt)
t=tt(ii);
Q=Q+rho_h*C*A*sqrt(2*uPos(ph-100)/ph)*h/2;
rho_h = ((a1*0.5^2+a2*0.5+a3)*(20+7.239*2.5^2*pi-2.413*2.5^2*pi)/(20+7.239*2.5^2*pi-2.5^2*pi*(b1*sin(omiga*t-pi/2)+b3)))-Q/(20+7.239*2.5^2*pi-2.5^2*pi*(b1*sin(omiga*t-pi/2)+b3));
ph = -(a2-sqrt(a2^2-4*a1*a3+4*a1*rho_h))/(2*a1);
Q=Q+rho_h*C*A*sqrt(2*uPos(ph-100)/ph)*h/2;
if ph < 0.5
    ph=0.5;
end
yy(ii)=ph;
end

p=@(t,p)(((a1*p1(t,To,tt,yy).^2+a2*p1(t,To,tt,yy)+a3)*C*A*sqrt(2*uPos(p1(t,To,tt,yy)-p)/(a1*p1(t,To,tt,yy).^2+a2*p1(t,To,tt,yy)+a3))-C.*AQout(t).*sqrt(2*uPos(p-0.1)./(a1*p.^2+a2*p+a3)).*(a1*p.^2+a2*p+a3))./(V*(2*a1*p+a2)));
Sol=RK(p,tspan,p0_2); 
t=tspan;
p=Sol;
Error(jj)=sum((p-100).^2)*h/TT;
end
plot(Omiga,Error)

function y=p1(t,To,tt,yy)
    tj=mod(t,To);
    [~, ll]=min(abs(tt(:)-tj));
    y=yy(ll);
end

function y=uPos(p)
    if p < 0
        y=0;
    else
        y=p;
    end
end

function y = AQout(t) 
if mod(t,100)<0.45
    t=mod(t,100);
    H=(-765.849117266945*t^5+562.147407600656*t^4-92.8313670316711*t^3+8.0818121167249*t^2-0.227266922978743*t);
else
    if mod(t,100)<=2
        H=2;
    else
        if mod(t,100)<=2.45
            t=mod(t,100);
            H=(-767.010424773208*(2.45-t)^5+565.099154813215*(2.45-t)^4-94.2388551044277*(2.45-t)^3+8.23892966960556*(2.45-t)^2-0.233733910408773*(2.45-t));
        else
            H=0;
        end
    end
end
Hm=(sqrt(1.25^2+0.7^2)-1.25)/tan(9/108*pi);
if H<Hm
    y=pi*(tan(9/180*pi))^2*H^2+2*pi*tan(9/180*pi)*H*1.25;
else
    y=pi*0.7^2;
end
end

function Sol=RK(f,xspan,y0) 
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