```
#结构动力学第三次作业
##简谐荷载下结构动力响应

m=545;                      %mass (kg)
K=5.9511111*10^6;           %（kn/m)
kesi=0.01;                  %damping ratio
theta=10*pi;                %the frequency of loads (rad/s)
w_n=(K/m)^0.5;          % natural frequency of vibration
w_d=w_n*(1-kesi^2)^0.5;     %natural frequency of damping system
P_0=267;                    %amplitude of loads (N)
% P_t=P_0.*sin(theta*t);      %function of loads
y_st=(P_0/K);               %displacement (mm)
y_0=0;
v_0=0;

clear
clc
[t_1,y_th]=thearitical_methodfor37(10,0.005);%theoretical solution
[t,y_du]=Duhamelintegralfor37(10,0.005);% solution in time domain
[t_3,y_fre]=fre_methodfor37(10,0.005);%fft and ifft
[t_4,y_digi]=digi_methodfor37(10,0.005);%numerical solution
[t_b,y_t_b] = considering_bitemfor37(10,0.005);%considering the components of vibration at the natural frequency induced by excitation force
%以上项考虑了由激励力引起的按自振频率振动的分量，称为伴随自由振动，随着时间很快衰减。
figure(1)
plot(t_1,y_th,'-b*');
 hold on
plot(t,y_du,':go');
hold on
plot(t_3,real(y_fre),'-.rp')
hold on
plot(t_4,y_digi,'--ks')
hold on
plot(t_b,y_t_b,'k>')
legend('closed-solution','Duhamel-interal','Frequency-domain-method','numerical-method','considering_b');
hold off
%%remark
%两种方法考虑了伴随自由振动，杜哈梅积分的方法和数值方法

function [t,y_du] = Duhamelintegralfor37(aa,bb)
%inputting paremeters.
% aa=10;up bound of time
%bb=0.02;% time interval
%%%%%
t=bb:bb:aa;
tau=t;% time series
%%%%%

m=545;                      %mass (kg)
K=5.9511111*10^6;           %（kn/m)
kesi=0.01;                  %damping ratio
theta=10*pi;                %the frequency of loads (rad/s)
w_n=(K/m)^0.5;              % natural frequency of vibration
w_d=w_n*(1-kesi^2)^0.5;     %natural frequency of damping system
P_0=267;                    %amplitude of loads (N)
P=P_0.*sin(theta*t);        %function of loads
y_st=(P_0/K);               %displacement (mm)
y_0=0;
v_0=0;


%%%%%%%%the main program
for i=1:(aa/bb)
    for j=1:i
        %paramet1(j)=p(j)/(m*w)*bb*sin(w*(t(i)-t1(j)));                               %Duhameil's integral function for non-damping system
        paramet2(j)=P(j)/(m*w_d)*bb*exp(-kesi*w_n*(t(i)-tau(j)))*sin(w_d*(t(i)-tau(j)));%%Duhameil's integral function for damping system
    end
    y_du(i)=sum(paramet2);%%displacement   
end
end


function [t_3,y_fre] = fre_methodfor37(t_p,diata_t )
%t_p=10;
%diata_t=0.01;
nn=0:1:t_p/diata_t-1;
t_3=nn*diata_t;
N=t_p/diata_t;

theta=10*pi;%the unit is radian/second load
p_0=267;
p=p_0*sin(theta*t_3);
p_f=fft(p,N);%fft变换后包含了正负频率的完整数据，正负频率分别存储在结果的前后半部分，matlab在振动信号中的应用一书中有讲。
p_ff=[p_f(N/2+1:end),p_f(1:N/2)];%after fft transform，改成负、正频率对应的幅值
f=(-N/2:1:N/2-1)/(N*diata_t);%频率是从-N/2 到正N/2
%plot(f(1:N/2),abs(y(1:N/2))*2/N)
k=5951.11111*10^3;%the standard unit
w_n=104.5;%the unit is radian/second natrual vibration
kesai=0.01;
beita=2*f*pi/w_n;
H_iw=(1/k)*((1-beita.^2)+i*2*kesai*beita).^-1;%频响函数也是负正频率对应的值
y_f=p_ff.*H_iw;%响应的幅值
y_ff=[y_f(N/2+1:end),y_f(1:N/2)];%正负频率对应的响应的幅值
y_fre=ifft(y_ff,N);%逆FFT变换

end

function [t_4,y_digi] = digi_methodfor37(t_p,diat_t)
%中心差分法计算结构动力响应
%the parameter of structure
%t_p=10;
%diat_t=0.005;

m=545;%mass
k=5951.11111*10^3;%the stiffness (standard )
w_n=104.5;%the unit is radian/second natrual vibration
kesai=0.01;%the ratio of damping
c=2*m*w_n*kesai;%the coefficiency of damping
n=0:1:t_p/diat_t;
t_4=n*diat_t;
p0=267;%amplitude of loads (N)
theta=10*pi;%the frequency of loads (radian/s)
p=p0.*sin(theta*t_4);%function of loads

%初始条件（初位移、初速度）
y_digi(1)=0;
v_0=0;
%基本数据准备和初始条件计算
a_0=(1/m)*(p(1)-c*v_0-k*y_digi(1));
y_nt=y_digi(1)-diat_t*v_0+diat_t^2*a_0/2;
%计算等效刚度和中心差分计算公式中的系数
k_p=m/diat_t^2+c/2/diat_t;
a=k-2*m/diat_t^2;
b=m/diat_t^2-c/2/diat_t;
% 根据ti及ti以前时刻的运动，计算ti+1时刻的运动。
p_p(1)=p(1)-a*y_digi(1)-b*y_nt;% corresponding to the time of zero point
y_digi(2)=p_p(1)/k_p;%corresponding to the time of the first diat_t
for i=2:t_p/diat_t
p_p(i)=p(i)-a*y_digi(i)-b*y_digi(i-1);
y_digi(i+1)=p_p(i)/k_p;
end
end

function [t_b,y_t_b] = considering_bitemfor37(aa,bb)
%解析法
aa=6;%time range (s)
bb=0.005;% time interval
m=545;                      %mass (kg)
K=5.9511111*10^6;           %（kn/m)
kesi=0.01;                  %damping ratio
theta=10*pi;                %the frequency of loads (rad/s)
w_n=(K/m)^0.5;          % natural frequency of vibration
w_d=w_n*(1-kesi^2)^0.5;     %natural frequency of damping system
P_0=267;                    %amplitude of loads (N)
% P_t=P_0.*sin(theta*t);      %function of loads
y_st=(P_0/K);               %displacement (mm)
y_0=0;
v_0=0;

t_b=0:bb:aa;% time series
A=y_0;
B=(v_0/w_n)-y_st*(theta/w_n)/(1-(theta/w_n)^2);
C=y_st*(1-(theta/w_n)^2)/((1-(theta/w_n)^2)^2+(2*kesi*(theta/w_n))^2);
D=y_st*((-2)*kesi*(theta/w_n))/((1-(theta/w_n)^2)^2+(2*kesi*(theta/w_n))^2);
E=exp(-kesi*w_n*t_b);
F=A.*cos(w_d*t_b)+B.*sin(w_d*t_b);
G=E.*F;
H=C*sin(theta*t_b)+D*cos(theta*t_b);
y_t_b=G+H;%displacement
end
```
