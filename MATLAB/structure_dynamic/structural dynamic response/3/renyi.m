%% �ṹ����ѧ3-7��������µ���Ӧ
%�Ź�÷���ֺ͸���Ҷ�任��ֻ������̬�⣬�������������ʼ�������ٶȺ�λ��Ϊ�㣬���ٶȿɲ�Ϊ�㣩�ɶ������˲̬��,����γ������⡣��ֵ����ֱ����������⡣
clear
clc
m=importdata('C:\Users\13718\Desktop\ffffffffffffff\jiegou\0 (Ve).ap');
inpu_dat=(m.data(:,1)*1000)';
[t,y_du]=Duhamelintegralfor37(50,1/625,inpu_dat);% solution in time domain Duhamel integral
[t_4,y_digi]=digi_methodfor37(50,1/625,inpu_dat);%numerical solution
[t_3,y_fre]=fre_methodfor37(60,1/625,inpu_dat);%fft and ifft
figure(1) 
subplot(2,1,1)
plot(t,inpu_dat(1:50*625+1));
hold off
subplot(2,1,2)
plot(t,y_du,':go');
hold on
plot(t_4,y_digi,'--ks')
hold on
plot(t_3,real(y_fre),'-.rp')

legend('Duhamel-interal','numerical-method','Frequency-domain-method');%,'Frequency-domain-half'
hold off



function [t,y_du] = Duhamelintegralfor37(aa,bb,inpu_dat)
%inputting paremeters.
% aa=10;up bound of time
%bb=0.02;% time interval
%%%%%
t=0:bb:aa;
t1=t;% time series
%%%%%
theta=10*pi;%the frequency of loads (radian/s)
w=104.5;% natural frequency of vibration (radian/s)
kesi=0.01;%damping ratio
w_d=w*(1-kesi^2)^0.5;% natural frequency of damping system
m=545;% mass (kg)
%%%%% 
p=inpu_dat;%function of loads
 
%%%%%%%%the main program
for i=1:(aa/bb+1)
    for j=1:i
        %paramet1(j)=p(j)/(m*w)*bb*sin(w*(t(i)-t1(j)));%Duhameil's integral function for non-damping system
        paramet2(j)=p(j)/(m*w_d)*bb*exp(-kesi*w*(t(i)-t1(j)))*sin(w_d*(t(i)-t1(j)));%%Duhameil's integral function for damping system
       
    end
    y_du(i)=sum(paramet2);%%displacement 
   
end

end

function [t_3,y_fre] = fre_methodfor37(t_p,diata_t,inpu_dat)
%t_p=10;
%diata_t=0.01;
nn=0:1:t_p/diata_t-1;
t_3=nn*diata_t;
N=t_p/diata_t;
w=10*pi;%the unit is radian/second load
p_t=inpu_dat;
p_f=fft(p_t,N);%fft�任�����������Ƶ�ʵ��������ݣ�����Ƶ�ʷֱ�洢�ڽ����ǰ��벿�֣�matlab�����ź��е�Ӧ��һ�����н���
p_ff=[p_f(N/2+1:end),p_f(1:N/2)];%after fft transform���ĳɸ�����Ƶ�ʶ�Ӧ�ķ�ֵ
f=(-N/2:1:N/2-1)/(N*diata_t);%Ƶ���Ǵ�-N/2 ����N/2��
%plot(f(1:N/2),abs(y(1:N/2))*2/N)
k=5951.11111*10^3;%the standard unit
w_n=104.5;%the unit is radian/second natrual vibration
kesai=0.01;
beita=2*f*pi/w_n;
fre_res=(1/k)*((1-beita.^2)+i*2*kesai*beita).^-1;%Ƶ�캯��Ҳ�Ǹ���Ƶ�ʶ�Ӧ��ֵ
y_f=p_ff.*fre_res;%��Ӧ�ķ�ֵ
y_ff=[y_f(N/2+1:end),y_f(1:N/2)];%����Ƶ�ʶ�Ӧ����Ӧ�ķ�ֵ
y_fre=ifft(y_ff,N);%��FFT�任
end
%%remark
%��������㹻�����Ƶ��õ��Ľ�����������ַ����õ��Ľ����ͬ������ʼ��һ��ʱ�䲻��ͬ��

function [t_4,y_digi] = digi_methodfor37(t_p,diat_t,inpu_dat)
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
% p0=267;%amplitude of loads (N)
theta=10*pi;%the frequency of loads (radian/s)
p=inpu_dat;%function of loads
%the initial condition
y_digi(1)=0;
v_0=0;
a_0=(1/m)*(p(1)-c*v_0-k*y_digi(1));
y_nt=y_digi(1)-diat_t*v_0+diat_t^2*a_0/2;
%the parameter
k_p=m/diat_t^2+c/2/diat_t;
a=k-2*m/diat_t^2;
b=m/diat_t^2-c/2/diat_t;
% the 
p_p(1)=p(1)-a*y_digi(1)-b*y_nt;% corresponding to the time of zero point
y_digi(2)=p_p(1)/k_p;%corresponding to the time of the first diat_t
for i=2:t_p/diat_t
p_p(i)=p(i)-a*y_digi(i)-b*y_digi(i-1);
y_digi(i+1)=p_p(i)/k_p;
end
end