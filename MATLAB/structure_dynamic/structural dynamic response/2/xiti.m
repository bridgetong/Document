%% 习题3.7
clear
clc
%稳态反应
w=31.42;%荷载频率
l=0.01;%阻尼比
p=267/2;%外荷载幅值
k=1/(3.36071*10^-7);%刚度
wn=104.4965;%自振频率
b=w/wn;%频率比
f=-atan((2*l*b)/(1-b^2));%象位角
yst=p/k;%静位移
Rd=1/sqrt((1-b^2)^2+(2*l*b)^2);%动力放大系数
interval_t=2*pi/w/20;%时间间隔
t=0:interval_t:10;
A=yst.*Rd*sin(w.*t-f); %位移
a=-A*w^2;%加速度
%瞬态反应
y0=0.05;v0=0;%初始条件
wd=wn.*sqrt(1-l^2);
C1=(v0+l.*wn.*y0)/wd;
C2=y0;
interval_t1=2.*pi/20;
t1=0:interval_t1:100;
y=(exp(-l.*w.*t1)).*(C1.*sin(wd.*t1)+C2.*cos(wd.*t1));
subplot(3,1,1);
plot(t,A);
xlabel('t/s');
ylabel('u/m');
set(gca,'ytick',-1*10^-4:5*10^-5:1*10^-4);
title('稳态位移图')
subplot(3,1,2);
plot(t,a);
xlabel('t/s');
ylabel('a/(m/s^2)');
set(gca,'ytick',-0.1:0.05:0.1);
title('稳态加速度图')
subplot(3,1,3);
plot(t1,y);
xlabel('t/s');
ylabel('u/m');
title('瞬态位移图')
