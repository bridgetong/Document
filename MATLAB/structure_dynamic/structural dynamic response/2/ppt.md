```matlab
%% 结构动力学作业(有阻尼简协振动）PPT上的例题
clear
clc
% 1 稳态解
w=83.78;%荷载频率
l=0.2;%阻尼比
p=30;%外荷载幅值
k=1314.5*10^3;%刚度
wn=91.79;%自振频率
b=w/wn;%频率比
f=-atan((2*l*b)/(1-b^2));%象位角
yst=p/k;%静位移
Rd=1/sqrt((1-b^2)^2+(2*l*b)^2);%动力放大系数
interval_t=2*pi/w/50;%时间间隔
t=0:interval_t:10;
A=yst.*Rd*sin(w.*t-f).*1000;%稳态解或特解
% 2 瞬态解(有阻尼自由振动）
y0=0.01;v0=0;%初始条件
wd=wn.*sqrt(1-l^2);
C1=(v0+l.*wn.*y0)/wd;
C2=y0;
interval_t1=2.*pi/wd/50
t1=0:interval_t1:10;
y=(exp(-l.*w.*t1)).*(C1.*sin(wd.*t1)+C2.*cos(wd.*t1))*1000;
% 3 绘图
hold on
plot(t,A);
xlabel('t/s')
ylabel('u/mm')
plot(t1,y,':');
legend('稳态反应','瞬态反应');
hold off
```
