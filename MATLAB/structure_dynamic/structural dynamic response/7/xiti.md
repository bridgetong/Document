```
%% 结构动力学作业4-7内力计算
clear
clc
%===========================地震作用输入和处理============================
%----------------地震时程输入------------------
m=importdata('C:\Users\13718\Desktop\ffffffffffffff\jiegou\TH.txt');
mm=importdata('C:\Users\13718\Desktop\ffffffffffffff\jiegou\TS.txt');
EQD=m(:,1)';%地震位移向量
EQT=mm(:,1)';%地震时间向量
diata=EQT(2)-EQT(1);%地震时间间隔
%-----------地震加速度计算(中心差分法)----------
EQD0=0;%定义地震初始条件(加速度、速度、位移为零）
EQA0=0;
EQV0=0;
EQDI=EQD0-diata*EQV0+diata^2*EQA0/2;%计算-1时刻的位移
EQA=(EQD(2)-2*EQD(1)+EQDI)/diata^2;%计算循环i=1时的地震加速度
for l=2:18000%计算循环i=1以后的地震加速度
EQA(l)=(EQD(l+1)-2*EQD(l)+EQD(l-1))/diata^2;%
end

%=============================结构模态计算==============================
%-----------各阶振型和固有频率计算------------
k1=1e8;
k2=k1;
k3=k1;
m1=160e3;
m2=m1;
m3=80e3;
M=[m1,0,0;0,m2,0;0,0,m3]; %输入质量矩阵
K=[k1+k2,-k2,0;-k2,k2+k3,-k3;0,-k3,k3]; %输入刚度矩阵
kesai=[0.05;0.05];  %输入前两阶阻尼比
n=length(M);  %维数
[v,d]=eig(inv(K)*M);  %对柔度平衡方程求特征值和特征向量，输出固有频率导数的平方（对角矩阵）和振型矩阵（列向量）
eigvalue_1=(1./diag(d)).^0.5; %求固有频率(圆频率)（列向量）
[eigvalue,IJ]=sort(eigvalue_1,1,'ascend'); %对固有频率向量排序,并输出原来位置序号
eigvalue
for l=1:n;
model(:,l)=v(:,IJ(l))./(v(n,IJ(l)));  %振型矩阵归一化,并排序与频率一一对应。
end
model
%----------------各阶阻尼比计算----------------
dam_a=inv([1/eigvalue(1),eigvalue(1);1/eigvalue(2),eigvalue(2)])*kesai*2;   %求瑞利阻尼矩阵系数
c=dam_a(1)*M+dam_a(2)*K;  %求出瑞利阻尼矩阵(目的是为了求各阶阻尼比,阻尼比才有用)
for j=1:n
    kesai(j)=dam_a(1)/2/eigvalue(j)+dam_a(2)*eigvalue(j)/2;  %已知前两阶阻尼比，求所有阻尼比
end
kesai

%=====================地震作用下的结构位移计算=======================
pe=-diag(M)*EQA;%地震等效作用力
for l=1:n;
MG(l)=model(:,l)'*M*model(:,l);%计算广义质量
KG(l)=MG(l).*eigvalue(l).^2;%计算广义刚度
end
PG=model'*pe;  %计算广义荷载（一个横向量对应一个阶数）
wd=eigvalue.*sqrt(1-kesai.^2);
%----------------------杜哈梅积分法求响应--------------------------
MGG=diag(MG);%广义质量矩阵
aa=diata*3000;%时间上限
t=0:diata:aa;
t1=t;% time series
for l=1:(aa/diata+1)
    for j=1:l
 paramet2(:,j)=PG(:,j)./(MG'.*wd).*diata.*exp(-kesai.*eigvalue.*(t(l)-t1(j))).*sin(wd.*(t(l)-t1(j)));%对应时间为t(i)时，计算杜哈曼单个离散时间段广义坐标
    end
    dd=0;
    for z=1:l
        dd=dd+paramet2(:,z); %对应时间为t（i)时，将所有离散时间段的广义坐标值相加
    end
    D(:,l)=dd; %(一个横向量对应一个阶数）
end
y=model*D; %振型叠加,一个横向量对应一个质量点的运动.

%=====================结构内力计算（剪力，弯矩)(忽略阻尼力）=======================
%---------------------剪力计算---------------------
f=K*y;%各质点的弹性力
for l=1:n;
    cc=0;
    for j=l:n
    cc=cc+f(j,:);%由各质点的弹性力计算各层楼的剪力过程
    end
    F(l,:)=cc;%各层楼的剪力(同一层楼不同高度处剪力都相同）
end
 %---------------------弯矩计算---------------------
 moment(1,:)=[3,6,9]*f;%为该层楼的最低点
 moment(2,:)=[0,3,6]*f;
 moment(3,:)=[0,0,3]*f;
 %=====================绘图=======================
subplot(3,2,1);
plot(t,F(1,:))
legend('剪力');
xlabel('t/s');
ylabel('剪力/N');
title('质点1剪力变化图')

subplot(3,2,2);
plot(t,moment(1,:))
legend('弯矩');
xlabel('t/s');
ylabel('弯矩/N*m');
title('质点1弯矩变化图')

subplot(3,2,3);
plot(t,F(2,:))
legend('剪力');
xlabel('t/s');
ylabel('剪力/N');
title('质点2剪力变化图')

subplot(3,2,4);
plot(t,moment(2,:))
legend('弯矩');
xlabel('t/s');
ylabel('弯矩/N*m');
title('质点2弯矩变化图')

subplot(3,2,5);
plot(t,F(3,:))
legend('剪力');
xlabel('t/s');
ylabel('剪力/N');
title('质点3剪力变化图')

subplot(3,2,6);
plot(t,moment(3,:))
legend('弯矩');
xlabel('t/s');
ylabel('弯矩/N*m');
title('质点3弯矩变化图')
```
