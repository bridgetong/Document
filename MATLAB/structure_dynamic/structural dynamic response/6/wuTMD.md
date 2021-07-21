```
%% 结构动力学作业4-7地震作用下结构响应
%杜哈梅积分和傅里叶变换都只能求稳态解，求出解后加入非零初始条件（速度和位移为零，加速度可不为零）可额外求出瞬态解,最后形成完整解。数值方法直接求的完整解。
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

%=====================地震作用下的结构位移计算(三种方法)=======================
pe=-diag(M)*EQA;%地震等效作用力
for l=1:n;
MG(l)=model(:,l)'*M*model(:,l);%计算广义质量
KG(l)=MG(l).*eigvalue(l).^2;%计算广义刚度
end
PG=model'*pe;  %计算广义荷载（一个横向量对应一个阶数）
wd=eigvalue.*sqrt(1-kesai.^2);
%----------------------1)杜哈梅积分法--------------------------
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

%----------------------2)快速傅里叶变换法-----------------------
nn=0:1:aa/diata-1;%时间步
t_3=nn*diata;%持续时间
N=aa/diata;%离散数
for l=1:n
p_f(l,:)=fft(PG(l,:),N);%fft变换后包含了正负频率的完整数据，正负频率分别存储在结果的前后半部分，matlab在振动信号中的应用一书中有讲。
%fft中N为采样频率，其为最大识别频率的2倍,所以fft变换后只有-0.5N~0.5N频率范围的数据.数据量n应大于N，这样频谱分辨率N/n小于1HZ从而更加精确，即为什么要在数据后面加很多零.
end
p_ff=[p_f(:,N/2+1:end),p_f(:,1:N/2)];%after fft transform，改成负、正频率对应的幅值
f=(-N/2:1:N/2-1)/(N*diata);%频率是从-N/2 到正N/2，
beita=2.*f.*pi./wd;%频率比 
for l=1:n
fre_res(l,:)=(1/KG(l))*((1-beita(l,:).^2)+i*2*kesai(l)*beita(l,:)).^-1;%频响函数也是负正频率对应的值
end
y_f=p_ff.*fre_res;%响应的频谱
y_ff=[y_f(:,N/2+1:end),y_f(:,1:N/2)];%调整正负频率对应的响应的频谱
for l=1:n
    y_fr(l,:)=ifft(y_ff(l,:),N);%逆FFT变换
end
D_fre=real(y_fr);%取实部
y_fre=model*D_fre;

%---------------------3)数值方法(中心差分法）_______________________
nnn=0:1:aa/diata;
t_4=nnn*diata;
%the initial condition
v_0=0;
for l=1:n;
Ddig(l,1)=0;
a_0(l)=(1/MG(l))*(PG(l,1)-2*kesai(l)*eigvalue(l)*v_0-eigvalue(l)^2*MG(l)*Ddig(l,1));
Ddigg(l)=Ddig(l,1)-diata*v_0+diata^2*a_0(l)/2;
end
% %the parameter
for l=1:n
k_p(l)=MG(l)/diata^2+kesai(l)*eigvalue(l)*MG(l)/diata;
a(l)=eigvalue(l)^2*MG(l)-2*MG(l)/diata^2;
b(l)=MG(l)/diata^2-kesai(l)*eigvalue(l)*MG(l)/diata;
p_p(l,1)=PG(l,1)-a(l)*Ddig(l,1)-b(l)*Ddigg(l);% corresponding to the time of zero point
Ddig(l,2)=p_p(l,1)/k_p(l);

for j=2:aa/diata
p_p(l,j)=PG(l,j)-a(l)*Ddig(l,j)-b(l)*Ddig(l,j-1);
Ddig(l,j+1)=p_p(l,j)/k_p(l);
end
end
y_dig=model*Ddig;

%=====================绘图=======================
subplot(3,1,1);
hold on
plot(t,y(1,:),':go')
plot(t_3,y_fre(1,:),'--ks')
plot(t_4,y_dig(1,:),'-.rp')
legend('Duhamel-interal','Frequency-domain-method','numerical-method');%,'Frequency-domain-half'
hold off
xlabel('t/s');
ylabel('位移/m');
title('质点1位移图')

subplot(3,1,2);
hold on
plot(t,y(2,:),':go')
plot(t_3,y_fre(2,:),'--ks')
plot(t_4,y_dig(2,:),'-.rp')
legend('Duhamel-interal','Frequency-domain-method','numerical-method');%,'Frequency-domain-half'
hold off
xlabel('t/s');
ylabel('位移/m');
title('质点2位移图')

subplot(3,1,3);
hold on
plot(t,y(3,:),':go')
plot(t_3,y_fre(3,:),'--ks')
plot(t_4,y_dig(3,:),'-.rp')
legend('Duhamel-interal','Frequency-domain-metsahod','numerical-method');%,'Frequency-domain-half'
hold off
xlabel('t/s');
ylabel('位移/m');
title('质点3位移图')
```
