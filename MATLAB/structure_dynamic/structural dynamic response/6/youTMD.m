%% �ṹ����ѧ��ҵ4-7��װ��TMD�ڵ��������µ���Ӧ
%�Ź�÷���ֺ͸���Ҷ�任��ֻ������̬�⣬�������������ʼ�������ٶȺ�λ��Ϊ�㣬���ٶȿɲ�Ϊ�㣩�ɶ������˲̬��,����γ������⡣��ֵ����ֱ����������⡣
clear
clc
%===========================������������ʹ���============================
%----------------����ʱ������------------------
m=importdata('C:\Users\13718\Desktop\ffffffffffffff\jiegou\TH.txt');
mm=importdata('C:\Users\13718\Desktop\ffffffffffffff\jiegou\TS.txt');
EQD=m(:,1)';%����λ������
EQT=mm(:,1)';%����ʱ������
diata=EQT(2)-EQT(1);%����ʱ����
%-----------������ٶȼ���(���Ĳ�ַ�)----------
EQD0=0;%��������ʼ����(���ٶȡ��ٶȡ�λ��Ϊ�㣩
EQA0=0;
EQV0=0;
EQDI=EQD0-diata*EQV0+diata^2*EQA0/2;%����-1ʱ�̵�λ��
EQA=(EQD(2)-2*EQD(1)+EQDI)/diata^2;%����ѭ��i=1ʱ�ĵ�����ٶ�
for l=2:18000%����ѭ��i=1�Ժ�ĵ�����ٶ�
EQA(l)=(EQD(l+1)-2*EQD(l)+EQD(l-1))/diata^2;%
end

%=============================�ṹģ̬����==============================
%-----------�������ͺ͹���Ƶ�ʼ���------------
k1=1e8;
k2=k1;
k3=k1;
k4=k1
m1=160e3;
m2=m1;
m3=80e3;
m4=(m1+m2+m3)*0.03;
kesai4=0.05;
M=[m1,0,0,0;0,m2,0,0;0,0,m3,0;0,0,0,m4]; %������������
K=[k1+k2,-k2,0,0;-k2,k2+k3,-k3,0;0,-k3,k3+k4,-k4;0,0,-k4,k4]; %����նȾ���
kesai=[0.05;0.05];  %����ǰ���׺�TMD�����
kesaiTMD=0.05
n=length(M);  %ά��
[v,d]=eig(inv(K)*M);  %�����ƽ�ⷽ��������ֵ�������������������Ƶ�ʵ�����ƽ�����ԽǾ��󣩺����;�����������
eigvalue_1=(1./diag(d)).^0.5; %�����Ƶ��(ԲƵ��)����������
[eigvalue,IJ]=sort(eigvalue_1,1,'ascend'); %�Թ���Ƶ����������,�����ԭ��λ�����
eigvalue
for l=1:n;
model(:,l)=v(:,IJ(l))./(v(n,IJ(l)));  %���;����һ��,��������Ƶ��һһ��Ӧ��
end 
model
%----------------��������ȼ���----------------
dam_a=inv([1/eigvalue(1),eigvalue(1);1/eigvalue(2),eigvalue(2)])*kesai*2;   %�������������ϵ��
c=dam_a(1)*M+dam_a(2)*K;  %��������������(�Ǳ�Ҫ��
for j=1:n
    kesai(j)=dam_a(1)/2/eigvalue(j)+dam_a(2)*eigvalue(j)/2;  %��֪ǰ��������ȣ������������(������Ϊ�Ǿ�ȷֵ��TMD������ȿ��þ�ȷֵ�滻��
end
kesai(4)=kesaiTMD;%TMD����Ⱦ�ȷֵ�滻�Ǿ�ȷֵ
kesai%Ϊʲô�����Ϊ�������뿴��34��
%=====================���������µĽṹλ�Ƽ���(���ַ���)=======================
pe=-diag(M)*EQA;%�����Ч������
for l=1:n;
MG(l)=model(:,l)'*M*model(:,l);%�����������
KG(l)=MG(l).*eigvalue(l).^2;%�������ն�
end
PG=model'*pe;  %���������أ�һ����������Ӧһ��������
wd=eigvalue.*sqrt(1-kesai.^2);
%----------------------1)�Ź�÷���ַ�--------------------------
MGG=diag(MG);%������������
aa=diata*3000;%ʱ������
t=0:diata:aa;
t1=t;% time series
for l=1:(aa/diata+1)
    for j=1:l
 paramet2(:,j)=PG(:,j)./(MG'.*wd).*diata.*exp(-kesai.*eigvalue.*(t(l)-t1(j))).*sin(wd.*(t(l)-t1(j)));%��Ӧʱ��Ϊt(i)ʱ������Ź���������ɢʱ��ι�������
    end
    dd=0;
    for z=1:l
        dd=dd+paramet2(:,z); %��Ӧʱ��Ϊt��i)ʱ����������ɢʱ��εĹ�������ֵ���
    end
    D(:,l)=dd; %(һ����������Ӧһ��������
end
y=model*D; %���͵���,һ����������Ӧһ����������˶�.

%----------------------2)���ٸ���Ҷ�任��-----------------------
nn=0:1:aa/diata-1;%ʱ�䲽
t_3=nn*diata;%����ʱ��
N=aa/diata;%��ɢ��
for l=1:n
p_f(l,:)=fft(PG(l,:),N);%fft�任�����������Ƶ�ʵ��������ݣ�����Ƶ�ʷֱ�洢�ڽ����ǰ��벿�֣�matlab�����ź��е�Ӧ��һ�����н���
%fft��NΪ����Ƶ�ʣ���Ϊ���ʶ��Ƶ�ʵ�2��,����fft�任��ֻ��-0.5N~0.5NƵ�ʷ�Χ������.������nӦ����N������Ƶ�׷ֱ���N/nС��1HZ�Ӷ����Ӿ�ȷ����ΪʲôҪ�����ݺ���Ӻܶ���.
end
p_ff=[p_f(:,N/2+1:end),p_f(:,1:N/2)];%after fft transform���ĳɸ�����Ƶ�ʶ�Ӧ�ķ�ֵ
f=(-N/2:1:N/2-1)/(N*diata);%Ƶ���Ǵ�-N/2 ����N/2��
beita=2.*f.*pi./wd;%Ƶ�ʱ�(���ף�
for l=1:n
fre_res(l,:)=(1/KG(l))*((1-beita(l,:).^2)+i*2*kesai(l)*beita(l,:)).^-1;%Ƶ�캯��Ҳ�Ǹ���Ƶ�ʶ�Ӧ��ֵ(������)
end
y_f=p_ff.*fre_res;%��Ӧ��Ƶ��
y_ff=[y_f(:,N/2+1:end),y_f(:,1:N/2)];%��������Ƶ�ʶ�Ӧ����Ӧ��Ƶ��
for l=1:n
    y_fr(l,:)=ifft(y_ff(l,:),N);%��FFT�任
end
D_fre=real(y_fr);%ȡʵ��
y_fre=model*D_fre;

%---------------------3)��ֵ����(���Ĳ�ַ���_______________________
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

%=====================��ͼ=======================
subplot(3,1,1);
hold on
plot(t,y(1,:),':go')
plot(t_3,y_fre(1,:),'--ks')
plot(t_4,y_dig(1,:),'-.rp')
legend('Duhamel-interal','Frequency-domain-method','numerical-method');%,'Frequency-domain-half'
hold off
xlabel('t/s');
ylabel('λ��/m');
title('�ʵ�1λ��ͼ(TMD)')

subplot(3,1,2);
hold on
plot(t,y(2,:),':go')
plot(t_3,y_fre(2,:),'--ks')
plot(t_4,y_dig(2,:),'-.rp')
legend('Duhamel-interal','Frequency-domain-method','numerical-method');%,'Frequency-domain-half'
hold off
xlabel('t/s');
ylabel('λ��/m');
title('�ʵ�2λ��ͼ(TMD)')

subplot(3,1,3);
hold on
plot(t,y(3,:),':go')
plot(t_3,y_fre(3,:),'--ks')
plot(t_4,y_dig(3,:),'-.rp')
legend('Duhamel-interal','Frequency-domain-metsahod','numerical-method');%,'Frequency-domain-half'
hold off
xlabel('t/s');
ylabel('λ��/m');
title('�ʵ�3λ��ͼ(TMD)')