## 2. 用MATLAB对杆系结构进行有限元分析。
[杆系结构的MATLAB有限元分析](https://stumail-my.sharepoint.cn/:b:/g/personal/19tgpu_stu_edu_cn/ES8RJ6a_P45GlhIHeq5qsugBRGK32HtCInoAZHAnnFCKxw?e=qvsXop) (附详细分析过程）

代码如下：
```Matlab
1、EXCEL表格输入结构信息
计算之前在同目录文件下的名为“imform_list.xls”的EXCEL表格中写入结构相关信息，
工作表中有四个工作簿，分别为node、element、force、constaint。

2、定义参数并读入EXCEL表格信息
%定义参数
E=2.06e5;%弹性模量
A=10;%截面面积
I=100;%截面惯性矩
[informlist_node]=xlsread('imform_list',1,'A2:I100')     
%从Excel表中读入节点信息。

[informlist_element]=xlsread('imform_list',2)                
%从Excel表中读入单元信息，包括单元号，单元连接的节点，节点坐标。

[informlist_force]=xlsread('imform_list',3)                    
%从Excel表中读入荷载信息。

[informlist_constraint]=xlsread('imform_list',4,'A2:A100')   
%从Excel表中读入有约束的位移列阵信息。

[informlist_unconstraint]=xlsread('imform_list',4,'B2:B100') 
%从Excel表中读入没有约束的位移列阵信息。

N_node=informlist_node(end,1)%节点总数
N_element=informlist_element(end,1)%单元总数
N_dis=max(max(informlist_node( : ,4:6)))%未知位移分量的总个数（包括边界约束处位移）
%绘制结构，以检查数据的正确性。
for i=1:N_element
    axis([min(min(informlist_node(1:end,2)))-0.2,...
        max(max(informlist_node(1:end,2)))+0.2,...
        min(min(informlist_node(1:end,3)))-0.2,...
        max(max(informlist_node(1:end,3)))+0.2]);
    line(informlist_element(i,[5,7]),informlist_element(i,[6,8]),'Color','k','LineStyle','-')
end

3、计算各单元转换矩阵
%计算各单元转换矩阵
T=cell(N_element,1) %定义一个元胞数组，用来存放转换矩阵。
for i=1:N_element
T(i)={Transtionmatrix(i,informlist_element)};
end

4、计算各单元的刚度矩阵
%计算各单元刚度矩阵（在整体坐标系下）
k=cell(N_element,1)       %首先定义一个元胞数组，用来装单元刚度矩阵。
for i=1:N_element
k(i)={element_Stiffness(E,A,I,informlist_element,i)}
end

5、建立整体刚度方程
%KK清零，然后调用函数Node _Assembly进行刚度矩阵的集成组装。
KK=zeros(N_dis)
for i=1:N_element
KK=element_Assembly(KK,k{i},informlist_element,i)
end

6、生成综合节点荷载列阵
Pe=EqLode(informlist_node,informlist_element) %整体坐标系下的非结点荷载
Pp=informlist_force  %结点荷载
P1=Pe+Pp%综合结点荷载

7、边界条件的处理及刚度方程求解
采用高斯消去法进行求解，注意：MATLAB中的反斜线符号“\”就是采用高斯消去法。
K=KK(informlist_unconstraint,informlist_unconstraint)
P=P1(informlist_unconstraint,1)
U=K\P     %位移
UU=zeros(N_dis,1);
UU(informlist_unconstraint)=U
q=zeros(N_node,3)     %所有的节点位移
for i=1:N_node
    if isnan(informlist_node(i,6))~=1
        q(i,1:3)=UU(informlist_node(i,4:6))
    else
        q(i,1:2)=UU(informlist_node(i,4:5))   
    end
end

8、位移变形图绘制
%变形之后的节点坐标
nodepos_xy=zeros(N_element,6)
for i=1:N_element
    nodepos_xy(i,1:3)=q(informlist_element(i,3),1:3)
    nodepos_xy(i,4:6)=q(informlist_element(i,4),1:3)
end
nodepos_xy1=nodepos_xy(:,[1,2,4,5])*(0.3/max(max(nodepos_xy)))+...
    informlist_element(1:end,5:8)  %变形之后的节点坐标列阵，其中变形位移放大了适当的倍数。


%绘制变形位移图
axis([min(min(informlist_node(1:end,2)))-1,max(max(informlist_node(1:end,2)))+...
    1,min(min(informlist_node(1:end,3)))-1,max(max(informlist_node(1:end,3)))+1]);
for i=1:N_element
line(nodepos_xy1(i,[1,3]),nodepos_xy1(i,[2,4]),'Color','r','LineStyle','--')
hold on
line(informlist_element(i,[5,7]),informlist_element(i,[6,8]),'Color','k','LineStyle','-')
end

8、计算各杆的内力绘制
注：各杆的杆端内力为局部坐标系下的。
FF=zeros(N_element,6) %定义一个矩阵，用来存储各单元内力
for i=1:N_element
    if informlist_element(i,2)==1
        FF(i,1:6)=T{i}*k{i}*nodepos_xy(i,1:6)' 
    else
        FF(i,[1,2,4,5])=T{i}*k{i}*nodepos_xy(i,[1,2,4,5])' 
    end      
end

函数定义部分
function T=Transtionmatrix(elementnumber,informlist_element)
%该函数计算单元的转换矩阵
%输入单元信息列阵informent_element、单元编号elementnumber
%输出单元转换矩阵
if informlist_element(elementnumber,2)==1
    L=sqrt((informlist_element(elementnumber,8)-informlist_element(elementnumber,6))^2+...
        (informlist_element(elementnumber,7)-informlist_element(elementnumber,5))^2);
    C=(informlist_element(elementnumber,7)-informlist_element(elementnumber,5))/L;
    S=(informlist_element(elementnumber,8)-informlist_element(elementnumber,6))/L;
    T=[C S 0 0 0 0; -S C 0 0 0 0;0 0 1 0 0 0;0 0 0 C S 0;0 0 0 -S C 0;0 0 0 0 0 1];
else
    L=sqrt((informlist_element(elementnumber,8)-informlist_element(elementnumber,6))^2+...
        (informlist_element(elementnumber,7)-informlist_element(elementnumber,5))^2);
    C=(informlist_element(elementnumber,7)-informlist_element(elementnumber,5))/L;
    S=(informlist_element(elementnumber,8)-informlist_element(elementnumber,6))/L;
    T=[C S 0 0; -S C 0 0;0 0 C S;0 0 -S C];
end
end

function k=element_Stiffness(E,A,I,informlist_element,elementnumber)
%该函数计算单元在整体坐标系下的刚度矩阵，根据单元类型分别计算。
%输入弹性模量E，横截面积A，单元信息列表informlist_element
%输出单元刚度矩阵，如果为链杆单元，单元刚度矩阵为k(4X4)，如果为梁单元，单元刚度矩阵为k(6X6)
L=sqrt((informlist_element(elementnumber,8)-informlist_element(elementnumber,6))^2+...
    (informlist_element(elementnumber,7)-informlist_element(elementnumber,5))^2);
if informlist_element(elementnumber,2)==1
    k=Transtionmatrix(elementnumber,informlist_element)'*...
        [E*A/L 0 0 -E*A/L 0 0;0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;...
        0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;-E*A/L 0 0 E*A/L 0 0;...
        0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;...
        0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L]*...
        Transtionmatrix(elementnumber,informlist_element);
else
    k=E*A/L*Transtionmatrix(elementnumber,informlist_element)'*...
        [1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0]*...
        Transtionmatrix(elementnumber,informlist_element);
end
end

function z=element_Assembly(KK,k,informlist_element,m)
%该函数进行单元刚度矩阵的组装
%输入单元刚度矩阵k
%输出整体刚度矩阵KK

if informlist_element(m,2)==1
    for n1=1:6
        for n2=1:6
        KK(informlist_element(m,n1+8),informlist_element(m,n2+8))=...
            KK(informlist_element(m,n1+8),informlist_element(m,n2+8))+k(n1,n2);
        end
    end
else
for n1=1:4
    for n2=1:4
    dof(1)=9;
    dof(2)=10;
    dof(3)=12;
    dof(4)=13;
    KK(informlist_element(m,dof(n1)),informlist_element(m,dof(n2)))=...
        KK(informlist_element(m,dof(n1)),informlist_element(m,dof(n2)))+k(n1,n2);
    end
end
end
z=KK;
end

function Pe=EqLode(informlist_node,informlist_element)
%该函数进行结构非节点荷载的等效节点荷载计算。
%输入单元信息列表，返回结构在整体坐标系下的等效节点荷载列阵。
Pe=zeros(max(max(informlist_node( : ,4:6))),1);
for i=1:informlist_element(end,1)
     P1=informlist_element(i,16);
     L1=sqrt((informlist_element(i,8)-informlist_element(i,6))^2+...
         (informlist_element(i,7)-informlist_element(i,5))^2);
     if informlist_element(i,2)==1 
         if informlist_element(i,15)==0
             Pe(informlist_element(i,9:14))=0;
         elseif informlist_element(i,15)==1
             Pe(informlist_element(i,9:14))=Pe(informlist_element(i,9:14))+...
                 Transtionmatrix(i,informlist_element)'*...
                 [0 P1/2 -P1*L1/8 0 P1/2 P1*L1/8]';
         else
             Pe(informlist_element(i,9:14))=Pe(informlist_element(i,9:14))+...
                 Transtionmatrix(i,informlist_element)'*...
                 [0 P1*L1/2 -P1*L1^2/12 0 P1*L1/2 P1*L1^2/12]';
         end
     else
         if informlist_element(i,15)==0
             Pe(informlist_element(i,[9 10 12 13]))=0;
         elseif informlist_element(i,15)==1
             Pe(informlist_element(i,[9 10 12 13]))=...
                 Pe(informlist_element(i,[9 10 12 13]))+...
                 Transtionmatrix(i,informlist_element)'*[0 P1/2 0 P1/2]';
         else
             Pe(informlist_element(i,[9 10 12 13]))=...
                 Pe(informlist_element(i,[9 10 12 13]))+...
                 Transtionmatrix(i,informlist_element)'*[0 P1*L1/2 0 P1*L1/2]';
         end
     end
end
end
```
