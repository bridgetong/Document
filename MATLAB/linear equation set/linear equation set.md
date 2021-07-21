```Matlab
##%% *计算方法作业/ 线性方程组求解 Matlab 代码*
% 
% %% 
%% 1、输入8个线性方程组的计算参数

A1=[4 -1;1 5];
B1=[15;9];

A2=[8 -3;-1 4];
B2=[10;6];

A3=[-1 3;6 -2];
B3=[1;2];

A4=[2 3;7 -2];
B4=[1;1];

A5=[5 -1 1;2 8 -1;-1 1 4];
B5=[10;11;3];

A6=[2 8 -1;5 -1 1;-1 1 4];
B6=[11;10;3];

A7=[1 -5 -1;4 1 -1;2 -1 -6];
B7=[-8;13;-2];

A8=[4 1 -1;1 -5 -1;2 -1 -6];
B8=[13;-8;-2];
%将计算参数放到一个数组，方便引用
A={A1,A2,A3,A4,A5,A6,A7,A8};
B={B1,B2,B3,B4,B5,B6,B7,B8};
%初始值
P1=[0;0];
P2=[0;0;0];

%% 2、线性方程组计算

for i=1:4
    disp(i)
    a=ifconvergence (A{i}) 
    X_J=Jacobi(A{i},B{i},P1,0.000001,100)  
    X_G=Gseid(A{i},B{i},P1,0.000001,100)
    X_LU=lufact(A{i},B{i})
    X_up=uptrbk(A{i},B{i})
end

for i=5:8
    disp(i)
    a=ifconvergence (A{i}) 
    X_J=Jacobi(A{i},B{i},P2,0.000001,100) 
    X_G=Gseid(A{i},B{i},P2,0.000001,100)
    X_LU=lufact(A{i},B{i})
    X_up=uptrbk(A{i},B{i})
end
%% 3、函数部分

function  a=ifconvergence (A)
 N=size(A,1);
 for j=1
 B_J(1,1)=0;
 B_J(1,j+1:N)=-A(1,j+1:N)/A(j,j);
 end
 
for j=1:N
     B_J(j,[1:j-1,j+1:N])=-A(j,[1:j-1,j+1:N])/A(j,j);
end
eig_value=eig(B_J)
rho=max(eig_value);

if (real(rho)^2+imag(rho)^2)^0.5<1
    a='此方程能用雅可比迭代和高斯赛德尔迭代';
else
    a='此方程不能用雅可比迭代和高斯赛德尔迭代';   
end
end
%% 
% 
function X=Jacobi(A,B,P,delta,max1)
% Input  A is an N*N nonsingular matrix
%        B is an N*1 matrix
%        P is an N*1 matrix, the initial guess
%        delta is the tolerance for P
%       max1 is the Maximum number of iterations

% Output X is an N*1 matrix containing the solution to AX=B, the Jacobi
% approximation to the solution of AX=B 

N=length(B);
for k=1:max1
    for j=1:N
     X(j)=(B(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);
    end
     
    err=abs(norm(X'-P));
    relerr=err/(norm(X)+eps);
    P=X';
    if (err<delta)|(relerr<delta)
        break
    end
end
X=X';
end
%% 
% 
function X=Gseid(A,B,P,delta,max1)
% Input  A is an N*N nonsingular matrix
%        B is an N*1 matrix
%        P is an N*1 matrix, the initial guess
%        delta is the tolerance for P
%       max1 is the Maximum number of iterations
% Output X is an N*1 matrix containing the solution to AX=B, the gauss-seid
%        approximation to the solution of AX=B 

N=length(B);

for k=1:max1
    for j=1:N
        if j==1
            X(1)=(B(1)-A(1,2:N)*P(2:N))/A(1,1);
        elseif j==N
             X(N)=(B(N)-A(N,1:N-1)*(X(1:N-1))')/A(N,N);
        else
            % X contains the kth approximations and P the (k-1)st
             X(j)=(B(j)-A(j,1:j-1)*X(1:j-1)'-A(j,j+1:N)*P(j+1:N))/A(j,j);
        end
    end
     
    err=abs(norm(X'-P));
    relerr=err/(norm(X)+eps);
    P=X';
    if (err<delta)|(relerr<delta)
        break
    end
end
X=X';
end
%% 
% 
function X=lufact(A,B)  %%带选主元的LU分解法

% Input  A is an N*N nonsingular matrix
%        B is an N*1 matrix
% Output X is an N*1 matrix containing the solution to AX=B

% Initialize X , Y , the temporary storage matrix C and the row permutation
% information matrix R

  [N,N]=size(A);
   X=zeros(N,1);
   Y=zeros(N,1);
   C=zeros(1,N+1);
   R=1:N;
   
 for p=1:N-1
  % Find the pivot row for column p
  [max1,j]=max(abs(A(p:N,p)));

 %Interchange row p and j
  C=A(p,:);
  A(p,:)=A(j+p-1,:);
  A(j+p-1,:)=C;
  d=(R(p));  %%记录变换行的顺序
  R(p)=R(j+p-1);
  R(j+p-1)=d;
  if  A(p,p)==0
      'A was singular No unique solution'
      break
  end
  
  %Calculate multiplier and place in subdiagomal portion of A
  
  for k=p+1:N
      mult=A(k,p)/A(p,p);
      A(k,p)=mult;                   %%向前求出一个l             
      A(k,p+1:N)=A(k,p+1:N)-mult*A(p,p+1:N);  %%先后求出临时的u
  end
 end
 
 %solve for Y
 
 Y(1)=B(R(1));
 
 for k=2:N
   Y(k)=B(R(k))-A(k,1:k-1)*Y(1:k-1);
 end
 
 
 % solve for X
 X(N)=Y(N)/A(N,N);
 
 for k=N-1:-1:1
   X(k)=(Y(k)-A(k,k+1:N)*X(k+1:N))/A(k,k);
 end
end
%% 
% 
function X=uptrbk(A,B)
% Input  A is an N*N nonsingular matrix
%        B is an N*1 matrix
% Output X is an N*1 matrix containing the solution to AX=B

% Initialize X and the temporary storage matrix C
  [N,N]=size(A);
   X=zeros(N,1);
   C=zeros(1,N+1);
 
%Form the augmented matric:Aug=[A|B]
   Aug=[A B];
 
 for p=1:N-1
     
     % Partial pivoting for column p
     [Y,j]=max(abs(Aug(p:N,p)));
     
     %Interchange row p and j
     C=Aug(p,:);
     Aug(p,:)=Aug(j+p-1,:);
     Aug(j+p-1,:)=C;
  if  Aug(p,p)==0
      %'A was singular No unique solution'
      break
  end
     
  % Elimination process for column p
    
     for k=p+1:N
         m=Aug(k,p)/Aug(p,p);
         Aug(k,p:N+1)=Aug(k,p:N+1)-m*Aug(p,p:N+1);
    % Aug(2,1:4)=Aug(2,1:4)-m*Aug(1,1:4)
     end
 end
 
%Back Substitution on [U|Y] using Program 3.1

  X=backsub(Aug(1:N,1:N),Aug(1:N,N+1));
end
```
