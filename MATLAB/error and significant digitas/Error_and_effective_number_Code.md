```matlab
&计算方法第二次作业—舍入误差与有效数字  程序示例| 
%%计算通用程序

clear;
N=input('请输入一个N值  N=:');
accurate_value=double((3/2-1/N-1/(N+1))/2);%精确值采用双精度
Sn_StL=single(0);
Sn_LtS=single(0);
for i=2:N
    Sn_StL=Sn_StL+1/(i*i-1);
end
for i=N:-1:2
    Sn_LtS=Sn_LtS+1/(i*i-1);
end
fprintf('精确值为：%17.15f\n',accurate_value);%精确值保留小数点后15位
fprintf('从大到小的顺序累加结果为SN=%f\n',Sn_StL);
fprintf('从小到大的顺序累加结果为SN=%f\n',Sn_LtS);
```
