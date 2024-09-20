function [P] = heat_change(n)
%HEAT_CHANGE 此处显示有关此函数的摘要
%   此处显示详细说明
%生成热流峰值的坐标
par = para;
% x = [4,4,9,9,14,14];
% y = [5,13,5,13,5,13];
x = [14];
y = [14];
%生成山峰的高度
h = ones(1,1);
%生成山峰的坡度
xs = 2*ones(1,1);
ys = 2*ones(1,1);
%生成地形数据
Z = ones(17,17);
[X,Y] = meshgrid(1:17);
for k1 = 1
    Z =  Z + h(k1)*exp(-((X - x(k1))/xs(k1)).^2 - ((Y - y(k1))/ys(k1)).^2);
end
P_total = par.P_chip*17*17;
P_tem = sum(Z(:));
Z = Z*P_total/P_tem;
%求行数
Power = [];
k = 1;
for i = 1:17
    for j = 1:17
        if mod(i,2)==1
            Power(k) = Z(i,j);
        else
            Power(k) = Z(i,18-j);
        end
        k = k+1;
    end
end
Power(290) = Power(289);


P = Power(n);
% 均匀热流更改
P = 6;
end




