function [P] = heat_change(n,mode,h,slope)
%HEAT_CHANGE 此处显示有关此函数的摘要
%   输入单元数n,非均匀模式mode，热流高度h，热流梯度slope
%生成热流峰值的坐标
par = para;
% x = [4,4,9,9,14,14];
% y = [5,13,5,13,5,13];
switch mode
    case 1
        %入口热点
        x = [4];
        y = [4];
    
    case 2
        %中心热点
        x = [9];
        y = [9];
    case 3
        %出口热点
        x = [14];
        y = [14];
    otherwise
        %无热点
        x = [];
        y = [];
end
num_point = length(x);
%生成山峰的高度
h = h*ones(num_point,1);
%生成山峰的坡度
xs = slope*ones(num_point,1);
ys = slope*ones(num_point,1);
%生成地形数据
Z = ones(17,17);
[X,Y] = meshgrid(1:17);

for k1 = 1:num_point
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
if mode == 0
    P = par.P_chip;
end

end




