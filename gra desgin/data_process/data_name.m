function [data_name] = data_name(par,exist,mode,h,slope)
%此函数用来给定任意工况下的数据名字
%   此处显示详细说明
%读取工况数据
v = par.mg;
Tin = par.mg;
material = par.material;
P = par.P_chip;
fprate = h+1;
if exist == 1
    current_mode = 'ctrl ';
else 
    current_mode = 'zero ';
end
switch mode
    case 0 
        mode_name = 'uni';
    case 1
        mode_name = 'entr';
    case 2
        mode_name = 'mid';
    case 3
        mode_name = 'exit';
end
data_name = [current_mode mode_name ' Mtl_' material ' v_' num2str(v) ' Tin_' num2str(Tin)...
    ' P_' num2str(P) ' Pkrate_' num2str(fprate) ' slope_' num2str(slope)];

%对名字中的小数点进行处理
data_name = replace(data_name,'.','_');
end