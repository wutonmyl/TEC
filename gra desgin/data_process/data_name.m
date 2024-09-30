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
    current_mode = '调控数据 ';
else 
    current_mode = '零电流数据 ';
end
switch mode
    case 0 
        mode_name = '均匀热流';
    case 1
        mode_name = '入口热点';
    case 2
        mode_name = '中心热点';
    case 3
        mode_name = '出口热点';
end
data_name = [current_mode mode_name ' 材料_' material ' v_' num2str(v) ' Tin_' num2str(Tin)...
    ' P_' num2str(P) ' 峰平比_' num2str(fprate) ' 热流坡度_' num2str(slope)];

%对名字中的小数点进行处理
data_name = replace(data_name,'.','_');
end