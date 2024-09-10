function [result] = flag_shift(flag_tem)
%FLAG_SHIFT 此处显示有关此函数的摘要
%   此处显示详细说明
a = 1-sign(flag_tem);
result = 500^a-1;
end

