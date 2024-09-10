function [x] = judge_wan(j)
%JUDGE_WAN 此处显示有关此函数的摘要
%   此处显示详细说明
flag_wan_record = [];
for i = 1:290
    flag_wan = 0;
    if mod(i,17)==1||mod(i,17)==0||mod(i,17)==2
        flag_wan = 1;
    end
    flag_wan_record = [flag_wan_record,flag_wan];
end
% x = flag_wan_record(j);
x = 0;
end

