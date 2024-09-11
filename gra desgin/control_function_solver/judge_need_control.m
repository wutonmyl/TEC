function [need_control,Tchip_target_recmd] = judge_need_control(Tchip_record,i,j,i_special)
%   判断是否需要采用独立调控,如果需要就给出合适的调控目标
%   此处显示详细说明
grad_bond = 56.83*0.3;
Tchip_record_check = Tchip_record(1:i,j);
grad_check = gradient(Tchip_record_check);
prop = 3;%abs(grad_check(i,j));
%这里存在隐患，如果入口就大于300咋整
if Tchip_record_check(i,j) >= 330
    need_control = 2;
    %采用上一单元的温度
    Tchip_target_recmd = Tchip_record(i-1,j);
elseif (grad_check(i,j)) > grad_bond
    need_control = 1;
    %采用两单元间的中间值降低梯度
    Tchip_target_recmd = (min(Tchip_record(i_special,j),Tchip_record(i,j)) + abs(Tchip_record(i,j)-Tchip_record(i_special,j))/prop);
elseif grad_check(i,j) < -grad_bond
    need_control = -1;
    Tchip_target_recmd = (max(Tchip_record(i_special,j),Tchip_record(i,j)) - abs(Tchip_record(i,j)-Tchip_record(i_special,j))/prop);
else
    need_control = 0;
    Tchip_target_recmd = Tchip_record(i,j);
end