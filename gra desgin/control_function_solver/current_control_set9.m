function [I] = current_control_set9(par,mode,h,slope,Tchip_zero_lb,Tchip_zero_ub)
%% 本函数用于main9中沿程电流的优化设置，从主函数中读取基础参数后，配置相应的优化变量，优化约束，形成优化问题
%   输入：基础工况参数
%   输出：优化得出的电流
rng shuffle
%% 定义优化变量
lb = zeros(289,1);
ub = 1.5*ones(289,1);
x = optimvar('x',289,'LowerBound',lb,'UpperBound',ub);
%% 构建优化表达式

[data,i_final,break_flag,Tchip_record,grad_T_chip,Q_tec_overall,COP_record] = fcn2optimexpr(@(x)current_control_function_solver_fixed9(par,Tchip_zero_lb,Tchip_zero_ub,mode,h,slope,x),x);
%% 定义优化约束
%排除不可行解和无意义解(初值对不对还要检验)
limit1 = i_final == 289;
limit2 = break_flag == 0;
% %给出温度约束
% limit3 = 
% limit4 = max(Tchip_record) <= Tchip_zero_ub;
%% 定义优化目标
target = dot(grad_T_chip,grad_T_chip);

%% 设定问题，最小化target
prob = optimproblem('Objective',target);
%% 等式约束
%% 不等式约束
prob.Constraints.limit1 = limit1;
prob.Constraints.limit2 = limit2;
% prob.Constraints.limit3 = limit3;
% prob.Constraints.limit4 = limit4;
% prob.Constraints
% prob.Constraints
%% 查看可用求解器
[~,validsolvers] = solvers(prob);
disp(validsolvers);
rng shuffle
x0.x = (ub-lb).*rand(289,1)+lb;
disp(x0.x(1))
Options = optimoptions('fmincon','Display','iter');
[sol,fval,exitflag] = solve(prob,x0,'Options',Options);

for i = 1:289
    I(i) = sol.x(i);
end


end