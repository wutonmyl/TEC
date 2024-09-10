function [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed3(par,T_g_in,i,j,T_chip_taget,I_record)
%   这里要把芯片温度变成电流
%   将方程组彻底变为优化问题,让功率器件温度变为优化目标
rng default
% 定义优化变量
lb = [-10,100,T_chip_taget,185,185,185,T_chip_taget-50];
ub = [10,400,800,700,600,600,600];
x = optimvar('x',7,'LowerBound',lb,'UpperBound',ub);
%第一个式子关于pchip还需要修改
exp1 = heat_change(i)-par.k_ct*par.a_te*(x(7)-x(2));
exp2 = par.n*(par.alpha*x(1)*x(2)-0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      par.k_ct*par.a_te*(x(7)-x(2));
exp3 = par.n*(par.alpha*x(1)*x(3)+0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      par.k_ct*par.a_te*(x(3)-x(4));
exp4 = par.n*(par.alpha*x(1)*x(3)+0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      par.k_channel*par.a_te*(x(4)-x(5))/par.delta_channel;
ls5 = fcn2optimexpr(@(x)Heat_transfer_coff_SCH4(x(5),x(6),i,par)*par.a_g*(x(5)-x(6)),x);
exp5 = par.n*(par.alpha*x(1)*x(3)+0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      ls5;
ls6 = fcn2optimexpr(@(x)par.mg*(refpropm('H','T',x(6),'P',par.P,'methane')-refpropm('H','T',T_g_in,'P',par.P,'methane')),x);
exp6 = ls6-ls5;
eqn1 = exp1 == 0;
eqn2 = exp2 == 0;
eqn3 = exp3 == 0;
eqn4 = exp4 == 0;
eqn5 = exp5 == 0;
eqn6 = exp6 == 0;
a = 1;
b = 1;
%注意，这里采用全局优化求解器要给约束1加上乘项，不然求解器识别约束会报错
limit1 = (x(3)-x(2))>=0 ;
limit2 = par.n*(par.alpha*x(1)*(x(3)-x(2))+par.R*x(1)^2) >= 0;
%增加梯度项
% target = a*(x(7)-T_chip_taget)^2;
if i == 1
    target = a*(x(7)-T_chip_taget)^2;
else
    %把连续性要求也放到里面去
    target = a*(x(7)-T_chip_taget)^2+b*(I_record(i-1,j)-x(1))^2;
end
%设定问题类型，是taget对象的最小化问题
prob = optimproblem("Objective",target);
% prob.Objective.temp = target;
% prob.Objective.power = par.n*(par.alpha*x(1)*(x(3)-x(2))+par.R*x(1)^2);
%等式约束
prob.Constraints.eqn1 = eqn1;
prob.Constraints.eqn2 = eqn2;
prob.Constraints.eqn3 = eqn3;
prob.Constraints.eqn4 = eqn4;
prob.Constraints.eqn5 = eqn5;
prob.Constraints.eqn6 = eqn6;
% prob.Constraints.limit1 = limit1;
%不等式约束
prob.Constraints.limit2 = limit2;

%查看可用求解器
[~,validsolvers] = solvers(prob);
% disp(validsolvers);
% cansave = x(3)-x(2) >= 0;
% 
% prob.Constraints.cansave = cansave;
x0.x = (ub-lb).*rand(1,7)+lb;
% x0.x = [230,185,198,202,243,300];
Options = optimoptions("lsqnonlin",'Display','none','MaxFunctionEvaluations',3e4,'ConstraintTolerance',1e-5);
% Options = optimoptions("gamultiobj",'Display','iter');
[sol,fval,exitflag] = solve(prob,x0,'Options',Options);
% exitflag
flag = exitflag;

I = sol.x(1);
Tc = sol.x(2);
Th = sol.x(3);
To = sol.x(4);
Ti = sol.x(5);
Tg = sol.x(6);
T_chip = sol.x(7);
% while flag < -1
%     x0.x = [T_chip,Tc,Th,To,Ti,Tg];
%     Options = optimoptions('lsqnonlin','MaxFunctionEvaluations',3e4,'MaxIterations',3e4,'Algorithm','interior-point','OptimalityTolerance',1e-8);
%     [sol,fval,exitflag] = solve(prob,x0,'Options',Options);
% 
%     T_chip = sol.x(1);
%     Tc = sol.x(2);
%     Th = sol.x(3);
%     To = sol.x(4);
%     Ti = sol.x(5);
%     Tg = sol.x(6);
%     flag = exitflag
% end
% T_chip = sol.x(1);
% Tc = sol.x(2);
% Th = sol.x(3);
% To = sol.x(4);
% Ti = sol.x(5);
% Tg = sol.x(6);
% flag = exitflag;
[h,n,choice] = Heat_transfer_coff_SCH4(Ti,Tg,i,par);
% i

end