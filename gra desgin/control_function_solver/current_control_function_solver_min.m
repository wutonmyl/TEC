function [I,Tc,Th,To,Ti,Tg,h,flag,n,choice] = current_control_function_solver_min(par,T_g_in,i,j,T_chip)
%   这里要把芯片温度变成电流
%   此处显示详细说明
rng default
% gs = GlobalSearch;
x = optimvar('x',6,'LowerBound',[-6,100,100,100,100,100],'UpperBound',[6,800,800,800,930,600]);
%第一个式子关于pchip还需要修改
exp1 = heat_change(i)-par.k_ct*par.a_te*(T_chip-x(2));
exp2 = par.n*(par.alpha*x(1)*x(2)-0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      par.k_ct*par.a_te*(T_chip-x(2));
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
% limit1 = x(1)*(x(3)-x(2)) >= 0;
Power = (par.n*(par.alpha*x(1)*(x(3)-x(2))+par.R*x(1)^2))^2;

prob = optimproblem('Objective',Power);
prob.Constraints.eqn1 = eqn1;
prob.Constraints.eqn2 = eqn2;
prob.Constraints.eqn3 = eqn3;
prob.Constraints.eqn4 = eqn4;
prob.Constraints.eqn5 = eqn5;
prob.Constraints.eqn6 = eqn6;
% prob.Constraints.limit1 = limit1;
% cansave = x(3)-x(2) >= 0;
% 
% prob.Constraints.cansave = cansave;
x0.x = [3,185,198,202,243,300];
% x0.x = [230,185,198,202,243,300];
Options = optimoptions('lsqnonlin','Display','final','Algorithm','interior-point');

[sol,fval,exitflag] = solve(prob,x0,'Options',Options);
% exitflag
flag = exitflag;

I = sol.x(1);
Tc = sol.x(2);
Th = sol.x(3);
To = sol.x(4);
Ti = sol.x(5);
Tg = sol.x(6);
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