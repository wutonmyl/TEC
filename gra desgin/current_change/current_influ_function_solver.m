function [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = current_influ_function_solver(par,T_g_in,i,j)
%FUNCTION_SOLVER 此处显示有关此函数的摘要
%   此处显示详细说明
rng default
x = optimvar('x',6,'LowerBound',[50,50,50,50,100,100],'UpperBound',[800,800,800,800,930,600]);
%第一个式子关于pchip还需要修改
eq1 = par.P_chip==par.k_ct*par.a_te*(x(1)-x(2));
eq2 = par.n*(par.alpha*current_change(j)*x(2)-0.5*current_change(j)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)==...
      par.k_ct*par.a_te*(x(1)-x(2));
eq3 = par.n*(par.alpha*current_change(j)*x(3)+0.5*current_change(j)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)==...
      par.k_ct*par.a_te*(x(3)-x(4));
eq4 = par.n*(par.alpha*current_change(j)*x(3)+0.5*current_change(j)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)==...
      par.k_channel*par.a_te*(x(4)-x(5))/par.delta_channel;
ls5 = fcn2optimexpr(@(x)Heat_transfer_coff_SCH4(x(5),x(6),i,par)*par.a_g*(x(5)-x(6)),x);
eq5 = par.n*(par.alpha*current_change(j)*x(3)+0.5*current_change(j)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)==...
      ls5;
ls6 = fcn2optimexpr(@(x)par.mg*(refpropm('H','T',x(6),'P',par.P,'methane')-refpropm('H','T',T_g_in,'P',par.P,'methane')),x);
eq6 = ls6==ls5;

prob = eqnproblem;
prob.Equations.eq1 = eq1;
prob.Equations.eq2 = eq2;
prob.Equations.eq3 = eq3;
prob.Equations.eq4 = eq4;
prob.Equations.eq5 = eq5;
prob.Equations.eq6 = eq6;
x0.x = [230,185,198,202,243,300];
% x0.x = [230,185,198,202,243,300];
Options = optimoptions('lsqnonlin','Display','none','MaxFunctionEvaluations',3e4,'MaxIterations',3e4,'Algorithm','levenberg-marquardt',OptimalityTolerance=1e-8);
[sol,fval,exitflag] = solve(prob,x0,'Options',Options);
% exitflag
flag = exitflag;

T_chip = sol.x(1);
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
