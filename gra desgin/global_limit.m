function [c,ceq] = global_limit(x,par,T_g_in,i,j,T_chip,heat_change,Heat_transfer_coff_SCH4)
%用于全局求解器中的约束函数
%


c = [];
ceq(1) = heat_change()-par.k_ct*par.a_te*(T_chip-x(2));
ceq(2) = par.n*(par.alpha*x(1)*x(2)-0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      par.k_ct*par.a_te*(T_chip-x(2));
ceq(3) = par.n*(par.alpha*x(1)*x(3)+0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      par.k_ct*par.a_te*(x(3)-x(4));
ceq(4) = par.n*(par.alpha*x(1)*x(3)+0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      par.k_channel*par.a_te*(x(4)-x(5))/par.delta_channel;
ls5 = fcn2optimexpr(@(x)Heat_transfer_coff_SCH4(x(5),x(6),i,par)*par.a_g*(x(5)-x(6)),x);
ceq(5) = par.n*(par.alpha*x(1)*x(3)+0.5*x(1)^2*par.R+par.k_p*par.a_copper*(x(2)-x(3))/par.delta_p)-...
      ls5;
ls6 = fcn2optimexpr(@(x)par.mg*(refpropm('H','T',x(6),'P',par.P,'methane')-refpropm('H','T',T_g_in,'P',par.P,'methane')),x);
ceq(6) = ls6-ls5;
 
end