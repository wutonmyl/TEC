function [F] = function_solver_huge_sparse_matrix(x,I,heat_mode,h,slope,par)
%UNTITLED 本函数适用于求解TEC问题的大型稀疏矩阵非线性方程组表征
%   输入：给定的电流，基础参数,方正组未知量（包含方程组个数）
%   输出：对289.*6个非线性方程组的求解结果

n = length(x);
x_1 = 1:6:n;
%Ti
x_2 = 2:6:n;
%To
x_3 = 3:6:n;
%Th
x_4 = 4:6:n;
%Tc
x_5 = 5:6:n;
%Tchip
x_6 = 6:6:n;
%Tg
x_7 = 7:6:n;
%Ti
x_8 = 8:6:n;
F = zeros(n,1);
F(1,1) = par.mg.*(refpropm('H','T',x(1),'P',par.P,'methane')-refpropm('H','T',par.T_g_in,'P',par.P,'methane'))-...
    Heat_transfer_coff_SCH4(x(2),x(1),x_6/6,par).*par.a_g.*(x(2)-x(1));
F(x_2,1) = Heat_transfer_coff_SCH4(x(x_2),x(x_1),x_6/6,par).*par.a_g.*(x(x_2)-x(x_1))-...
    par.n.*(par.alpha.*I(x_6/6).*x(x_4)+0.5.*I(x_6/6).^2.*par.R+par.k_p.*par.a_copper.*(x(x_5)-x(x_4))/par.delta_p);
F(x_3,1) = par.n.*(par.alpha.*I(x_6/6).*x(x_4)+0.5.*I(x_6/6).^2.*par.R+par.k_p.*par.a_copper.*(x(x_5)-x(x_4))/par.delta_p)-...
    par.k_channel.*par.a_te.*(x(x_3)-x(x_2))/par.delta_channel;
F(x_4,1) = par.n.*(par.alpha.*I(x_6/6).*x(x_4)+0.5.*I(x_6/6).^2.*par.R+par.k_p.*par.a_copper.*(x(x_5)-x(x_4))/par.delta_p)-...
    par.k_ct.*par.a_te.*(x(x_4)-x(x_3));
F(x_5,1) = par.n.*(par.alpha.*I(x_6/6).*x(x_5)-0.5.*I(x_6/6).^2.*par.R+par.k_p.*par.a_copper.*(x(x_5)-x(x_4))/par.delta_p)-...
    par.k_ct.*par.a_te.*(x(x_6)-x(x_5));
F(x_6,1) = heat_change(x_6/6,heat_mode,h,slope)-par.k_ct.*par.a_te.*(x(x_6)-x(x_5));
F(x_7,1) = par.mg.*(refpropm('H','T',x(x_7),'P',par.P,'methane')-refpropm('H','T',x(x_1),'P',par.P,'methane'))-...
    Heat_transfer_coff_SCH4(x(x_8),x(x_7),x_6/6+1,par).*par.a_g.*(x(x_8)-x(x_7));
end