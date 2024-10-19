%% 本脚本是测试大规模稀疏非线性方程组求解的主程序
%% 设定基础参数
par = para_fixed;
heat_mode = 4;
h = par.h;
slope = par.slope;
num_unit = 2;
n = 6*num_unit;
I = zeros(num_unit,1);
%% 构建雅各比稀疏模式
e = ones(n,1);
Jstr = spdiags([e,e,e,e,e,e,e,e,e,e,e],-5:5,n,n);
% spy(Jstr);
%% 进行优化设置
%顺序Tg,Ti,To,Th,Tc,Tchip
lb = [100,100,50,50,50,50];
lb = repmat(lb,num_unit,1);
ub = [600,930,800,800,800,800];
ub = repmat(ub,num_unit,1);
% options = optimoptions('fsolve','Algorithm','trust-region','JacobPattern',Jstr);

options = optimoptions("lsqnonlin","Algorithm","trust-region-reflective",'JacobPattern',Jstr);
x0 = [300,243,202,198,195,210]
x0  = repmat(x0,num_unit,1);

[x,resnorm,residual,exitflag,output] = lsqnonlin(@(x)function_solver_huge_sparse_matrix(x,I,heat_mode,h,slope,par),x0,lb,ub);





