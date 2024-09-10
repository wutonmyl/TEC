addpath('current_change\');
addpath("basic_fuction\");
par = para;

%构造数据记录
i_record = zeros(1,289);
Dz_record = zeros(1,289);
Tchip_record = zeros(1,289);
Tc_record   = zeros(1,289);
Th_record   = zeros(1,289);
To_record   = zeros(1,289);
Ti_record   = zeros(1,289);
Tg_record   = zeros(1,289);
h_record    = zeros(1,289);
q1_r_record = zeros(1,289);
q2_l_record = zeros(1,289);
q2_r_record = zeros(1,289);
q3_l_record = zeros(1,289);
q3_r_record = zeros(1,289);
q4_r_record = zeros(1,289);
q5_r_record = zeros(1,289);
flag_record = zeros(1,289);
n_record = zeros(1,289);
choice_record = zeros(1,289);
mean_T_chip = zeros(1,289);
median_T_chip = zeros(1,289);
geomean_T_chip = zeros(1,289);
harmmean_T_chip = zeros(1,289);
range_T_chip = zeros(1,289);
var_T_chip = zeros(1,289);
current_self_record = zeros(1,289);
Q_tec_record = zeros(1,289);
COP_record = zeros(1,289);
exitflag_record = zeros(1,289);
output_record = cell(1,289);

dz = par.Dz;
Dz = 0.003;
i = 1;
T_g_in = par.T_g_in;
j = 1;
current_self = 0;
diary('myDiaryFile');
while Dz<=par.length
    T_chip_target = 300;
    % [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = current_influ_function_solver(par,T_g_in,i,j);
    
   
    %优化电流i使得温度调平
    %建立优化问题
    prob = optimproblem('ObjectiveSense','minimize');
    %给出优化变量电流
    x = optimvar('x',1,'LowerBound',0,'UpperBound',5);
    %规范化要用到的函数
    [T_chip_tem,Tc_tem,Th_tem,To_tem,Ti_tem,Tg_tem,h_tem,flag_tem,n_tem,choice_tem] = fcn2optimexpr(@current_influ_function_solver,par,T_g_in,i,x);
    %建立优化表达式（这里为了防止收敛于不可行点）
    Const = fcn2optimexpr(@flag_shift,flag_tem);
    % prob = optimproblem("Objective",(T_chip_tem+Const-T_chip_target)^2);
    prob.Objective = (T_chip_tem+Const-T_chip_target)^2;
    %给出优化约束
    cansolve = flag_tem>=1;
    prob.Constraints.cansolve = cansolve;
    %//全局最小值方法
    %优化选项
    options = optimoptions(prob);
    options.Display ="iter";
    x0.x = 5*rand();
    [sol,fval,exitflag,output] = solve(prob,x0,"Solver",'patternsearch','Options',options);
    %//局部最小值方法
    %  %优化选项
    % options = optimoptions(prob);
    % options.Display ="iter";
    % options.MaxIterations = 800;
    % options.Algorithm = "interior-point";
    
    % %给出优化初值
    % exitflag = -1;
    % num = 0;
    % while exitflag < 0
    %     x0.x = rand() ;
    % 
    %     % options.OutputFcn = @LMoutfun;
    %     [sol,fval,exitflag] = solve(prob,x0,"Options",options);
    %     num = num+1;
    %     if num>=2
    %          disp('处于调制初值阶段');
    % 
    %     end
    %     if num>=5
    %         break
    %     end
    % 
    % end
    current_self = sol.x(1);
    [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = current_influ_function_solver(par,T_g_in,i,current_self);
    
    %计算并记录对应结果
    q2_r = par.k_ct*par.a_te*(T_chip-Tc);
    q2_l = par.n*(par.alpha*current_change(current_self)*Tc-0.5*current_change(current_self)^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
    q3_r = par.k_ct*par.a_te*(Th-To);
    q3_l = par.n*(par.alpha*current_change(current_self)*Th+0.5*current_change(current_self)^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
    q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
    q5_r= h*par.a_g*(Ti-Tg);
    Q_tec = q3_l-q2_l;
    COP = q2_l/Q_tec;
    
    i_record(j,i) = i;
    Tchip_record(j,i) = T_chip;
    Tc_record(j,i) = Tc;
    Th_record(j,i) = Th;
    To_record(j,i) = To;
    Ti_record(j,i) = Ti ;
    Tg_record(j,i) = Tg;
    h_record(j,i) = h;
    
    q2_l_record(j,i) = q2_l;
    q2_r_record(j,i) = q2_r;
    q3_l_record(j,i) = q3_l;
    q3_r_record(j,i) = q3_r;
    q4_r_record(j,i) = q4_r;
    q5_r_record(j,i) = q5_r;
    flag_record(j,i) = flag;
    n_record(j,i) = n;
    choice_record(j,i) = choice;
    current_self_record(j,i) = current_self;
    COP_record(j,i) = COP;
    Q_tec_record(j,i) = Q_tec;
    exitflag_record(j,i) = exitflag;
    Dz_record(j,i) = Dz;
    output_record{j,i} = output;
    Dz = Dz+dz;
    
    i = i+1;
    i;
    disp(i);
    rate = Dz/par.length;
    disp(rate);
    
   
    T_g_in = Tg;
    disp(exitflag);
    
    


end
diary("off");





