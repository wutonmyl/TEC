%独立调控的主程序
%放宽优化要求，平时不开TEC，当温度发展或者梯度发展不满足要求时使用TEC
diary off
%引入热量变化函数库，确保可以计算不同热量分布的工况
addpath(['heat_change\'])
%引入基本方程，确保可以计算对流换热系数
addpath('basic_fuction\')
%引入参数类文件夹
addpath('parament\')
%引入求解器文件夹
addpath("control_function_solver\")

%导入零电流工况基础的沿程功率器件温度数据作为控制的上下限依据
Tchip_zero_record = load('material_0_SCH4_v_0.007_Tin_185K_初值随机_I_0_uniform.mat','Tchip_record');

Tchip_zero_ub = max(Tchip_zero_record.Tchip_record);
Tchip_zero_lb = min(Tchip_zero_record.Tchip_record);

par = para;
search_num = 1;
mode = 1;
h = par.h;
slope = par.slope;
% %设置迭代环节基础参数
% dz = par.Dz;
% Dz = dz;
% i = 1;
% j = 1 ;
% T_g_in = par.T_g_in;
% search_num = 1;


%构造数据记录
diary('test_DiaryFile');
% Q_tec_record = zeros(289,search_num);
% COP_record = zeros(289,search_num);
% i_record = zeros(289,search_num);
% Dz_record = zeros(289,search_num);
% Tchip_record = zeros(289,search_num);    
% Tc_record   = zeros(289,search_num);
% Th_record   = zeros(289,search_num);
% To_record   = zeros(289,search_num);
% Ti_record   = zeros(289,search_num);
% Tg_record   = zeros(289,search_num);
% h_record    = zeros(289,search_num);
% q1_r_record = zeros(289,search_num);
% q2_l_record = zeros(289,search_num);
% q2_r_record = zeros(289,search_num);
% q3_l_record = zeros(289,search_num);
% q3_r_record = zeros(289,search_num);
% q4_r_record = zeros(289,search_num);
% q5_r_record = zeros(289,search_num);
% flag_record = zeros(289,search_num);
% n_record = zeros(289,search_num);
% choice_record = zeros(289,search_num);
% I_record = zeros(289,search_num);
% grad_record = zeros(search_num);
% range_record = zeros(search_num);
% Cp_g_record = zeros(search_num);
% density_g_record = zeros(search_num);
% H_g_record = zeros(search_num);
% lamda_g_record = zeros(search_num);
% Pr_g_record = zeros(search_num);
% mu_g_record = zeros(search_num);
% T_chip_target_record = zeros(search_num);
% grad = zeros(289,search_num);
% need_control_record = zeros(289,search_num);
% i_special_record = zeros(289,search_num);

% i_special =1;

%% 构建优化问题（存疑，怎么大的嵌套小的问题还是两说）
% 优化目标：梯度限制最小
% 优化约束：能跑完整个沿程流程，中间不出现梯度过低，函数陷入循环的情况
% 难点：约束函数怎么构建，需要将我们的整个main函数重构一遍

%定义随机种子
rng shuffle

%定义优化变量
grad_ub = 56.83*0.3 ;
grad_lb = 0;
x = optimvar('x',1,'LowerBound',grad_lb,'UpperBound',grad_ub);

%给出优化约束，这里需要原本main5的主程序部分，程序命名为
% data = current_control_set5(par,Tchip_zero_ub,Tchip_zero_lb,x,search_num,mode,h,slope);
[data,wrong_num_record,i_final] = fcn2optimexpr(@(x)current_control_set5(par,Tchip_zero_ub,Tchip_zero_lb,x,search_num,mode,h,slope),x,'ReuseEvaluation',true);


limit1 = wrong_num_record == zeros(289,search_num);
limit2 = i_final == 290;
%给出优化目标
target = x;
%构建问题
prob = optimproblem('Objective',target);
% prob.Constraints.limit1 = limit1;
prob.Constraints.limit2 = limit2;

%查看可用的求解器
[~,validsolvers] = solvers(prob);
disp(validsolvers);
x0.x = (Tchip_zero_ub-Tchip_zero_lb).*rand()+Tchip_zero_lb;
disp(x0.x);

%创建优化选型
Options = optimoptions("fmincon",'Display','iter','ConstraintTolerance',1e-8);
[sol,fval,exitflag] = solve(prob,x0,'Options',Options,'Solver','fmincon');
% exitflag
grad_flag = exitflag;
grad_bond = sol.x;

[data,~,~] = current_control_set5(par,Tchip_zero_ub,Tchip_zero_lb,grad_bond,search_num,mode,h,slope);


%% 开始沿程循环
% while Dz<=par.length
    
%% 默认用无电流模式开始工作
%     [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = function_solver(par,T_g_in,i,mode,h,slope);
%     %先开始记录数据
%     q2_r = par.k_ct*par.a_te*(T_chip-Tc);
%     q2_l = par.n*(par.alpha*par.I*Tc-0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
%     q3_r = par.k_ct*par.a_te*(Th-To);
%     q3_l = par.n*(par.alpha*par.I*Th+0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
%     q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
%     q5_r= h*par.a_g*(Ti-Tg);
%     Q_tec = q3_l-q2_l;
%     Q_tec = abs(Q_tec);
%     % COP = 1/(Q_tec/q2_l);
%     [Cp_g,density_g,H_g,mu_g,lamda_g,Pr_g] = refpropm('CDHVL^','T',Tg,'P',par.P,'methane');
% 
% 
%     Cp_g_record(i,j) = Cp_g;
%     density_g_record(i,j) = density_g;
%     H_g_record(i,j) = H_g;
%     mu_g_record(i,j) = mu_g;
%     lamda_g_record(i,j) = lamda_g;
%     Pr_g_record(i,j) = Pr_g;
%     i_record(i,j) = i;
%     Tchip_record(i,j) = T_chip;
%     I_record(i,j) = 0;
%     Tc_record(i,j) = Tc;
%     Th_record(i,j) = Th;
%     To_record(i,j) = To;
%     Ti_record(i,j) = Ti;
%     Tg_record(i,j) = Tg;
%     h_record(i,j) = h;
% 
%     q2_l_record(i,j) = q2_l;
%     q2_r_record(i,j) = q2_r;
%     q3_l_record(i,j) = q3_l;
%     q3_r_record(i,j) = q3_r;
%     q4_r_record(i,j) = q4_r;
%     q5_r_record(i,j) = q5_r;
%     flag_record(i,j) = flag;
%     n_record(i,j) = n;
%     choice_record(i,j) = choice;
%     COP_record(i,j) = 0;
%     Q_tec_record(i,j) = Q_tec;
%     Dz_record(i,j) = Dz;
% 
%     %% 开始判断是否满足控温标准
% 
%     [need_control,Tchip_target_recmd] = judge_need_control(Tchip_record,i,j,Tchip_zero_lb,Tchip_zero_ub);
%     need_control_record(i,j) = need_control;
%     % % %存储突变前的点
%     if need_control == 0
%         i_special = i;
%     end
%     i_special_record(i,j) = i_special;
%     T_chip_target_record(i,j) = Tchip_target_recmd;
%     %命令行显示进度
%     disp(i);
%     disp(j);
%     disp('推荐温度为');
%     disp(Tchip_target_recmd)
%     disp('实际温度为');
%     disp(T_chip);
%     %% 如果需要控温，采用独立调控程序
%     contol_repeat = 0;
%     wrong_num = 0;
%     while need_control ~= 0
%         disp('进入独立调控环节')
%         contol_repeat = contol_repeat + 1;
%         disp('当前重复调控次数')
%         disp(contol_repeat);
%         [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed5(par,T_g_in,i,j,Tchip_target_recmd,I_record,mode,h,slope);
%         disp(i);
%         disp(j);
%         disp('推荐温度为');
%         disp(Tchip_target_recmd)
%         disp('实际温度为');
%         disp(T_chip);
%         disp('不满足记号为');
%         disp(need_control);
%         initial_change = 0;
% 
% 
%         %求解器抱错误解进入循环
%         while flag <= 0 
% 
%             %如果控不住，换随机初始值再来一次
%             % T_chip_target = T_chip_target + 1;
%             initial_change = initial_change + 1;
%             [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed5(par,T_g_in,i,j,Tchip_target_recmd,I_record,mode,h,slope);
%             disp(i);
%             disp(j);
%             disp('更换随机初始值数');
%             disp(initial_change);
%             disp('实际温度为');
%             disp(T_chip);
%         end
%         %开始数据记录
% 
% 
%         % %如果一直循环，比如循环超过5次
%         % if initial_change >= 5
%         % 
% 
% 
% 
%         q2_r = par.k_ct*par.a_te*(T_chip-Tc);
%         q2_l = par.n*(par.alpha*I*Tc-0.5*I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
%         q3_r = par.k_ct*par.a_te*(Th-To);
%         q3_l = par.n*(par.alpha*I*Th+0.5*I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
%         q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
%         q5_r= h*par.a_g*(Ti-Tg);
%         Q_tec = q3_l-q2_l;
%         Q_tec = abs(Q_tec);
%         COP = 1/(Q_tec/q2_l);
%         [Cp_g,density_g,H_g,mu_g,lamda_g,Pr_g] = refpropm('CDHVL^','T',Tg,'P',par.P,'methane');
% 
% 
%         Cp_g_record(i,j) = Cp_g;
%         density_g_record(i,j) = density_g;
%         H_g_record(i,j) = H_g;
%         mu_g_record(i,j) = mu_g;
%         lamda_g_record(i,j) = lamda_g;
%         Pr_g_record(i,j) = Pr_g;
%         i_record(i,j) = i;
%         Tchip_record(i,j) = T_chip;
%         I_record(i,j) = I;
%         Tc_record(i,j) = Tc;
%         Th_record(i,j) = Th;
%         To_record(i,j) = To;
%         Ti_record(i,j) = Ti;
%         Tg_record(i,j) = Tg;
%         h_record(i,j) = h;
% 
%         q2_l_record(i,j) = q2_l;
%         q2_r_record(i,j) = q2_r;
%         q3_l_record(i,j) = q3_l;
%         q3_r_record(i,j) = q3_r;
%         q4_r_record(i,j) = q4_r;
%         q5_r_record(i,j) = q5_r;
%         flag_record(i,j) = flag;
%         n_record(i,j) = n;
%         choice_record(i,j) = choice;
%         COP_record(i,j) = COP;
%         Q_tec_record(i,j) = Q_tec;
%         Dz_record(i,j) = Dz;
%         %结束数据记录
% 
% 
%         %再检查一遍是否需要调控
%         [need_control,Tchip_target_recmd] = judge_need_control(Tchip_record,i,j,Tchip_zero_lb,Tchip_zero_ub);
%         T_chip_target_record(i,j) = Tchip_target_recmd;
%     end
%     %结束一个单元的计算
%     i = i+1;
%     rate = Dz/par.length;
%     Dz = Dz+dz;
% 
%     T_g_in = Tg;
% 
% end
%% 数据处理
% grad_T_chip = gradient(Tchip_record(:,j));
% mean_T_chip = mean(Tchip_record(:,j));
% median_T_chip = median(Tchip_record(:,j));
% geomean_T_chip = geomean(Tchip_record(:,j));
% harmmean_T_chip = harmmean(Tchip_record(:,j));
% range_T_chip = range(Tchip_record(:,j));
% var_T_chip = var(Tchip_record(:,j));
% std_T_chip = std(Tchip_record(:,j));
diary off;
% show
