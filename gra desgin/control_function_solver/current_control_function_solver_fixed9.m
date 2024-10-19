function [data,i_final,break_flag,Tchip_record,grad_T_chip,Q_tec_overall,COP_record] = current_control_function_solver_fixed9(par,Tchip_zero_lb,Tchip_zero_ub,mode,h,slope,I)
%% 本函数是用于计算main9的优化表达式
%   输入：给定的电流序列
%   输出：约束量-单元可行总长（不可行的话直接跳出循环结束计算了）；温度范围
%         优化量-梯度，总功率，COP
dz = par.Dz;
Dz = dz;
pressure = par.P;
mg = par.mg;
% p_chip = heat_change(i,mode,h,slope);

T_g_in = par.T_g_in;
j = 1;
% mode = 0;
% slope = par.slope;
% h = par.h;
% current_exist = 0;
search_num =1;

num_I = 289;
%构造数据记录
Q_tec_record = zeros(289,search_num);
COP_record = zeros(289,search_num);
i_record = zeros(289,search_num);
Dz_record = zeros(289,search_num);
Tchip_record = zeros(289,search_num);
Tc_record   = zeros(289,search_num);
Th_record   = zeros(289,search_num);
To_record   = zeros(289,search_num);
Ti_record   = zeros(289,search_num);
Tg_record   = zeros(289,search_num);
h_record    = zeros(289,search_num);
q1_r_record = zeros(289,search_num);
q2_l_record = zeros(289,search_num);
q2_r_record = zeros(289,search_num);
q3_l_record = zeros(289,search_num);
q3_r_record = zeros(289,search_num);
q4_r_record = zeros(289,search_num);
q5_r_record = zeros(289,search_num);
flag_record = zeros(289,search_num);
n_record = zeros(289,search_num);
choice_record = zeros(289,search_num);
I_record = zeros(289,search_num);
grad_record = zeros(search_num);
range_record = zeros(search_num);
Cp_g_record = zeros(search_num);
density_g_record = zeros(search_num);
H_g_record = zeros(search_num);
lamda_g_record = zeros(search_num);
Pr_g_record = zeros(search_num);
mu_g_record = zeros(search_num);
T_chip_target_record = zeros(289,search_num);
grad_flag = [];
grad_bond = [];
%% 选择电流模式

for i = 1:289
    [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = function_solver9(par,T_g_in,pressure,mg,i,mode,h,slope,I(i));
    disp(I(i))
    q2_r = par.k_ct*par.a_te*(T_chip-Tc);
    q2_l = par.n*(par.alpha*I(i)*Tc-0.5*I(i)^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
    q3_r = par.k_ct*par.a_te*(Th-To);
    q3_l = par.n*(par.alpha*I(i)*Th+0.5*I(i)^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
    q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
    q5_r= h*par.a_g*(Ti-Tg);
    Q_tec = q3_l-q2_l;
    % Q_tec = abs(Q_tec);
    COP = 1/(Q_tec/q2_l);
    [Cp_g,density_g,H_g,mu_g,lamda_g,Pr_g] = refpropm('CDHVL^','T',Tg,'P',pressure,'methane');
    
    
    Cp_g_record(i,j) = Cp_g;
    density_g_record(i,j) = density_g;
    H_g_record(i,j) = H_g;
    mu_g_record(i,j) = mu_g;
    lamda_g_record(i,j) = lamda_g;
    Pr_g_record(i,j) = Pr_g;
    i_record(i,j) = i;
    i_final = i;
    Tchip_record(i,j) = T_chip;
    I_record(i,j) = I(i);
    Tc_record(i,j) = Tc;
    Th_record(i,j) = Th;
    To_record(i,j) = To;
    Ti_record(i,j) = Ti;
    Tg_record(i,j) = Tg;
    %对流换热系数
    h_record(i,j) = h;
    
    q2_l_record(i,j) = q2_l;
    q2_r_record(i,j) = q2_r;
    q3_l_record(i,j) = q3_l;
    q3_r_record(i,j) = q3_r;
    q4_r_record(i,j) = q4_r;
    q5_r_record(i,j) = q5_r;
    flag_record(i,j) = flag;
    n_record(i,j) = n;
    choice_record(i,j) = choice;
    COP_record(i,j) = COP;
    Q_tec_record(i,j) = Q_tec;
    Dz_record(i,j) = Dz;
    disp(i);
    break_flag = 0;
    if flag<=0
        break_flag = 1
        break
    elseif COP < 0
        break_flag = -1
        break
    elseif T_chip<=Tchip_zero_lb||T_chip>=Tchip_zero_ub
        break_flag = -2
        break
    elseif T_chip>=Tchip_zero_ub
        break_flag = 2
        break
 

    end

    
    rate = Dz/par.length;
    disp(rate)
    
    
    T_g_in = Tg;
end

%构造统计数据

grad_T_chip = gradient(Tchip_record(:,j));
Q_tec_overall = sum(Q_tec_record(:,j));
mean_T_chip = mean(Tchip_record(:,j));
median_T_chip = median(Tchip_record(:,j));
geomean_T_chip = geomean(Tchip_record(:,j));
harmmean_T_chip = harmmean(Tchip_record(:,j));
range_T_chip = range(Tchip_record(:,j));
var_T_chip = var(Tchip_record(:,j));
std_T_chip = std(Tchip_record(:,j));

%记录数据
data.grad_T_chip = grad_T_chip;
data.Q_tec_overall = Q_tec_overall;
data.mean_T_chip = mean_T_chip;
data.median_T_chip = median_T_chip;
data.geomean_T_chip = geomean_T_chip;
data.harmmean_T_chip = harmmean_T_chip;
data.range_T_chip = range_T_chip;
data.var_T_chip = var_T_chip;
data.std_T_chip = std_T_chip;

data.Q_tec_record = Q_tec_record;
data.COP_record = COP_record;
data.i_record = i_record;
data.i_final = i_final;
data.Dz_record = Dz_record;
data.Tchip_record = Tchip_record;
data.Tc_record = Tc_record;
data.Th_record = Th_record;
data.To_record = To_record;
data.Ti_record = Ti_record;
data.Tg_record = Tg_record;
data.h_record = h_record;
data.q1_r_record = q1_r_record;
data.q2_r_record = q2_r_record;
data.q2_l_record = q2_l_record;
data.q3_r_record = q3_r_record;
data.q3_l_record = q3_l_record;
data.q4_r_record = q4_r_record;
data.q5_r_record = q5_r_record;
data.flag_record = flag_record;
data.n_record = n_record;
data.choice_record = choice_record;
data.I_record = I_record;


data.Cp_g_record = Cp_g_record;
data.density_g_record = density_g_record;
data.H_g_record = H_g_record;
data.lamda_g_record = lamda_g_record;
data.Pr_g_record = Pr_g_record;
data.mu_g_record = mu_g_record;
data.T_chip_target_record = T_chip_target_record;


