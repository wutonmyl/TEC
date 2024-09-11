%独立调控的主程序
%放宽优化要求，平时不开TEC，当温度发展或者梯度发展不满足要求时使用TEC

%引入热量变化函数库，确保可以计算不同热量分布的工况
addpath(['heat_change\'])
%引入基本方程，确保可以计算对流换热系数
addpath('basic_fuction\')
%引入参数类文件夹
addpath('parament\')
%引入求解器文件夹
addpath("control_function_solver\")
par = para;
dz = par.Dz;
Dz = dz;
i = 1;
T_g_in = par.T_g_in;
search_num = 1;
%构造数据记录
% diary('myDiaryFile');
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
T_chip_target_record = zeros(search_num);
grad = zeros(289,search_num);
need_control_record = zeros(289,search_num);
j = 1 ;
i_special =1;
%开始沿程循环
while Dz<=par.length
    %默认用无电流模式开始工作
    [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = function_solver(par,T_g_in,i);
    %先开始记录数据
    q2_r = par.k_ct*par.a_te*(T_chip-Tc);
    q2_l = par.n*(par.alpha*par.I*Tc-0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
    q3_r = par.k_ct*par.a_te*(Th-To);
    q3_l = par.n*(par.alpha*par.I*Th+0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
    q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
    q5_r= h*par.a_g*(Ti-Tg);
    Q_tec = q3_l-q2_l;
    Q_tec = abs(Q_tec);
    % COP = 1/(Q_tec/q2_l);
    [Cp_g,density_g,H_g,mu_g,lamda_g,Pr_g] = refpropm('CDHVL^','T',Tg,'P',par.P,'methane');
    
    
    Cp_g_record(i,j) = Cp_g;
    density_g_record(i,j) = density_g;
    H_g_record(i,j) = H_g;
    mu_g_record(i,j) = mu_g;
    lamda_g_record(i,j) = lamda_g;
    Pr_g_record(i,j) = Pr_g;
    i_record(i,j) = i;
    Tchip_record(i,j) = T_chip;
    I_record(i,j) = 0;
    Tc_record(i,j) = Tc;
    Th_record(i,j) = Th;
    To_record(i,j) = To;
    Ti_record(i,j) = Ti;
    Tg_record(i,j) = Tg;
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
    COP_record(i,j) = 0;
    Q_tec_record(i,j) = Q_tec;
    Dz_record(i,j) = Dz;

    %开始判断是否满足控温标准
    
    [need_control,Tchip_target_recmd] = judge_need_control(Tchip_record,i,j,i_special);
    need_control_record(i,j) = need_control;
    %存储突变前的点
    if need_control == 0
        i_special = i;
    end
    T_chip_target_record(i,j) = Tchip_target_recmd;
    %命令行显示进度
    disp(i);
    disp(j);
    disp('推荐温度为');
    disp(Tchip_target_recmd)
    disp('实际温度为');
    disp(T_chip);
    %如果需要控温，采用独立调控程序
    contol_repeat = 0;
    while need_control ~= 0
        disp('进入独立调控环节')
        contol_repeat = contol_repeat + 1;
        disp('当前重复调控次数')
        disp(contol_repeat);
        [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed5(par,T_g_in,i,j,Tchip_target_recmd,I_record);
        disp(i);
        disp(j);
        disp('推荐温度为');
        disp(Tchip_target_recmd)
        disp('实际温度为');
        disp(T_chip);
        disp('不满足记号为');
        disp(need_control);
        initial_change = 0;
        while flag <= 0 
    
            %如果控不住，换随机初始值再来一次
            % T_chip_target = T_chip_target + 1;
            initial_change = initial_change + 1;
            [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed5(par,T_g_in,i,j,Tchip_target_recmd,I_record);
            disp(i);
            disp(j);
            disp('更换随机初始值数');
            disp(initial_change);
            disp('实际温度为');
            disp(T_chip);
        end
        %开始数据记录
        q2_r = par.k_ct*par.a_te*(T_chip-Tc);
        q2_l = par.n*(par.alpha*I*Tc-0.5*I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
        q3_r = par.k_ct*par.a_te*(Th-To);
        q3_l = par.n*(par.alpha*I*Th+0.5*I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
        q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
        q5_r= h*par.a_g*(Ti-Tg);
        Q_tec = q3_l-q2_l;
        Q_tec = abs(Q_tec);
        COP = 1/(Q_tec/q2_l);
        [Cp_g,density_g,H_g,mu_g,lamda_g,Pr_g] = refpropm('CDHVL^','T',Tg,'P',par.P,'methane');
        
        
        Cp_g_record(i,j) = Cp_g;
        density_g_record(i,j) = density_g;
        H_g_record(i,j) = H_g;
        mu_g_record(i,j) = mu_g;
        lamda_g_record(i,j) = lamda_g;
        Pr_g_record(i,j) = Pr_g;
        i_record(i,j) = i;
        Tchip_record(i,j) = T_chip;
        I_record(i,j) = I;
        Tc_record(i,j) = Tc;
        Th_record(i,j) = Th;
        To_record(i,j) = To;
        Ti_record(i,j) = Ti;
        Tg_record(i,j) = Tg;
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
        %再检查一遍是否需要调控
        [need_control,Tchip_target_recmd] = judge_need_control(Tchip_record,i,j,i_special);
        T_chip_target_record(i,j) = Tchip_target_recmd;
    end
    %结束一个单元的计算
    i = i+1;
    rate = Dz/par.length;
    Dz = Dz+dz;
    
    T_g_in = Tg;

end
grad_T_chip = gradient(Tchip_record(:,j));
mean_T_chip = mean(Tchip_record(:,j));
median_T_chip = median(Tchip_record(:,j));
geomean_T_chip = geomean(Tchip_record(:,j));
harmmean_T_chip = harmmean(Tchip_record(:,j));
range_T_chip = range(Tchip_record(:,j));
var_T_chip = var(Tchip_record(:,j));
std_T_chip = std(Tchip_record(:,j));


