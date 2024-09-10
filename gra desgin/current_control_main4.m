%独立调控的主程序
%放宽优化要求，以单一单元作为优化对象

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
search_num = 2;
%构造数据记录
% diary('myDiaryFile');
Q_tec_record = zeros(289,search_num);
COP_record = zeros(289,search_num);
i_record = zeros(289,search_num);
Dz_record = zeros(289,search_num);

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
j = 1 ;

if Tchip_record(289,1) == 0
Tchip_record = zeros(289,search_num);    
    while j == 1
        Dz = dz;
        i = 1;
        T_g_in = par.T_g_in;
        %先跑一遍零电流工况
        while Dz<=par.length
            [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = function_solver(par,T_g_in,i);
       
        
            q2_r = par.k_ct*par.a_te*(T_chip-Tc);
            q2_l = par.n*(par.alpha*par.I*Tc-0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
            q3_r = par.k_ct*par.a_te*(Th-To);
            q3_l = par.n*(par.alpha*par.I*Th+0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
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
            COP_record(i,j) = COP;
            Q_tec_record(i,j) = Q_tec;
            Dz_record(i,j) = Dz;
        
            i = i+1;
            rate = Dz/par.length;
            disp(rate)
            Dz = Dz+dz;
            
           
            T_g_in = Tg;
        end
    
        %计算零电流情况的统计数据
        
        mean_T_chip_zero = mean(Tchip_record(:,j));
        median_T_chip_zero = median(Tchip_record(:,j));
        geomean_T_chip_zero = geomean(Tchip_record(:,j));
        harmmean_T_chip_zero = harmmean(Tchip_record(:,j));
        range_T_chip_zero = range(Tchip_record(:,j));
        var_T_chip_zero = var(Tchip_record(:,j));
        std_T_chip_zero = std(Tchip_record(:,j));
        %进入下个工况环节
        j = j + 1;
    end
else
    j = j+1;
    disp('已有零电流工况数据，跳过一段计算环节')
end
%有了零电流的数据后开始优化

Dz = dz;
i = 1;
T_g_in = par.T_g_in;
%进入沿程循环 

while Dz<=par.length
    %采用方程求解函数4
    [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed4(par,T_g_in,i,j,Tchip_record,I_record);
    disp(i);
    disp(j);
    disp('实际温度为');
    disp(T_chip);
    initial_change = 0;
    while flag <= 0 

        %如果控不住，换随机初始值再来一次
        % T_chip_target = T_chip_target + 1;
        initial_change = initial_change + 1;
        [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed4(par,T_g_in,i,j,Tchip_record,I_record);
        disp(i);
        disp(j);
        disp('更换随机初始值数');
        disp(initial_change);
        disp('实际温度为');
        disp(T_chip);
    end
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
    % T_chip_target_record(i,j) = T_chip_target;

    i = i+1;
    rate = Dz/par.length;
    Dz = Dz+dz;
    
    T_g_in = Tg;
end
mean_T_chip = mean(Tchip_record(:,j));
median_T_chip = median(Tchip_record(:,j));
geomean_T_chip = geomean(Tchip_record(:,j));
harmmean_T_chip = harmmean(Tchip_record(:,j));
range_T_chip = range(Tchip_record(:,j));
var_T_chip = var(Tchip_record(:,j));
std_T_chip = std(Tchip_record(:,j));