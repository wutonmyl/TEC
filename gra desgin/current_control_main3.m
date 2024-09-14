%独立调控的主程序,添加了两段蛇形通道（2的修改部分）
%将主程序调整为适宜fixed3的格式

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
%给出在控温区间内大致搜索的范围
% search_num = 35;
search_num = 1;
%构造数据记录
diary('myDiaryFile');
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

j = 1;
k = 1;
for j = 1:search_num
    Dz = dz;
    i = 1;
    T_g_in = par.T_g_in;
    %进入沿程循环
    while Dz<=par.length
        
        % if i>=2 && Tchip_record(i-1,j) ~= T_chip_target
        %     T_chip_target = max(T_chip_target,Tchip_record(i-1,j))-3;
        % else
        %     %对于每一个j，都确定一个基准控温目标，在此设定控温目标
        %     % T_chip_target = 273-3+3*j+-40;
        %     T_chip_target = 280-3+3*j;
        % end
        T_chip_target = 294-3+3*j;
        %有了控温目标后开始非线性求解
        
        [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed3(par,T_g_in,i,j,T_chip_target,I_record);
        
        %在命令行输出i，j表征当前进度
        disp(i);
        disp(j);
        disp('当前控温目标为');
        disp(T_chip_target);
        disp('实际温度为');
        disp(T_chip);
        initial_change = 0;
        while flag <= 0 

            %如果控不住，换随机初始值再来一次
            % T_chip_target = T_chip_target + 1;
            initial_change = initial_change + 1;
            [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed3(par,T_g_in,i,j,T_chip_target,I_record);
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
        T_chip_target_record(i,j) = T_chip_target;
    
        i = i+1;
        rate = Dz/par.length;
        Dz = Dz+dz;
        
        T_g_in = Tg;
        % if i==154
        %     T_g_in = par.T_g_in;
        % else
        %     T_g_in = Tg;
        % end

    %结束沿程循环    
    end
    grad(:,j) = gradient(Tchip_record(:,j));
    
    if max(grad.^2) > (56.83*0.3)^2
        grad_record(j) = 0;
    else
        grad_record(j) = 1;

    end

    if max(Tchip_record(:,j)) <=273+65 && min(Tchip_record(:,j)) >= 273-40
        range_record(j) = 1;
    else
        range_record(j) = 0;
    end
    
end

mean_T_chip = mean(Tchip_record,1);
median_T_chip = median(Tchip_record,1);
geomean_T_chip = geomean(Tchip_record,1);
harmmean_T_chip = harmmean(Tchip_record,1);
range_T_chip = range(Tchip_record,1);
var_T_chip = var(Tchip_record,[],1);
std_T_chip = std(Tchip_record,0,1);
diary("off");