%独立调控的主程序

%引入热量变化函数库，确保可以计算不同热量分布的工况
addpath(['heat_change\'])
%引入基本方程，确保可以计算对流换热系数
addpath('basic_fuction\')
%引入参数类文件夹
addpath('parament\')
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
        T_chip_target = 290-3+3*j;
        %有了控温目标后开始非线性求解
        [I,Tc,Th,To,Ti,Tg,h,flag,n,choice] = current_control_function_solver_fixed2(par,T_g_in,i,j,T_chip_target);
        %在命令行输出i，j表征当前进度
        disp(i);
        disp(j);
        disp('当前控温目标为');
        disp(T_chip_target);
       
        while flag <= 0
    
            %如果控不住，得不到解放宽控温目标
            T_chip_target = T_chip_target + 3;
            [I,Tc,Th,To,Ti,Tg,h,flag,n,choice] = current_control_function_solver_fixed2(par,T_g_in,i,j,T_chip_target);
            disp(i);
            disp(j);
            disp('当前控温目标为');
            disp(T_chip_target);
        end
            
    
        q2_r = par.k_ct*par.a_te*(T_chip_target-Tc);
        q2_l = par.n*(par.alpha*I*Tc-0.5*I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
        q3_r = par.k_ct*par.a_te*(Th-To);
        q3_l = par.n*(par.alpha*I*Th+0.5*I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
        q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
        q5_r= h*par.a_g*(Ti-Tg);
        Q_tec = q3_l-q2_l;
        Q_tec = abs(Q_tec);
        COP = 1/(Q_tec/q2_l);
        
  
        i_record(i,j) = i;
        Tchip_record(i,j) = T_chip_target;
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
    
        i = i+1;
        rate = Dz/par.length;
        Dz = Dz+dz;
        
       
        T_g_in = Tg;
    %结束沿程循环    
    end
    grad = gradient(Tchip_record(:,j));
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