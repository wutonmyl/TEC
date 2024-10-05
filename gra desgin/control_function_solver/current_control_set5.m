function [data,wrong_num_record,i_final] = current_control_set5(par,ub,lb,grad_bond,search_num,mode,h,slope)
%%  本函数是用于整合原main5和main6函数的所有求解外设置部分，将求解流程整流为一整个输入输出函数，目的是用于新main5整体优化求取最小梯度限制值
%   此处显示详细说明


%%设置迭代环节基础参数
dz = par.Dz;
Dz = dz;
i = 1;
j = 1 ;
i_final = i;
T_g_in = par.T_g_in;
% mode = 1;
% h = par.h;
% slope = par.slope;
% search_num = 1;


%% 构造数据记录
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
wrong_num_record = zeros(289,search_num);

% i_special_record = zeros(289,search_num);


%% 开始沿程循环
while Dz<=par.length
    
    %% 默认用无电流模式开始工作
    [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = function_solver(par,T_g_in,i,mode,h,slope);
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

    %% 开始判断是否满足控温标准
    
    [need_control,Tchip_target_recmd] = judge_need_control(Tchip_record,i,j,lb,ub,grad_bond);
    need_control_record(i,j) = need_control;
   
    T_chip_target_record(i,j) = Tchip_target_recmd;
    %命令行显示进度
    disp(i);
    disp(j);
    disp(['推荐温度为:',num2str(Tchip_target_recmd)]);
    disp(['梯度限制为：',num2str(grad_bond)]);
    % disp(Tchip_target_recmd)
    disp(['实际温度为:',num2str(T_chip)]);
    % disp(T_chip);
    %% 如果需要控温，采用独立调控程序
    contol_repeat = 0;
    initial_change = 0;
    while need_control ~= 0
        % disp('进入独立调控环节')
        contol_repeat = contol_repeat + 1;
        disp(['进入独立调控环节,当前重复调控次数:',num2str(contol_repeat)])
        % disp(contol_repeat);
        [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed5(par,T_g_in,i,j,Tchip_target_recmd,I_record,mode,h,slope);
        disp(i);
        disp(j);
        disp(['推荐温度为:',num2str(Tchip_target_recmd)]);
        disp(['梯度限制为：',num2str(grad_bond)]);
        % disp(Tchip_target_recmd)
        disp(['实际温度为:',num2str(T_chip)]);
        % disp(T_chip);
        disp(['不满足记号为:',num2str(need_control)]);
        % disp(need_control);
        initial_change = 0;


        %求解器抱错误解进入循环
        while flag <= 0 
    
            %如果控不住，换随机初始值再来一次
            % T_chip_target = T_chip_target + 1;
            initial_change = initial_change + 1;
            [I,Tc,Th,To,Ti,Tg,h,flag,n,choice,T_chip] = current_control_function_solver_fixed5(par,T_g_in,i,j,Tchip_target_recmd,I_record,mode,h,slope);
            disp(i);
            disp(j);
            disp(['更换随机初始值数:',num2str(initial_change)]);
            % disp(initial_change);
            disp(['实际温度为:',num2str(T_chip)]);
            % disp(T_chip);
            %如果一直控不住，大于五次循环的时候记录了一下wrong_num,然后跳出循环
            if initial_change >= 5
                wrong_num_record(i,j) = wrong_num_record(i,j) + 200;
                disp('求解器一直无解，跳出更换随机初值循环')
                break;
            end
        end

        %% 开始数据记录
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
        %结束数据记录
        %求解器一直无解，就跳出判定是否需要调控循环
        if initial_change >= 5 
            disp('求解器一直无解，跳出循环')
            break;
        elseif contol_repeat >= 10
            disp('没办法把梯度调到限度一下，跳出判定是否需要调控循环')
            wrong_num_record(i,j) = wrong_num_record(i,j) + 985;
            break;
            
        end

        %再检查一遍是否需要调控
        [need_control,Tchip_target_recmd] = judge_need_control(Tchip_record,i,j,lb,ub,grad_bond);
        T_chip_target_record(i,j) = Tchip_target_recmd;
    end
    %% 结束一个单元的计算
    %求解器一直无解，就跳出沿程循环
    if initial_change >= 5 
        disp('求解器一直无解，跳出循环')
        break;
    elseif contol_repeat >= 10
        disp('没办法把梯度调到限度一下，跳出判定是否需要调控循环')
        
        break;
        
    end
    i = i+1;
    rate = Dz/par.length;
    Dz = Dz+dz;
    
    T_g_in = Tg;
    

end
i_final = i;
%% 数据处理
grad_T_chip = gradient(Tchip_record(:,j));
mean_T_chip = mean(Tchip_record(:,j));
median_T_chip = median(Tchip_record(:,j));
geomean_T_chip = geomean(Tchip_record(:,j));
harmmean_T_chip = harmmean(Tchip_record(:,j));
range_T_chip = range(Tchip_record(:,j));
var_T_chip = var(Tchip_record(:,j));

std_T_chip = std(Tchip_record(:,j));
grad_Power_rate = (56.83*0.3-grad_bond)/sum(Q_tec_record(:,j));
grad_Power_product = grad_bond*sum(Q_tec_record(:,j));

%将数据放入结构体
data.Q_tec_record = Q_tec_record;
data.COP_record = COP_record;
data.i_record = i_record;
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
data.grad_record = grad_record;
% data.range_record = range_record;
data.Cp_g_record = Cp_g_record;
data.density_g_record = density_g_record;
data.H_g_record = H_g_record;
data.lamda_g_record = lamda_g_record;
data.Pr_g_record = Pr_g_record;
data.mu_g_record = mu_g_record;
data.T_chip_target_record = T_chip_target_record;
data.need_control_record = need_control_record;
data.grad_T_chip = grad_T_chip;
data.wrong_num_record = wrong_num_record;
data.mean_T_chip = mean_T_chip;
data.median_T_chip = median_T_chip;
data.geomean_T_chip = geomean_T_chip;
data.harmmean_T_chip = harmmean_T_chip;
data.range_T_chip = range_T_chip;
data.var_T_chip = var_T_chip;
data.std_T_chip = std_T_chip;
data.i_final = i_final;
data.grad_Power_rate = grad_Power_rate;
data.grad_Power_product = grad_Power_product;
data.grad_bond = grad_bond;



end