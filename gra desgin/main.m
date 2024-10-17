%计算基本工况的主函数，要计算包括基础工况的沿程热性能参数和梯度限制
%引入基本方程
addpath('basic_fuction\')
%引入数据类文件夹
addpath('parament\')
%设定基础参数
par = para_fixed;
dz = par.Dz;
Dz = dz;
i = 1;
T_g_in = par.T_g_in;
j = 1;
mode = 2;
slope = par.slope;
h = par.h;

current_exist = 0;
search_num =1;

% data_zero = struct();
%需要改数据昵称
data_zero.fixed.case4_2_0 = current_control_set1(par,mode,h,slope,search_num);






% data_name = data_name(par,current_exist,mode,h,slope)
% %构造数据记录
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
% data_tem = struct();
% 
% while Dz<=par.length
%     [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = function_solver(par,T_g_in,i,mode,h,slope);
% 
% 
%     q2_r = par.k_ct*par.a_te*(T_chip-Tc);
%     q2_l = par.n*(par.alpha*par.I*Tc-0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
%     q3_r = par.k_ct*par.a_te*(Th-To);
%     q3_l = par.n*(par.alpha*par.I*Th+0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
%     q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
%     q5_r= h*par.a_g*(Ti-Tg);
%     Q_tec = q3_l-q2_l;
%     Q_tec = abs(Q_tec);
%     COP = 1/(Q_tec/q2_l);
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
%     COP_record(i,j) = COP;
%     Q_tec_record(i,j) = Q_tec;
%     Dz_record(i,j) = Dz;
% 
%     i = i+1;
%     rate = Dz/par.length;
%     disp(rate)
%     Dz = Dz+dz;
% 
%     T_g_in = Tg;
% end
% % Tc_record = Tc_record';
% % Th_record = Th_record';
% % To_record = To_record';
% % Ti_record = Ti_record';
% % Tg_record = Tg_record';
% 
% % Tchip_record = Tchip_record';
% % h_record = h_record';
% mean_T_chip = mean(Tchip_record,1);
% median_T_chip = median(Tchip_record,1);
% geomean_T_chip = geomean(Tchip_record,1);
% harmmean_T_chip = harmmean(Tchip_record,1);
% range_T_chip = range(Tchip_record,1);
% var_T_chip = var(Tchip_record,[],1);
% std_T_chip = std(Tchip_record,0,1);
% 
% 
% 

% %无tec部分
% dz_notec = par.Dz;
% Dz_notec = 0;
% i_notec = 1;
% T_g_in_notec = par.T_g_in;
% %构造数据记录
% i_record_notec = [];
% Dz_record_notec = [];
% Tchip_record_notec = [];
% Tc_record_notec   = [];
% Th_record_notec   = [];
% To_record_notec   = [];
% Ti_record_notec   = [];
% Tg_record_notec   = [];
% h_record_notec    = [];
% q1_r_record_notec = [];
% q2_l_record_notec = [];
% q2_r_record_notec = [];
% q3_l_record_notec = [];
% q3_r_record_notec = [];
% q4_r_record_notec = [];
% q5_r_record_notec = [];
% flag_record_notec = [];
% n_record_notec = [];
% choice_record_notec = [];
% while Dz_notec<=par.length
%     [T_chip_notec,Tc_notec,Th_notec,To_notec,Ti_notec,Tg_notec,h_notec,flag_notec,n_notec,choice_notec] = function_solver_nocurrent(par,T_g_in_notec,i_notec);
% 
%     q2_r_notec = par.k_ct*par.a_te*(T_chip_notec-Tc_notec);
%     q2_l_notec = par.n*(par.alpha*par.I*Tc_notec-0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc_notec-Th_notec)/par.delta_n);
%     q3_r_notec = par.k_ct*par.a_te*(Th_notec-To_notec);
%     q3_l_notec = par.n*(par.alpha*par.I*Th_notec+0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc_notec-Th_notec)/par.delta_n);
%     q4_r_notec= par.k_channel*par.a_te*(To_notec-Ti_notec)/par.delta_channel;
%     q5_r_notec= h_notec*par.a_g*(Ti_notec-Tg_notec);
% 
%     i_record_notec = [i_record_notec,i_notec];
%     Tchip_record_notec = [Tchip_record_notec,T_chip_notec];
%     Tc_record_notec = [Tc_record_notec,Tc_notec];
%     Th_record_notec = [Th_record_notec,Th_notec];
%     To_record_notec = [To_record_notec,To_notec];
%     Ti_record_notec = [Ti_record_notec,Ti_notec];
%     Tg_record_notec = [Tg_record_notec,Tg_notec];
%     h_record_notec = [h_record_notec,h_notec];
% 
%     q2_l_record_notec = [q2_l_record_notec,q2_l_notec];
%     q2_r_record_notec = [q2_r_record_notec,q2_r_notec];
%     q3_l_record_notec = [q3_l_record_notec,q3_l_notec];
%     q3_r_record_notec = [q3_r_record_notec,q3_r_notec];
%     q4_r_record_notec = [q4_r_record_notec,q4_r_notec];
%     q5_r_record_notec = [q5_r_record_notec,q5_r_notec];
%     flag_record_notec = [flag_record_notec, flag_notec];
%     n_record_notec = [n_record_notec,n_notec];
%     choice_record_notec = [choice_record_notec,choice_notec];
% 
%     i_notec = i_notec+1;
%     rate = Dz_notec/par.length
%     Dz_notec = Dz_notec+dz_notec;
%     Dz_record_notec = [Dz_record_notec,Dz_notec];
% 
%     T_g_in_notec = Tg_notec;
% 
% 
% 
% end
