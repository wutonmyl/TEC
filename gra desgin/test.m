par = para;
% % x = [4,4,9,9,14,14];
% % y = [5,13,5,13,5,13];
% % %生成山峰的高度
% % h = ones(1,6);
% % %生成山峰的坡度
% % xs = 2*ones(1,6);
% % ys = 2*ones(1,6);
% %生成地形数据
% x = [14];
% y = [14];
% %生成山峰的高度
% h = ones(1,1);
% Z = ones(171,171);
% [X,Y] = meshgrid(0:0.1:17);
% for k1 = 1
%     Z =  Z + h(k1)*exp(-((X - x(k1))/xs(k1)).^2 - ((Y - y(k1))/ys(k1)).^2);
% end
% P_total = par.P_chip*171*171;
% P_tem = sum(Z(:));
% Z = Z*P_total/P_tem;
% mesh(Z)
% 
% 
% % par = para;
% I_max = par.alpha*Tc_record/par.R;
% Q_max = par.n*(par.alpha*I_max.*Tc_record-par.R*I_max.^2/2);
% Q_real = par.n*(par.alpha*I_record.*Tc_record-par.R*I_record.^2/2);
% diff = Q_max-Q_real;
% dletaT_max = (diff-6)/(par.n*par.k_p*par.a_copper/par.delta_p);
% par = para;
% 
% % [Cp_g,density_g,H_g,mu_g,lamda_g,Pr_g]=refpropm('CDHVL^','T',Tg_record,'P',par.P,'methane')
% a = [13;3;3];
% var(a(1:2,1))
% b = [13,3];
% % var(b)
% 
% a = [1 3 5 6 8 3 34 23 1];
% gradient(a)



% Tchip_zero_record = load('material_0_SCH4_v_0.007_Tin_185K_初值随机_I_0_entrance.mat','Tchip_record');
% 
% Tchip_zero_ub = max(Tchip_zero_record.Tchip_record);
% Tchip_zero_lb = min(Tchip_zero_record.Tchip_record);
% [Q_tec_record,COP_record,i_record,Dz_record,Tchip_record,Tc_record,Th_record,To_record,Ti_record,Tg_record,...
% h_record,q1_r_record,q2_l_record,q2_r_record,q3_l_record,q3_r_record,q4_r_record,q5_r_record,flag_record,...
% n_record,choice_record,I_record,grad_record,range_record,Cp_g_record,density_g_record,H_g_record,...
% lamda_g_record,Pr_g_record,mu_g_record,T_chip_target_record,grad,need_control_record,grad_T_chip,...
% wrong_num_record,mean_T_chip,median_T_chip,geomean_T_chip,harmmean_T_chip,range_T_chip,...
% var_T_chip,std_T_chip,i_final] = current_control_set5(par,Tchip_zero_ub,Tchip_zero_lb,grad_bond,search_num);

% x = [4];
% y = [4];
% %生成山峰的高度
% h = 2*ones(1,1);
% %生成山峰的坡度
% xs = 2*ones(1,1);
% ys = 2*ones(1,1);
% %生成地形数据
% Z = ones(17,17);
% [X,Y] = meshgrid(1:17);
% for k1 = 1
%     Z =  Z + h(k1)*exp(-((X - x(k1))/xs(k1)).^2 - ((Y - y(k1))/ys(k1)).^2);
% end
% P_total = par.P_chip*17*17;
% P_tem = sum(Z(:));
% Z = Z*P_total/P_tem;

% %读取结构体字段的程序
% grad_list = [];
% product_list = [];
% Qtec_list = [];
% a = data_v0_Tin_185K_P_6W_exit;
% list = a.grad_17_014;
% fileds = fieldnames(a);
% for i = 1:length(fileds)
%     fileds_i = fileds(i);
%     key = fileds_i{1};
%     list(i) = a.(key);
%     Qtec_list(i) = sum(a.(key).Q_tec_record);
%     grad_list(i) = a.(key).grad_bond ;
%     product_list(i) = a.(key).grad_Power_product;
% end
% tiledlayout(2,1)
% nexttile(1)
% plot(grad_list,Qtec_list);
% title('tec功率');
% nexttile(2)
% plot(grad_list,product_list);
% title('功率梯度');
% 

%比较两种梯度限制下的控温目标
% data_1 = data_ctl.case1_1_1.grad_8_5849;
% data_2 = data_ctl.case1_1_1.grad_4_8849;
% Tchip_1 = data_ctl.case1_1_1.grad_8_5849.T_chip_target_record;
% Tg_1 = data_1.Tg_record;
% Tg_2 = data_2.Tg_record;
% h_1 = data_1.h_record;
% h_2 = data_2.h_record
% Tchip_2 = data_ctl.case1_1_1.grad_4_8849.T_chip_target_record;
% tiledlayout(2,4);
% nexttile(1)
% plot(Tchip_1);
% title('梯度限制8.5849控温目标');
% nexttile(2)
% plot(Tg_1(1:65));
% title('梯度限制8.5849流体温升');
% nexttile(3);
% plot(h_1(1:65));
% title('梯度限制8.5849对流换热系数');
% nexttile(4)
% plot(data_1.Cp_g_record);
% title('梯度限制8.5849比热容');
% nexttile(5)
% plot(Tchip_2(1:65));
% title('梯度限制4.8849控温目标');
% nexttile(6);
% plot(Tg_2(1:65));
% title('梯度限制4.8849流体温度');
% nexttile(7);
% plot(h_2(1:65));
% title('梯度限制4.8849对流换热系数');
% nexttile(8);
% plot(data_2.Cp_g_record(1:65));
% title('梯度限制4.8849比热容');

%% 测试gpu

% i_record = zeros(289,1,'gpuArray');
% % i_record = gpuArray(i_record);
% for i = 1:289
% 
%     i_record(i,1) = i;
% end

disp('dddd\bf')