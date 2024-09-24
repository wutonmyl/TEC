par = para;
% x = [4,4,9,9,14,14];
% y = [5,13,5,13,5,13];
% %生成山峰的高度
% h = ones(1,6);
% %生成山峰的坡度
% xs = 2*ones(1,6);
% ys = 2*ones(1,6);
%生成地形数据
x = [14];
y = [14];
%生成山峰的高度
h = ones(1,1);
Z = ones(171,171);
[X,Y] = meshgrid(0:0.1:17);
for k1 = 1
    Z =  Z + h(k1)*exp(-((X - x(k1))/xs(k1)).^2 - ((Y - y(k1))/ys(k1)).^2);
end
P_total = par.P_chip*171*171;
P_tem = sum(Z(:));
Z = Z*P_total/P_tem;
mesh(Z)


% par = para;
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
% h = ones(1,1);
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
