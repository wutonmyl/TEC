% plot(Dz_record,Tchip_record);
% plot(Dz_record,Tc_record);
% plot(Dz_record,Th_record);
% plot(Dz_record,To_record);
% plot(Dz_record,Ti_record);
%  plot(Dz_record,Tg_record);
% plot(Dz_record,h_record);
% plot(Dz_record,n_record);
% plot(Dz_record,choice_record);

% plot(Dz_record,Tchip_record,Dz_record,Tc_record,'g',Dz_record,Th_record,'*',Dz_record,To_record,'b--o',Dz_record,Ti_record,'b',Dz_record,Tg_record,'r')
% plot(Dz_record,q2_l_record,'g',Dz_record,q2_r_record,'g',Dz_record,q3_l_record,'b--o',Dz_record,q3_r_record,'b--o',Dz_record,q4_r_record,'c*',Dz_record,q5_r_record,'b');
% plot(Dz_record,q2_l_record,'g',Dz_record,q2_r_record,'g')
% Tchip_record = Tchip_record';
% Tc_record = Tc_record';
% Th_record = Th_record';
% To_record = To_record';
% Ti_record = Ti_record';
% Tg_record = Ti_record';
% h_record = h_record';
data_single = data_v0_Tin_185K_P_6W_uniform.grad_0_26465
tiledlayout(8,3);
% nexttile
% plot(exitflag_record)
% title('外部求解器返回值')
nexttile(1)
%左一
num = 1;
plot(data_single.h_record(:,num))
title('对流换热系数h (W\cdotm^{-2}\cdotK^{-1})')
nexttile(2)
%右一
plot(data_single.T_chip_target_record(:,num))
title('功率器件控温目标T-chip-target (K)')
nexttile(3)
plot(data_single.Ti_record(:,num)-data_single.Tg_record(:,num))
title(['内侧壁面与流体温度差 （K）'])
nexttile(4)

%左二
plot(data_single.Cp_g_record(:,num))
title("流体比热容Cp (J\cdotkg^{-1}\cdotK^{-1})")
nexttile(5)
%右二
plot(data_single.Tchip_record(:,num))
title('功率器件实际温度Tchip (K)')
nexttile(6)
plot(gradient(data_single.Tchip_record(:,num))/0.3)
title('功率器件实际温度梯度 Tgrad (K/cm)')
nexttile(7)

%左三
plot(data_single.lamda_g_record(:,num))
title('流体热导率\lambda (W\cdotm^{-1}\cdotK^{-1})')
nexttile(8)
%右三
plot(data_single.Tc_record(:,num))
title('冷端温度Tc (K)')
nexttile(10)

%左四
plot(data_single.choice_record(:,num))
title('传热关联式选项n')
nexttile(11)
%右四
plot(data_single.Th_record(:,num))
title("热端温度Th (K)")
nexttile(13)

%左五
plot(data_single.I_record(:,num))
title('电流I (A)')
nexttile(14)
%右五
plot(data_single.To_record(:,num))
title('热沉外侧温度To (K)')
nexttile(16)

%左六
plot(data_single.Q_tec_record(:,num))
title('沿程TEC功率W (W)')
nexttile(17)
%右六
plot(data_single.Ti_record(:,num))
title("热沉内侧温度Ti (A)")
nexttile(19)

%左七
plot(data_single.COP_record(:,num))
title('COP')
nexttile(20)
%右七
plot(data_single.Tg_record(:,num))
title('流体温度Tf (K)')
nexttile(22)

%左八
plot(data_single.q4_r_record(:,num))
title('热端热流Qh (W)')

nexttile(23)
%右八

plot(data_single.flag_record(:,num))
title('求解器返回值')



% nexttile
% plot(q3_l_record)
% title('3式左端热量计算')
% nexttile
% plot(q3_r_record)
% title('3式右端热量计算')
% 
% nexttile
% plot(q4_r_record)
% title('4式右端热量计算')
% nexttile
% plot(q5_r_record)
% title('5式右端热量计算')





% nexttile
% plot(Re_record)
% title('雷诺数')
% nexttile
% plot(Nu_record)
% title('努塞尔数')

% plot(Dz_record,Tchip_record,Dz_record,Tc_record,'g',Dz_record,Th_record,'*',Dz_record,To_record,'b--o',Dz_record,Ti_record,'b',Dz_record,Tg_record,'r')
% title('各个层级温度比较图')
% nexttile
% plot(Dz_record,q2_r_record,Dz_record,q2_l_record);
% title('二式两端热量比较图')
% nexttile
% Tchip_matrix = zeros(17,17);
% for j =1:17
% 
% 
%     if mod(j,2)==1
%         for k=1:17
%             Tchip_matrix(j,k) = Tchip_record(17*(j-1)+k);
%         end
%     else
%         for k =1:17
%             Tchip_matrix(j,18-k) = Tchip_record(17*(j-1)+k);
%         end
%     end
% end
% [X,Y] = meshgrid(1:17,1:17);
% meshz(X,Y,Tchip_matrix);
% title('功率器件温度三维图');
% nexttile
% for j =1:17
% 
% 
%     if mod(j,2)==1
%         for k=1:17
%             Tchip_matrix_notec(j,k) = Tchip_record_notec(17*(j-1)+k);
%         end
%     else
%         for k =1:17
%             Tchip_matrix_notec(j,18-k) = Tchip_record_notec(17*(j-1)+k);
%         end
%     end
% end
% [X,Y] = meshgrid(1:17,1:17);
% meshz(X,Y,Tchip_matrix_notec);
% title('无tec功率器件温度三维图');