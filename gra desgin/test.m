% par = para;
% x = [4,4,9,9,14,14];
% y = [5,13,5,13,5,13];
% %生成山峰的高度
% h = ones(1,6);
% %生成山峰的坡度
% xs = 2*ones(1,6);
% ys = 2*ones(1,6);
% %生成地形数据
% Z = ones(171,171);
% [X,Y] = meshgrid(0:0.1:17);
% for k1 = 1:6
%     Z =  Z + h(k1)*exp(-((X - x(k1))/xs(k1)).^2 - ((Y - y(k1))/ys(k1)).^2);
% end
% P_total = par.P_chip*17*17;
% P_tem = sum(Z(:));
% Z = Z*P_total/P_tem;
% mesh(Z)


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
% var(b)

a = [1 3 5 6 8 3 34 23 1];
gradient(a)