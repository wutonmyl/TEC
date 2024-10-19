function [alpha_i,n,choice] = Heat_transfer_coff_SCH4(Ti,Tg,i,par)
%HEAT 计算流体侧对流传热系数
%   此处显示详细说明
[Cp_g,density_g,H_g,mu_g,lamda_g,Pr_g]=refpropm('CDHVL^','T',Tg,'P',par.P,'methane');
[density_i,H_i]=refpropm('DH','T',Ti,'P',par.P,'methane');
de = 4.*par.l.*par.w./(2.*(par.l+par.w));
u_g = (par.mg./density_g)./(par.a_cross);%计算流速

Re_g = density_g.*u_g.*de./mu_g;
Cp_ = (H_i-H_g)./(Ti-Tg);

% if (Tg<Ti&&Ti<=par.Tpc)||(1.2.*par.Tpc<=Tg&&Tg<Ti)
%     n = 0.4;
%     choice = 1;
% elseif Tg<par.Tpc&&par.Tpc<=Ti
%     n = 0.4+0.2.*(Ti./par.Tpc-1);
%     choice = 2;
% elseif (par.Tpc<Tg&&Tg<=1.2.*par.Tpc)
%     n = 0.4+0.2.*(Ti./par.Tpc-1).*(1-5.*(Tg./par.Tpc-1));
%     choice = 3;
% else
%     n = 0.4;
%     choice = 4;
% end
            
n = 0.4;
% choice = 0;
Nu_in = 0.0183.*Re_g.^0.82.*Pr_g.^0.5.*(density_i./density_g).^0.3.*(Cp_./Cp_g).^n;%.*(1+1.7.*(de./par.R_channel)^3)^judge_wan(i);%努塞尔数
alpha_i = Nu_in.*lamda_g./de;%对流传热系数
end

