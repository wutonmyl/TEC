% mean_T_chip = mean(Tchip_record,2);
% median_T_chip = median(Tchip_record,2);
% geomean_T_chip = geomean(Tchip_record,2);
% harmmean_T_chip = harmmean(Tchip_record,2);
% range_T_chip = range(Tchip_record,2);
% var_T_chip = var(Tchip_record,[],2);
% % std_T_chip = std(Tchip_record,0,2);
% Tc_record = Tc_record';
% Th_record = Th_record';
% To_record = To_record';
% Ti_record = Ti_record';
% Tg_record = Tg_record';
% 
% Tchip_record = Tchip_record';
% h_record = h_record';
% % Dz_record = Dz_record';
% % COP_record = COP_record';
% Q_tec_record = Q_tec_record';

% COP_record = 1./COP_record;
% Dz_record = 1000*Dz_record
% cc =  mean(I_record,1)
Th_record_matrx = reshape(Th_record,[17,17]);