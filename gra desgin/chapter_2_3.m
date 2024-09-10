Ic = 50;
Vce0_25_IGBT = 0.9034;
Vce0_25_FWD = 0.0033;
Rce0_25_IGBT = 0.01093;
Rce0_25_FWD = 0.030846;
Kv_IGBT = -0.0010339;
Kv_FWD = 0.0105752;
Kr_IGBT = 5.0677e-5;
Kr_FWD = -0.0002215;
Eon_IGBT = 0.46e-3;
Eoff_IGBT = 1.2e-3;
Eoff_FWD = 0.42;
sigma = 0.6;

P_cond_Tr = Ic*(Vce0_25_IGBT+Kv_IGBT*(Tj-25))*sigma/3+Ic^2*(Rce0_25_IGBT+Kr_IGBT*(Tj-25))*sigma/3;
P_cond_
P_total_tr1 = 