par = para;
dz_notec = par.Dz;
Dz_notec = 0;
i_notec = 1;
T_g_in_notec = par.T_g_in;
%构造数据记录
i_record_notec = [];
Dz_record_notec = [];
Tchip_record_notec = [];
Tc_record_notec   = [];
Th_record_notec   = [];
To_record_notec   = [];
Ti_record_notec   = [];
Tg_record_notec   = [];
h_record_notec    = [];
q1_r_record_notec = [];
q2_l_record_notec = [];
q2_r_record_notec = [];
q3_l_record_notec = [];
q3_r_record_notec = [];
q4_r_record_notec = [];
q5_r_record_notec = [];
flag_record_notec = [];
n_record_notec = [];
choice_record_notec = [];
while Dz_notec<=par.length
    [T_chip_notec,Tc_notec,Th_notec,To_notec,Ti_notec,Tg_notec,h_notec,flag_notec,n_notec,choice_notec] = function_solver(par,T_g_in_notec,i_notec);
   
    q2_r_notec = par.k_ct*par.a_te*(T_chip_notec-Tc_notec);
    q2_l_notec = par.n*(par.alpha*par.I*Tc_notec-0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc_notec-Th_notec)/par.delta_n);
    q3_r_notec = par.k_ct*par.a_te*(Th_notec-To_notec);
    q3_l_notec = par.n*(par.alpha*par.I*Th_notec+0.5*par.I^2*par.R+par.k_n*par.a_copper*(Tc_notec-Th_notec)/par.delta_n);
    q4_r_notec= par.k_channel*par.a_te*(To_notec-Ti_notec)/par.delta_channel;
    q5_r_notec= h_notec*par.a_g*(Ti_notec-Tg_notec);
    
    i_record_notec = [i_record_notec,i_notec];
    Tchip_record_notec = [Tchip_record_notec,T_chip_notec];
    Tc_record_notec = [Tc_record_notec,Tc_notec];
    Th_record_notec = [Th_record_notec,Th_notec];
    To_record_notec = [To_record_notec,To_notec];
    Ti_record_notec = [Ti_record_notec,Ti_notec];
    Tg_record_notec = [Tg_record_notec,Tg_notec];
    h_record_notec = [h_record_notec,h_notec];
    
    q2_l_record_notec = [q2_l_record_notec,q2_l_notec];
    q2_r_record_notec = [q2_r_record_notec,q2_r_notec];
    q3_l_record_notec = [q3_l_record_notec,q3_l_notec];
    q3_r_record_notec = [q3_r_record_notec,q3_r_notec];
    q4_r_record_notec = [q4_r_record_notec,q4_r_notec];
    q5_r_record_notec = [q5_r_record_notec,q5_r_notec];
    flag_record_notec = [flag_record_notec, flag_notec];
    n_record_notec = [n_record_notec,n_notec];
    choice_record_notec = [choice_record_notec,choice_notec];
    
    i_notec = i_notec+1;
    rate = Dz_notec/par.length
    Dz_notec = Dz_notec+dz_notec;
    Dz_record_notec = [Dz_record_notec,Dz_notec];
   
    T_g_in_notec = Tg_notec;
    

    
end