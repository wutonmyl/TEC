addpath('current_change\');
addpath("basic_fuction\");
par = para;
diary('myDiaryFile');

%构造数据记录
i_record = [];
Dz_record = [];
Tchip_record = [];
Tc_record   = [];
Th_record   = [];
To_record   = [];
Ti_record   = [];
Tg_record   = [];
h_record    = [];
q1_r_record = [];
q2_l_record = [];
q2_r_record = [];
q3_l_record = [];
q3_r_record = [];
q4_r_record = [];
q5_r_record = [];
flag_record = [];
n_record = [];
choice_record = [];
mean_T_chip = [];
median_T_chip = [];
geomean_T_chip = [];
harmmean_T_chip = [];
range_T_chip = [];
var_T_chip = [];


for j = 1:20
    dz = par.Dz;
    Dz = 0;
    i = 1;
    T_g_in = par.T_g_in;
    while Dz<=par.length
        [T_chip,Tc,Th,To,Ti,Tg,h,flag,n,choice] = current_influ_function_solver(par,T_g_in,i,j);
       
        q2_r = par.k_ct*par.a_te*(T_chip-Tc);
        q2_l = par.n*(par.alpha*current_change(j)*Tc-0.5*current_change(j)^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
        q3_r = par.k_ct*par.a_te*(Th-To);
        q3_l = par.n*(par.alpha*current_change(j)*Th+0.5*current_change(j)^2*par.R+par.k_n*par.a_copper*(Tc-Th)/par.delta_n);
        q4_r= par.k_channel*par.a_te*(To-Ti)/par.delta_channel;
        q5_r= h*par.a_g*(Ti-Tg);
        Q_tec = q3_l-q2_l;
        COP = q2_l/Q_tec;
        
        i_record(j,i) = i;
        Tchip_record(j,i) = T_chip;
        Tc_record(j,i) = Tc;
        Th_record(j,i) = Th;
        To_record(j,i) = To;
        Ti_record(j,i) = Ti ;
        Tg_record(j,i) = Tg;
        h_record(j,i) = h;
        
        q2_l_record(j,i) = q2_l;
        q2_r_record(j,i) = q2_r;
        q3_l_record(j,i) = q3_l;
        q3_r_record(j,i) = q3_r;
        q4_r_record(j,i) = q4_r;
        q5_r_record(j,i) = q5_r;
        flag_record(j,i) = flag;
        n_record(j,i) = n;
        choice_record(j,i) = choice;
        COP_record(j,i) = COP;
        Q_tec_record(j,i) = Q_tec;
        Dz_record(j,i) = Dz;
        i = i+1;
        rate = Dz/par.length;
        Dz = Dz+dz;
        
       
        T_g_in = Tg;
        disp(i);
        disp(j);
        disp(rate);
        
    
    
    end
end





diary("off");