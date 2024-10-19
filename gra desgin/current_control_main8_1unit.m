%% 本函数是用于计算一个tec传热单元下的工况表现
%要求本函数拥有足够的输入参量可定义性，同时具有大规模并行计算的能力
par = para_fixed;
%%设置迭代环节基础参数
dz = par.Dz;
Dz = dz;
i = 1;
j = 1 ;
i_final = i;


%% 设置计算环境参数
%功率器件热流
p_chip = 6;
%质量流率
mg = 0.007*0.5;
%压力
pressure = 5e3;
%入口温度
T_g_in = 195;
%电流调控模式
current_mode = 0;

%%将数据导入设置函数，由设置函数进行导入计算
data = current_control_set8(current_mode,pressure,mg,p_chip,par,T_g_in);

data_one_unit.case1_0_0 = data;
% save('data_material_SCH4_ctl.mat','data_ctl','-append');
