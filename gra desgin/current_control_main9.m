%% 本主函数是用于尝试大规模优化变量优化tec沿程电流的问题
%% 读取基础参数，基础函数
%引入热量变化函数库，确保可以计算不同热量分布的工况
addpath(['heat_change\'])
%引入基本方程，确保可以计算对流换热系数
addpath('basic_fuction\')
%引入参数类文件夹
addpath('parament\')
%引入求解器文件夹
addpath("control_function_solver\")
par = para_fixed;
search_num = 1;
mode = 2;
h = par.h;
slope = par.slope;
%导入基础工况，给出温度上下界
data_zero = load('data_material_SCH4_zero.mat','data_zero');
data_zero = data_zero.data_zero;
data_zero_use = data_zero.fixed.case1_2_0;
Tchip_zero_ub = data_zero_use.Tchip_record(289,1);
Tchip_zero_lb = min(data_zero_use.Tchip_record);
grad_lb = 0;
grad_ub = par.grad_ub;

I = current_control_set9(par,mode,h,slope,Tchip_zero_lb,Tchip_zero_ub);

[data,i_final,break_flag,Tchip_record,grad_T_chip,Q_tec_overall,COP_record] = current_control_function_solver_fixed9(par,I);