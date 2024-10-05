%用于计算宽松调控下梯度限制内的各梯度下沿程热性能参数
%% 读取基础参数，基础函数
%引入热量变化函数库，确保可以计算不同热量分布的工况
addpath(['heat_change\'])
%引入基本方程，确保可以计算对流换热系数
addpath('basic_fuction\')
%引入参数类文件夹
addpath('parament\')
%引入求解器文件夹
addpath("control_function_solver\")
par = para;
search_num = 1;
mode = 1;
h = par.h;
slope = par.slope;

%导入零电流工况基础的沿程功率器件温度数据作为控制的上下限依据
%改数据计算要确认
%1. heat_change函数对应修改
Bond = load('material_0_SCH4_v_0.007_Tin_185K_初值随机_I_0_uniform.mat','Tchip_record','grad_bond');
Tchip_zero_ub = max(Bond.Tchip_record);
Tchip_zero_lb = min(Bond.Tchip_record);

%导入梯度限制
grad_lb = Bond.grad_bond;
grad_ub = par.grad_ub;

%% 构造搜索计算的目标梯度，在梯度上下限中取值，并运用for循环记录相应计算数据

grad_target = grad_lb:0.1:grad_ub;

%计算对应梯度下数据
data = struct();

%注v0 = 0.007
data_v0_Tin_185K_P_6W_uniform_tem = struct();

%创造一个跟数据统一格式的结构体，用于创造替身结构体数组
origin_stc = current_control_set5(par,Tchip_zero_ub,Tchip_zero_lb,grad_ub,search_num,mode,h,slope);

%获取结构体的字段数据
%遍历结构体，清空数据，保留字段，减少运行内存
fileds = fieldnames(origin_stc);
for i = 1:length(fileds)
    fileds_i = fileds(i);
    key = fileds_i{1};
    origin_stc.(key) = [];
end
%empty = load('Empty.mat');
% data_v0_Tin_185K_P_6W_uniform_tem = repmat(origin_stc,[length(grad_target),1]);
data_v0_Tin_185K_P_6W_uniform_tem = repmat(origin_stc,[length(grad_target),1]);


%采用并行计算for
parfor i = 1:length(grad_target)
    [data_tem,~,~] = current_control_set5(par,Tchip_zero_ub,Tchip_zero_lb,grad_target(i),search_num,mode,h,slope);
    %随循环改变输出变量名
    data_v0_Tin_185K_P_6W_uniform_tem(i) = data_tem;

end

%创造正主结构体
data_v0_Tin_185K_P_6W_uniform = struct();
%将替身的数据导入独立命名的正主中
for i = 1:length(grad_target)
    expr_grad_name = replace(num2str(grad_target(i)),'.','_');
    expr = ['data_v0_Tin_185K_P_6W_uniform.grad_' expr_grad_name ' = data_v0_Tin_185K_P_6W_uniform_tem(i);'];
    eval (expr);
   
end


