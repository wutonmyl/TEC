% 用于计算给定零电流工况下的相应梯度调控工况
diary('test_DiaryFile');
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
mode = 2;
h = par.h;
slope = par.slope;
% data_ctl = struct();
%导入零电流工况
data_zero = load('data_material_SCH4_zero.mat','data_zero');
data_zero = data_zero.data_zero;
%遍历工况结构体,构造一个结构体数组储存零电流工况
% data_zero_tem = struct();
% data_zero.case2_1_0.data_name = 'case2_1_0';
% data_zero_tem = data_zero.case2_1_0;

%构造data_ctl内部字段
% fileds = fieldnames(data_zero);
% for i = 1:length(fileds)
%     fileds_i = fileds(i);
%     key = fileds_i{1};
%     %给单个数据添加数据名字段
%     data_zero.(key).data_name = key;
%     % data_zero_tem(i) = data_zero.(key); 
%     key_ctl = replace(key,'0','1');
%     data_ctl.(key_ctl) = [];
% end

%% 采用并行计算求解不同梯度下对应工况的热性能
grad_lb = 0;
grad_ub = par.grad_ub;
grad_target = grad_lb:0.1:grad_ub;
%% 需要修改，使用的基础工况
data_zero_use = data_zero.case1_2_0;
%这里根据需要自行修改温度上界和温度下界
Tchip_zero_ub = mean(data_zero_use.Tchip_record);
% Tchip_zero_ub = max(data_zero_use.Tchip_record);
Tchip_zero_lb = min(data_zero_use.Tchip_record);
%% 创造一个和正常data数据同字段的空结构体用于构造结构体数组
%创造一个跟数据统一格式的结构体，用于创造替身结构体数组
origin_stc = current_control_set7(par,Tchip_zero_ub,Tchip_zero_lb,grad_ub,search_num,mode,h,slope);

%获取结构体的字段数据
%遍历结构体，清空数据，保留字段，减少运行内存
fileds = fieldnames(origin_stc);
for i = 1:length(fileds)
    fileds_i = fileds(i);
    key = fileds_i{1};
    origin_stc.(key) = [];
end
data_ctl_tem = repmat(origin_stc,[length(grad_target),1]);
%% 并行计算，将数据存入结构体数组中去
parfor i = 1:length(grad_target)
    [data_tem,~,~] = current_control_set7(par,Tchip_zero_ub,Tchip_zero_lb,grad_target(i),search_num,mode,h,slope);
    %随循环改变输出变量名
    data_ctl_tem(i) = data_tem;
end
%% 需要修改，将结构体数组数据导入到命名好的结构体数据中去
for i = 1:length(grad_target)
    expr_grad_name = replace(num2str(grad_target(i)),'.','_');
    expr = ['data_ctl.test.case2_1_1_strict.grad_' expr_grad_name ' = data_ctl_tem(i);'];
    eval (expr);
end

diary off;


