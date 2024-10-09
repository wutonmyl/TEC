%本脚本用于整合不同版本的数据
data_tem = load("material_0_SCH4_v_0.007_Tin_185K_初值随机_宽松调控_multigrad.mat","data_v0_Tin_185K_P_6W_uniform");
data_ctl.case1_0_1 = data_tem.data_v0_Tin_185K_P_6W_uniform;
