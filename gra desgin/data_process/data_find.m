%本脚本用于获取梯度数据中可行的数据范围和最小的可行梯度，并用图像表示出来
data = data_ctl.test.COP.case1_2_1;
%获取data的字段用于遍历
fileds = fieldnames(data);
data_grad = [];
data_i_final = [];
data_flag = [];
for i = 1:length(fileds)
    fileds_i = fileds(i);
    key = fileds_i{1};
    data_i_final(i) = data.(key).i_final;
    data_grad(i) = data.(key).grad_bond;
    if data_i_final(i) == 290
        data_flag(i) = data_grad(i);
    else
        data_flag(i) = 20;
    end

end

min_grad = min(data_flag);
if min_grad == 20
    min_grad = 'fail';
end
plot(data_flag);

data_ctl.grad_process.case1_3_1.min_grad = min_grad;

data_ctl.grad_process.case1_3_1.data_flag = data_flag;

