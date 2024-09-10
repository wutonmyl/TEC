classdef para_water
    %PARA 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        % alpha = 2*301e-6
        alpha = 2*186e-6
        %电流基本工况为0
        I = 5
        % R = 2
        R = 0.028949942922041 
        % R = 0.023512940483788
        n = 14*7

        k_p = 1.2
        k_n = 1.2
        % k_p = 0.63
        % k_n = 0.63
        k_copper = 400
        k_ceramic = 1.5
        %这里算了综合热导率，而且把厚度算进去了，方程里面不需要厚度了
        k_ct = 118366.541724206

        delta_p = 80e-6
        delta_n = 80e-6
        delta_copper =46e-6
        delta_ceramic = 12.5e-6
        delta_ct = 5.85000000000000e-05
        lou_p = 1/(9.259e4)
        lou_n = 1/(9.259e4)
        % lou_p = 1/(11.4e4)
        % lou_n = 1/(11.4e4)
        lou_copper = 1/(5.988e7)

        a_p = 160.71*160.71e-12
        a_n = 160.71*160.71e-12
        a_copper = 371.42*160.71e-12
        a_te = 9e-6
        a_g = 7.5e-5
        a_cross = 2.5*10e-6
        l =10e-3
        w = 2.5e-3
        length = 870e-3

        %功率不确定
        % P_chip = 6 做一些修改
        P_chip = 6
        k_channel = 400
        delta_channel = 0.25e-3
        %流体运行参数
        P = 23e3
        mg = 0.007
        Dz = 3e-3
        % T_g_in = 185
        T_g_in = 273
        Tpc = 648
        R_channel = 1.5e-3
        



    end
    
    methods
        function obj = para()
            %PARA 构造此类的实例
            %   此处显示详细说明
           
        end
        
        % function x = judge_wan(j)
        %     flag_wan_record = [];
        %     for i = 1:290
        %         flag_wan = 0;
        %         if mod(i,17)==1||mod(i,17)==0
        %             flag_wan = 1;
        %         end
        %         flag_wan_record = [flag_wan_record,flag_wan];
        %     end
        %     x = flag_wan_record(j);
        % end
    end
end

