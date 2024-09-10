%计算对流传热系数
function [alpha,Nu,Cpf]=Heat_transfer_cofficiance_LNG(Tw,Tf)
global Pb mg d Tpc
Pb = 5000;
mg = 0.007;
d = 0.004;
Tpc = 193.3;
[Cpf,densityf,Hf,muf,lambdaf,Prf]=refpropm('CDHVL^','T',Tf,'P',Pb,'methane');%主流体的物性参数
[densityw,Hw]=refpropm('DH','T',Tw,'P',Pb,'methane');%内壁面的物性参数
u=(mg/densityf)/(pi*(d/2)^2);%流速
Ref=densityf*u*d/muf;%雷诺数
Cp_=(Hw-Hf)/(Tw-Tf);
if (Tf<Tw&&Tw<=Tpc)||(1.2*Tpc<=Tf&&Tf<Tw)
n=0.4;
else
if Tf<Tpc&&Tpc<=Tw
n=0.4+0.2*(Tw/Tpc-1);
else
if Tpc<Tf&&Tf<=1.2*Tpc
n=0.4+0.2*(Tw/Tpc-1)*(1-5*(Tf/Tpc-1));
else
n=0.4;
end
end
end
Nu=0.0183*Ref^0.82*Prf^0.5*(densityw/densityf)^0.3*(Cp_/Cpf)^n;%努塞尔数
alpha=Nu*lambdaf/d;%对流传热系数
if isnan(alpha)
alpha = 0;
end
a = isreal(alpha);
if a == 0
alpha = 0;
end
