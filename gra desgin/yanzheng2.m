alpha = 0.0566;
k = 1.392;
R = 1.857;
Th = 50+273;
I = [2,4,6,8,10];
Qc_30 = 30;
Qc_40 = 40;
Qc_50 = 50;
for i =1:5
    delta_30(i) = (alpha*I(i)*Th-0.5*I(i)^2*R-Qc_30)/(k+alpha*I(i));
end
for i =1:5
    delta_40(i) = (alpha*I(i)*Th-0.5*I(i)^2*R-Qc_40)/(k+alpha*I(i));
end
for i =1:5
    delta_50(i) = (alpha*I(i)*Th-0.5*I(i)^2*R-Qc_50)/(k+alpha*I(i));
end
delta_30 = delta_30';
delta_40 = delta_40';
delta_50 = delta_50';