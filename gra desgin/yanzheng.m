N = 31;
S = 0.0002;
r = 0.00001;
K = 1.5;
G = 0.01196;
I = [74.75 79.45 84.5];
Tc = [40 60 80];
deltaT = [-20 -10 0 10 20 30 40 50];
COP = [];
for i = 1:3
    for j =1:8
        Qc(j,i) = 2*N*(S*I(i)*Tc(i)-0.5*I(i)^2*r/G-K*G*deltaT(j));
        Qp(j,i) = 2*N*(I(i)^2*r*G+S*I(i)*deltaT(j));
        COP(j,i) = Qc(j,i)/(Qc(j,i)+Qp(j,i));
    end
end
