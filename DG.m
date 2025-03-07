function DG = DG(P, T, Po, To, DCp, DSo, DHo, DVo, Da, Dk)
% Calculate the change in Gibbs free energy across activation barrier
    DG = DCp*(T-To)-DCp*T.*log(T/To)-T*DSo+Da*(T-To).*(P-Po)+DVo*(P-Po)-Dk/2*(P-Po).^2+DHo;
end