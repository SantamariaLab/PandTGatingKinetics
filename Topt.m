function Topt = Topt(P, Po, To, DCp, DHo, DVo, Da, Dk, R)
% Calculate the optimal tempertature as a function of pressure
    Topt = (DCp*To-DHo+Da*To*(P-Po)-DVo*(P-Po)+Dk/2*(P-Po).^2)/(DCp+R);
end
