function k = k(T, G, kb, h, R)
% Calculate the rate coefficient
    k = (kb*T/h).*exp(-G./(R*T));
end
