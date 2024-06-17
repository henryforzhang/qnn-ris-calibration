function y = init_estPsSet(Mris,Nb)
y = zeros(Mris,Nb);
stdPsSet = (0:Nb-1)/Nb*2*pi; % set as the initial value 
for ll = 1 : Mris
    y(ll,:) = stdPsSet;
end
end