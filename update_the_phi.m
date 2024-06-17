function estPsSet = update_the_phi(estPsSet,g,phi)
gLen = length(g);

for gg = 1 : gLen
    estPsSet(gg,g(gg)) = phi(gg);
end
end