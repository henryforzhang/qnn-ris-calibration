function ris = gear2phase(risPsSet,gearIdx)
len = length(gearIdx);
ris = zeros(len,1);
for ll = 1 : len
    ris(ll) = risPsSet(ll,gearIdx(ll));
end
end