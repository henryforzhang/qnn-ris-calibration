function [risEmplyed,gearEmplyed] = nn_gen_multi_d(risPhsSet,permtrSet)
[O,Mris,Nb,~] = size(permtrSet);
oNb=O*Nb;
risEmplyed = zeros(Mris,oNb);
gearEmplyed = zeros(Mris,oNb);
const = (1:Nb)';
for oo = 1:O
    permtr = squeeze(permtrSet(oo,:,:,:));
    tmp1 = zeros(Mris,Nb);
    tmp2 = zeros(Mris,Nb);
    for mm = 1:Mris
       tmpPer = squeeze(permtr(mm,:,:));
       tmp1(mm,:) =  tmpPer*risPhsSet(mm,:).';
       tmp2(mm,:) =  tmpPer*const;
    end
    risEmplyed(:,(oo-1)*Nb+(1:Nb)) = tmp1;
    gearEmplyed(:,(oo-1)*Nb+(1:Nb)) = tmp2;
end
end