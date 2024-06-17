%% full duplex communication system
%
%
clear;
% basic setting
%% RIS setting
MrisSet = [16 32 64];
MrisSetLen = length(MrisSet);
MrisMax = max(MrisSet);
%% Bs setting
Mr = 4;
% MrSetMax = max(MrSet);
Mt = 1;
% pathloss
pilotLen = 100;
pilotGain = pow2db(pilotLen);
noiseGain = 3; % dB
gain = pilotGain-noiseGain;

snrDbSet = [0:5:30];
snrDbSetLen = length(snrDbSet);

% load('testdata0804.mat','HrBt','HbrR','Had')

b = 4;
errStd = 20;
Nb = 2^b;
O = 15;
oNb = O*Nb;
stdPsSet = (0:Nb-1)/Nb*2*pi;
stdPhsSet = exp(1j*stdPsSet);


risGeniePhsSetMax = gen_ris_phaseSet(MrisMax, b, errStd);
permtrSetMax = gen_permtr_mtr(MrisMax,Nb,O);
% risEmplyedMax = gen_multi_d(risGeniePhsSetMax,permtrSetMax);
[risEmplyedMax,gearEmplyedMax] = nn_gen_multi_d(risGeniePhsSetMax,permtrSetMax);
%% channel
% HrBt = sqrt(1/2)*(randn(Mris,Mt)+1j*randn(Mris,Mt));
% HbrR = sqrt(1/2)*(randn(Mr,Mris)+1j*randn(Mr,Mris));
% Had = sqrt(1/2)*(randn(Mr,Mt)+1j*randn(Mr,Mt));
load('channelData030623_1.mat');
Nc = 3;
Nray = 5;

mseSet = zeros(MrisSetLen,snrDbSetLen);
crbSet = zeros(MrisSetLen,snrDbSetLen);

nnMonte = 1e2;
for nn = 1 : nnMonte
    unitNoise = sqrt(1/2)*(randn(Mr,oNb)+1j*randn(Mr,oNb));
    for mm = 1 : MrisSetLen
        Mris = MrisSet(mm);
        if nn == 1 || mod(nn,5) == 0
            fprintf(['\nThe Iteration Number: %d, the number of RIS: %d, at ', datestr(now,"HH:MM"),'\n'], nn, Mris);
        end
        HrBt = HrBtMax{mm};
        HbrR = HbrRMax{mm};
        Hp = diag(HrBt);
        Hbreve = HbrR*Hp;
        risGeniePhsSet = risGeniePhsSetMax(1:Mris,:);
        permtrSet = permtrSetMax(:,1:Mris,:,:);
        risEmplyed = risEmplyedMax(1:Mris,:);
        gearEmplyed = gearEmplyedMax(1:Mris,:);
        HeffDeNoise = Hbreve*risEmplyed;
        estPsSet = init_estPsSet(Mris,Nb);
        estPhsSet = exp(1j*estPsSet);
        risEstEmplyed = gen_multi_d(estPhsSet,permtrSet);    
        for ss = 1 : snrDbSetLen
            snrDb = snrDbSet(ss);
            noisePow = db2pow(-(snrDb+gain));
            if nn == 1
                fisher = cal_CRB_of_RIS_cali(permtrSet,risGeniePhsSet,Hbreve);
                fisherInv = fisher^(-1);
                fisherInvDiag = diag(fisherInv);
                crbSet(mm,ss) = noisePow/2*mean(fisherInvDiag(1:(Mris*(Nb-1))));
            end
            Heff = HeffDeNoise+sqrt(noisePow)*unitNoise;
            risEstPhsSet =  qnn_cali_ris(Heff,gearEmplyed,risEstEmplyed,estPsSet);
            risEstPhsSetRmAm = remove_ambiguity(risEstPhsSet);
            mse = norm(angle(risEstPhsSetRmAm./risGeniePhsSet),'fro')^2/(Mris*(Nb-1));
            mseSet(mm,ss) = mseSet(mm,ss)+ mse;
        end
    end
end
mseSetAve = mseSet/nn;
rmseSet = sqrt(mseSetAve)*180/pi;
rcrbSet = sqrt(crbSet)*180/pi;
%% plot
plt1 = {'ro','kx','b^'};
plt2 = {'r-','k--','b:'};
figure;
for mm = 1 : MrisSetLen
    semilogy(snrDbSet,rmseSet(mm,:),plt1{mm},'linewidth',2,'markersize',8);
    hold on
    semilogy(snrDbSet,rcrbSet(mm,:),plt2{mm},'linewidth',2,'markersize',8);
end
grid on
xlabel('SNR (dB)');
ylabel('RMSE (^\circ)');
title([num2str(b),'-bit RIS, ', 'Phase Deviations [-',num2str(errStd),'\circ',', ',num2str(errStd),'\circ',']'])
legend('M_{ris} = 16, Estimate','M_{ris} = 16, CRB',...
    'M_{ris} = 32, Estimate','M_{ris} = 32, CRB',...,
    'M_{ris} = 64, Estimate','M_{ris} = 64, CRB')









