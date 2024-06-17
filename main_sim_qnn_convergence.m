%% full duplex communication system
%
%
clear;
% basic setting
%% RIS setting
MrisSet = [16,32,64];
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

% snrDbSet = [0 10 20 30];
snrDbSet = [30];
snrDbSetLen = length(snrDbSet);

b = 4;
errStd = 20;
Nb = 2^b;
O = 15;
oNb = O*Nb;
stdPsSet = (0:Nb-1)/Nb*2*pi;
stdPhsSet = exp(1j*stdPsSet);
risGeniePhsSetMax = gen_ris_phaseSet(MrisMax, b, errStd);
permtrSetMax = gen_permtr_mtr(MrisMax,Nb,O);
[risEmplyedMax,gearEmplyedMax] = nn_gen_multi_d(risGeniePhsSetMax,permtrSetMax);
%% channel
Nc = 3; Nray = 5;
[HbrRMax,HrBtMax] = save_channel(MrisSet,Mr,Nc,Nray,8);
plt2 = {'r-','k--','b:'};
plt = {'a','b','c'};
figure;
unitNoise = sqrt(1/2)*(randn(Mr,oNb)+1j*randn(Mr,oNb));
for ss = 1 : snrDbSetLen
    snrDb = snrDbSet(ss);
    noisePow = db2pow(-(snrDb+gain));
    for mm = 1 : MrisSetLen
        Mris = MrisSet(mm);
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
        Heff = HeffDeNoise+sqrt(noisePow)*unitNoise;
        [risEstPhsSet, convSet]= qnn_cali_ris_conv(Heff,gearEmplyed,risEstEmplyed,estPsSet);
        if ss == 1
            plt{mm} = semilogy(convSet,plt2{mm},'linewidth',2,'markersize',8);
            hold on
        else
            semilogy(convSet,plt2{mm},'linewidth',2,'markersize',8);
            hold on
        end
    end
end
%% plot
xlim([1, 2000])
grid on
xlabel('Number of Epochs');
ylabel('Cost Function: C_{ave}');
title([num2str(b),'-bit RIS, ', 'Phase Deviations [-',num2str(errStd),'\circ',', ',num2str(errStd),'\circ',']'])
legend(['a','b','c'],{'M_{ris} = 16','M_{ris} = 32','M_{ris} = 64'}, 'FontSize', 12)









