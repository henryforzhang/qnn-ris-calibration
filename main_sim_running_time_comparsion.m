%% full duplex communication system
%
%
clear;
% basic setting
%% RIS setting
MrisSet = [8,16,32,64];
% MrisSet = [64];
MrisSetLen = length(MrisSet);
MrisMax = max(MrisSet);
%% Bs setting
Mr = 8;
% MrSetMax = max(MrSet);
Mt = 1;
% pathloss
pilotLen = 100;
pilotGain = pow2db(pilotLen);
noiseGain = 3; % dB
gain = pilotGain-noiseGain;

snrDb = 20;
noisePow = db2pow(-(snrDb+gain));
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
Nc = 3;
Nray = 5;
[HbrRMax,HrBtMax] = save_channel(MrisSet,Mr,Nc,Nray,8);
% save('HbrRMaxHrBtMax.mat','HbrRMax','HbrRMax')
% load('HbrRMaxHrBtMax.mat');
%% simulation
nnCostTimeSet = zeros(1,MrisSetLen);
nnRmseSet = zeros(1,MrisSetLen);
rcgLsCostTimeSet = zeros(1,MrisSetLen);
rcgLsRmseSet = zeros(1,MrisSetLen);

nnMonte = 1e1;
for nn = 1 : nnMonte
    if nn == 1 || mod(nn,10) == 0
        fprintf(['\nThe number of Iteration: %d, at ', datestr(now,"HH:MM"),'\n'], nn);
    end
    for mm = 1 : MrisSetLen
        Mris = MrisSet(mm);
        HrBt = HrBtMax{mm};
        HbrR = HbrRMax{mm};
        Hp = diag(HrBt);
        Hbreve = HbrR*Hp;
        risGeniePhsSet = risGeniePhsSetMax(1:Mris,:);
        estPsSet = init_estPsSet(Mris,Nb);
        estPhsSet = exp(1j*estPsSet);
        
        unitNoise = sqrt(1/2)*(randn(Mr,oNb)+1j*randn(Mr,oNb));
        permtrSet = permtrSetMax(1:O,1:Mris,:,:);
        risEmplyed = risEmplyedMax(1:Mris,1:oNb);
        risEstEmplyed = gen_multi_d(estPhsSet,permtrSet);
        gearEmplyed = gearEmplyedMax(1:Mris,1:oNb);
        HeffDeNoise = Hbreve*risEmplyed;
        %% running time test       
        Heff = HeffDeNoise+sqrt(noisePow)*unitNoise;
        [nnCostTime,nnRmse] =  qnn_cali_ris_test_running_time(Heff,gearEmplyed,risEstEmplyed,estPsSet,risGeniePhsSet);
        [rcgLsCostTime,rcgLsRmse] = rcgLs_cali_ris_test_running_time(Heff,O,Nb,Mr,Mris,stdPhsSet,permtrSet,risGeniePhsSet);
        nnCostTimeSet(mm) = nnCostTimeSet(mm) + nnCostTime; nnRmseSet(mm) = nnRmseSet(mm) + nnRmse;
        rcgLsCostTimeSet(mm) = rcgLsCostTimeSet(mm) + rcgLsCostTime; rcgLsRmseSet(mm) = rcgLsRmseSet(mm) + rcgLsRmse;
    end
end
rcgLsCostTimeSetAve = rcgLsCostTimeSet/nn; rcgLsRmseSetAve = rcgLsRmseSet/nn;

%% plot
plt1 = {'ro','kx','b^'};
plt2 = {'r-','k--','b:'};
figure;

xlim([8,64])
plot(MrisSet,nnCostTimeSetAve,'r-x','linewidth',2,'markersize',8);
hold on
plot(MrisSet,rcgLsCostTimeSetAve,'k-o','linewidth',2,'markersize',8);
grid on
legend('QNN','Algorithm Proposed in [10]');
xlabel('Number of RIS Elements');
ylabel('Running Time (s)');
xlim([min(MrisSet) max(MrisSet)])








