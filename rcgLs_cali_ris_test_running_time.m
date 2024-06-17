function [nnCostTime,rmse] = rcgLs_cali_ris_test_running_time(Heff,O,Nb,Mr,Mris,stdPhsSet,permtrSet,risGeniePhsSet)
Yvec = gen_hvec(Heff,O,Nb,Mr);
%% calibration
risStdPhsSet = ones(Mris,1)*stdPhsSet;
risEstPhsSet = risStdPhsSet;
iterNum = 1e3;
eps1 = 1e-2;
tloss2 = Inf;
tic;
for ii = 1 : iterNum
    risEstEmplyed = gen_multi_d(risEstPhsSet,permtrSet);
    %% channel estimation
    HEst = wls_chan_est_fast_version(risEstEmplyed,Heff);
    E = gen_E(HEst,permtrSet);
    %% phase estimation
    F = risEstPhsSet.';
    f = F(:);
    [f, cost] = my_manopt_method_correlate_chan(Yvec,E,f,Mr,Nb,O);
    F = reshape(f,Nb,Mris);
    risEstPhsSet = F.';
    
    risEstPhsSetRmAm = remove_ambiguity(risEstPhsSet);
    mse = norm(angle(risEstPhsSetRmAm./risGeniePhsSet),'fro')^2/(Mris*(Nb-1));
    rmseTmp = sqrt(mse)*180/pi;
    if abs(tloss2-rmseTmp)/rmseTmp < eps1
        rmse = rmseTmp;
        nnCostTime = toc;
        break;
    end
    tloss2 = rmseTmp;
end
% figure;semilogy(testData)
end