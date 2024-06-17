function [nnCostTime,rmse] = qnn_cali_ris_test_running_time(hLabel,gearSet,risEstEmplyed,estPsSet,risGeniePhsSet)
[~,trainSize] = size(hLabel);
Mris = size(risEstEmplyed,1);
Nb = size(estPsSet,2);
%% NN training
times = 1e5;
learningRate = 0.5e-2;
%% initialization
tloss2 = Inf;
eps1 = 1e-2;
H = hLabel*risEstEmplyed'/(risEstEmplyed*risEstEmplyed');
%% simulation begin
tic;
for tt = 1 : times
    loss = 0;
    test = zeros(1,trainSize);
    for rr = 1:trainSize
        %% feed forward
        g = gearSet(:,rr); % input
        phi = gear2phase(estPsSet,g);
        d = exp(1j*phi);
        h = hLabel(:,rr);
        nnOutput = H*d;
        
        cDerOut = h - nnOutput;
        test(rr) = norm(cDerOut,'fro');
        loss = loss + test(rr)^2;
        %% back propagation
        % H
        dH = -h*d'+H*(d*d');
        % phi
        dPhi = -2*imag((H'*cDerOut).*conj(d));
        % 
        %% update the gradient
        H = H - learningRate*dH;
        phi = phi-learningRate*dPhi;
        estPsSet = update_the_phi(estPsSet,g,phi);
      
    end
    estPhsSet = exp(1j*estPsSet);
    risEstPhsSetRmAm = remove_ambiguity(estPhsSet);
    mse = norm(angle(risEstPhsSetRmAm./risGeniePhsSet),'fro')^2/(Mris*(Nb-1));
    %% test
    rmseTmp = sqrt(mse)*180/pi; 

%     tlossTmp = loss/trainSize;
%     tloss(tt) = tlossTmp;
%     if mod(tt, 2e3) == 0 %|| tt == 1
%         fprintf(['Training Steps: ', num2str(tt), '/', num2str(times), ';\t', 'loss: ', num2str(tlossTmp), '\n']);
%     end
    if mod(tt, 1e1) == 0
        if abs(tloss2-rmseTmp)/rmseTmp < eps1 
            rmse = rmseTmp;
            nnCostTime = toc;
            break;
        end
        tloss2 = rmseTmp;
    end
end
% figure; semilogy(tloss);
end


