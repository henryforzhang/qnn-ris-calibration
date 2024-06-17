function [estPhsSet,tloss] = qnn_cali_ris_conv(hLabel,gearSet,risEstEmplyed,estPsSet)
[~,trainSize] = size(hLabel);
%% NN training
times = 1e5;
learningRate = 0.5e-2;
%% initialization
tloss2 = Inf;
eps1 = 1e-6;
H = hLabel*risEstEmplyed'/(risEstEmplyed*risEstEmplyed');
% H = (randn(Mr,Mris)+randn(Mr,Mris))/sqrt(2);
tloss = zeros(1,times);
%% simulation begin
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
    tlossTmp = loss/trainSize;
    tloss(tt) = tlossTmp;
    if mod(tt, 2e3) == 0 %|| tt == 1
        fprintf(['Training Steps: ', num2str(tt), '/', num2str(times), ';\t', 'loss: ', num2str(tlossTmp), '\n']);
    end
    if mod(tt, 1e2) == 0
        if abs(tloss2-tlossTmp) < eps1 && tt > 2000
            break;
        end
        tloss2 = tlossTmp;
    end
end
estPhsSet = exp(1j*estPsSet);

% figure; semilogy(tloss);
end


