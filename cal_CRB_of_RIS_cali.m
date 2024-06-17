function fisher = cal_CRB_of_RIS_cali(permtrSet,risGeniePhsSet,Hbreve)
Mr = size(Hbreve,1);
[O,Mris,Nb,~] = size(permtrSet);
mrisNb = Mris*Nb;
mrNb = Mr*Nb;
tmp1 = risGeniePhsSet.';
g = tmp1(:);
tmp2 = kron(Hbreve, eye(Nb));

tolNum1 = Mris*(Nb-1);
tolNum2 = Mr*Mris;
tolNum = tolNum1+2*tolNum2;

fisher = 0;
for oo = 1 : O
    %% muDerPhi
    muDerEta = zeros(mrNb,tolNum);
    permtrSeto = squeeze(permtrSet(oo,:,:,:));
    muDerEta(:,1:tolNum1) = cal_muDerPhi();
    %% muDerHbreve
    muDerEta(:,tolNum1+(1:tolNum2)) = cal_muDerh();
    muDerEta(:,end-tolNum2+1:end) = 1j*muDerEta(:,tolNum1+(1:tolNum2));
    fisher = fisher+muDerEta'*muDerEta;
end
fisher = real(fisher);
    function y = cal_muDerPhi()
        y = zeros(mrNb,Mris*(Nb-1));
        Pio = zeros(mrisNb);
        for mm = 1:Mris
            Pio((mm-1)*Nb+(1:Nb),(mm-1)*Nb+(1:Nb)) = squeeze(permtrSeto(mm,:,:));
        end
        tmp3 = 1j*(tmp2*Pio)*diag(g);
        for mm = 1 : Mris
            y(:,(mm-1)*(Nb-1)+(1:Nb-1)) = tmp3(:,(mm-1)*Nb+(2:Nb));
        end
    end
    function y = cal_muDerh()
        y = zeros(mrNb,tolNum2);
        for ii = 1 : Mris
            for jj = 1 : Mr
                Pioj = squeeze(permtrSeto(ii,:,:));
                gj = tmp1(:,ii);
                y((jj-1)*Nb+(1:Nb),(jj-1)*Mris+ii) = Pioj*gj;
            end
        end
    end
end