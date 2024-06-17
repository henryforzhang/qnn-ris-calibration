function [H]=channel_generation_ura(NtSet,NrSet,Nc,Nray,p0)
% ***************************************
%  generate URA channel
%  author - Yimeng Feng
%  input- Nt: transmitter antenna number
%            Nr: receiver antenna number
%            Nc: cluster number
%            Nray: ray number
%            dx_t: inter-element distance for transmitter
%            dx_r: inter-element distance for receiver
%            lambda: wavelength
%            p0: line number for transmitter
%            p1: line number for receiver
%  output-H: URA Channel
%             alpha: multipath gain
%             Ar: antenna response for receiver
%             At: antenna response for transmitter
%copyright - CSRL@Fudan,2021/01/12
%  ************************************
angle_sigma = 20/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
Np = Nc*Nray;
gamma = sqrt(1/Np); %normalization factor
sigma = 1; %according to the normalization condition of the H
AoD = zeros(2,Np);
AoA = zeros(1,Np);
for c = 1:Nc
    AoD_m = unifrnd(0,2*pi,1,2);
    AoA_m = unifrnd(0,2*pi,1,1);
    
    AoD(1,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoD_m(1),angle_sigma);
    AoD(2,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoD_m(2),angle_sigma);
    AoA((c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoA_m(1),angle_sigma);
end
NtSetLen = length(NtSet);
NrSetLen = length(NrSet);
H = cell(NtSetLen,NrSetLen);
for mm = 1 : NtSetLen
    Nt = NtSet(mm);
    for nn = 1 : NrSetLen
        Nr = NrSet(nn);
        Htmp = zeros(Nr,Nt);
        for j = 1:Np
            At = array_response_UPA(AoD(1,j),AoD(2,j),Nt,p0); %UPA array response
            Ar = array_response_ULA(AoA(1,j),Nr);
            alpha = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
            Htmp = Htmp + alpha*Ar*At';
        end
        H{mm,nn} = gamma * Htmp;
    end
end