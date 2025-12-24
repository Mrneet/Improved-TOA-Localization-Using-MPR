function [mprSol,carSol] = TOA_MPR_MLE(iniMPR,senPos,rNsed,Q)
% MLE for TOA-based localization method in MPR
%   Detailed explanation goes here
warning off

for it = 1:30
    u1 = [cos(iniMPR(1))*cos(iniMPR(2)); sin(iniMPR(1))*cos(iniMPR(2)); sin(iniMPR(2))];
    k1 = u1 - senPos*iniMPR(end);
    K1 = k1./sqrt(sum(k1.^2,1));
    
    K2 = 1/iniMPR(3)* ...
        [-sin(iniMPR(1))*cos(iniMPR(2)),  -cos(iniMPR(1))*sin(iniMPR(2)),   -cos(iniMPR(1))*cos(iniMPR(2))/iniMPR(3);
        cos(iniMPR(1))*cos(iniMPR(2)),  -sin(iniMPR(1))*sin(iniMPR(2)),   -sin(iniMPR(1))*cos(iniMPR(2))/iniMPR(3);
        0,                               cos(iniMPR(2)),                  -sin(iniMPR(2))/iniMPR(3)];
    K = K1'*K2;
    
    r_bar = sqrt(sum((u1-senPos*iniMPR(3)).^2,1))'/iniMPR(3);
    itMPR2 = iniMPR + (K'/Q*K)\K'/Q*(rNsed-r_bar);
    %             if norm(itMPR1-itMPR2)<1e-10;break;end
    iniMPR = itMPR2;
end
mprSol = iniMPR;
carSol = [cos(iniMPR(1))*cos(iniMPR(2)); sin(iniMPR(1))*cos(iniMPR(2)); sin(iniMPR(2))]/iniMPR(3);
            
end

