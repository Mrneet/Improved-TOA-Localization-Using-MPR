function [mprSol,carSol] = TOA_MPR_SDR(senPos,rNsed,W)
% SDR method for TOA-based localization method in MPR
%   Detailed explanation goes here

[N,M] = size(senPos);

G = [senPos', 0.5*(rNsed.^2-sum(senPos.^2,1)'), -0.5*ones(M,1)];

M1 = blkdiag(eye(N),zeros(2));
% M2 = blkdiag(zeros(N),[0,1;1,0]);
M2 = blkdiag(zeros(N),[0,1;0,0]);
M3 = blkdiag(zeros(N),[0,0;1,0]);

cvx_begin quiet
    variable Vx(5,5) symmetric
    minimize(trace(G'*W*G*Vx))
    subject to
    1-trace(M1*Vx) == 0
%     2-trace(M2*Vx) == 0
    1-trace(M2*Vx) == 0
    1-trace(M3*Vx) == 0
    Vx == semidefinite(5)
cvx_end

[V,D]=eig(Vx);
[lambda,ind]=sort(sum(D));
psi=V(:,ind(end))*sqrt(lambda(end));

mprSol = [atan2(psi(2),psi(1)); atan2(psi(3),norm(psi(1:2))); psi(4)];
carSol = psi(1:N)/psi(4);

if psi(4)<0 || psi(5)<0
    pause;
end

end

