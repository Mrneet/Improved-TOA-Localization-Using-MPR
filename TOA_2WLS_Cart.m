function [varargout] = TOA_2WLS_Cart(senPos,r,varargin)
% TOA localization for precise sensors and erronous sensors.

% Z. Ma and K. C. Ho, "TOA localization in the presence of random sensor 
% position errors," in 2011 IEEE International Conference on Acoustics, 
% Speech and Signal Processing (ICASSP), 2011, pp. 2468-2471.

% by Y. Sun, UESTC, 2018.12.15

%   varargout: position in Cartesian (and Modified Polar Representation)
%   senPos: sensor positions
%   r: TOA measurements
%   varargin: covariance matrix of TOA (and covariance matrix of sensor errors)


[N,M] = size(senPos);
L = length(varargin);

h1 = r.^2 - diag(senPos'*senPos);
G1 = -2*[senPos',-0.5*ones(M,1)];
psi1 = (G1'*G1)\G1'*h1;

V1 = 2*diag(r);
Qr = varargin{1};
if L == 1
    W1 = inv(V1*Qr*V1');
elseif L == 2
    Qs = varargin{2};
    for s = 1:M
        O1(s,(1:N)+(s-1)*N) = (psi1(1:N)-senPos(:,s))';
    end
    W1 = inv(V1*Qr*V1'+O1*Qs*O1');
else
    error('too many input');
end
psi2 = (G1'*W1*G1)\G1'*W1*h1;
C1 = inv(G1'*W1*G1);

h2 = [psi2(1:N).^2;psi2(4)];
G2 = [eye(N);ones(1,N)];
V2 = diag([2*psi2(1:N);1]);
W2 = inv(V2*C1*V2');
psi3 = (G2'*W2*G2)\G2'*W2*h2;

u1e = sign(psi2(1:N)).*sqrt(psi3);
u1ef = sign(real(u1e)).*abs(u1e);
varargout{1} = u1ef;
varargout{2} = [atan2(u1ef(2)-senPos(2,1),u1ef(1)-senPos(1,1));atan2(u1ef(3)-senPos(2,1),norm(u1ef(1:2)-senPos(1:2,1),2));1/norm(u1ef-senPos(:,1),2)];

end

