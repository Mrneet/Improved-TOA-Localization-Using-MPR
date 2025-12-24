function [CRB_MPR,CRB_Cart] = CRLB_TOA(senPos,souLoc,Q)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

theta = atan2(souLoc(2), souLoc(1));
phi = atan2(souLoc(3), norm(souLoc(1:2),'fro'));
u0 = [cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)];
r = norm(souLoc);

p = u0 - senPos/r;
P = p./sqrt(sum(p.^2,1));

MM = r* ...
    [-sin(theta)*cos(phi),  -cos(theta)*sin(phi),   -cos(theta)*cos(phi)*r;
    cos(theta)*cos(phi),   -sin(theta)*sin(phi),   -sin(theta)*cos(phi)*r;
    0,                     cos(phi),               -sin(phi)*r];

U = P'*MM;
CRB_MPR = inv(U'/Q*U);

U2 = P';
CRB_Cart = inv(U2'/Q*U2);
        
end

