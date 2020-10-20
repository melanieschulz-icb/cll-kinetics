function [x_out, y_out, z_out] = solutionModel12(par, tspanX, tspanY, tspanZ)
% Integraded function, integration established by Lea
%Outputs:
%   time            time span used in the integral
%   x/y/z _out        output from integral for respectively x, y and z

%Input:
%   par             parameters for integral
%   tspan X/Y?Z       time span for respectively x, y and z

d2 = par(1);

m = par(2);
x0 = par(3);
y0 = par(4);


m2=par(5);
z0=par(6);
c2=par(7);
m3=par(8);


B = (m+m3);

A2=c2;
B2=m2;


CxB1=x0;

Cz=A2/B2;
CzB1=m3*CxB1/(B2-B);
CzB2=z0-Cz-CzB1;

Cy=(m2*Cz)/d2;
CyB1=(m*CxB1+m2*CzB1)/(d2-B);
CyB2=m2*CzB2/(d2-B2);
Cyd2=y0-Cy-CyB1-CyB2;

x_out = CxB1*exp(-B*tspanX);

z_out=Cz + CzB1*exp(-B*tspanZ)+CzB2*exp(-B2*tspanZ);

y_out = Cy+CyB1*exp(-B*tspanY)+CyB2*exp(-B2*tspanY)+Cyd2*exp(-d2*tspanY);
end