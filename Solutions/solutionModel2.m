function [x_out, y_out, z_out] = solutionModel2(par, tspanX, tspanY, tspanZ)
%	Solution to Model 2
%Outputs:
%   x/y/z _out        output from integral for respectively x, y and z

%Input:
%   par             parameters for integral
%   tspan X/Y/Z     time span for respectively x, y and z

d2 = par(1);
d1 = par(2);
m = par(3);
x0 = par(4);
y0 = par(5);
c = par(6);
d3=par(7);
m2=par(8);
z0=par(9);
c2=par(10);

A = d1*c;
B = (m+d1);

A2=d3*c2;
B2=(m2+d3);

x_out = A/B + (x0 - A/B)*exp(-B*tspanX);

z_out=A2/B2 + (z0 - A2/B2)*exp(-B2*tspanZ);

y_out = (m*A)/(B*d2) + (m2*A2)/(B2*d2) +...
   ( m/(d2 - B) * (x0 - A/B)*exp(-B*tspanY) ) + ( m2/(d2 - B2) * (z0 - A2/B2)*exp(-B2*tspanY) )+...
    (y0 - (m*A)/(B*d2) - (m/(d2 - B) * (x0 - A/B)) - (m2*A2)/(B2*d2) - (m2/(d2 - B2) * (z0 - A2/B2)) )*exp(-d2 * tspanY);

end