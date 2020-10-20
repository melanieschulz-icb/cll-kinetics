function [featX, featY, featZ] = normedDifferencesMultiplicativeNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill)
%Calculates the normed differences between the data and simulation output
%normalizes different datasources with l2 norm
%outputs:                          
%   [featX, featY, featZ] = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for measurements
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y and z respectively

% determine differences
% [featX, featY, featZ] = differencesMultiplicativeNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill);
% 
% 
% %normalize
% featX = featX./norm(log(xdata));
% featY = featY./norm(log(ydata)); 
% featZ = featZ./norm(log(zdata));

[x_out, y_out, z_out] = solution(par, tspanX, tspanY, tspanZ, model);

%calculate total number of cells within the BM
BM_cells=getBMTotalCount(zdata, z_out(1), refill);
%calculate CLL percentage on all cells
z_percentage=z_out./BM_cells;

%log data
x_out_logged=log(x_out/norm(xdata));
y_out_logged=log(y_out/norm(ydata));
z_out_logged=log(z_percentage/norm(zdata));

xdata_logged=log(xdata/norm(xdata));
ydata_logged=log(ydata/norm(ydata));
zdata_logged=log(zdata/norm(zdata));

 
%calculate differences 
featX = (xdata_logged - x_out_logged);
featY = (ydata_logged - y_out_logged); 
featZ=(zdata_logged - z_out_logged);