function [featX, featY, featZ] = differencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill)
%Calculates the normed differences between the data and simulation output
%normalizes different datasources with l2 norm
%outputs:                          
%   [featX, featY, featZ] = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for measurements
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y and z respectively

% determine output for adapted parameters
[x_out, y_out, z_out] = solution(par, tspanX, tspanY, tspanZ, model);

%calculate total number of cells within the BM
BM_cells=getBMTotalCount(zdata, z_out(1), refill);

%calculate CLL percentage on all cells
z_percentage=z_out./BM_cells;

%calculate differences 
featX = (xdata - x_out);
featY = (ydata - y_out); 
featZ = (zdata - z_percentage);

end