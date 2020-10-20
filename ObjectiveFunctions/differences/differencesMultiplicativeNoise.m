function [featX, featY, featZ] = differencesMultiplicativeNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill)
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

%log data
x_out_logged=log(x_out);
y_out_logged=log(y_out);
z_out_logged=log(z_percentage);

xdata_logged=log(xdata);
ydata_logged=log(ydata);
zdata_logged=log(zdata);

 
%calculate differences 
featX = (xdata_logged - x_out_logged);
featY = (ydata_logged - y_out_logged); 
featZ=(zdata_logged - z_out_logged);


end