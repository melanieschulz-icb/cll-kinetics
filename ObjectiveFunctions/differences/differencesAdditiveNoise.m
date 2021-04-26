function [featX, featY, featZ, featHB] = differencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,hbdata, tspanYPlats, hb_model,model_pb_pre)
%Calculates the normed differences between the data and simulation output
%normalizes different datasources with l2 norm
%outputs:                          
%   [featX, featY, featZ] = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for measurements
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y and z respectively

if tspanY(1)<0
    tspanY_pre=tspanY(tspanY<=0);
    tspanY=tspanY(tspanY>0);
else
    tspanY_pre=[];
end

if  ~isempty(tspanY_pre)
    y_pre_out=solution_PB_pre(par(end-1:end),tspanY_pre,model_pb_pre);
else
    y_pre_out=[];
end
% determine output for adapted parameters
if ~isempty([tspanY, tspanX,tspanZ])
    l=length(getRelevantHBParams(hb_model))+length(getRelevantParams_PB_pre(model_pb_pre));
    if ~isempty(tspanY_pre)
        par_nums=getRelevantParams(model);
        ind_pb0=find(par_nums==5,1);
        par=[par(1:ind_pb0-1),y_pre_out(end),par(ind_pb0:end)];
    end
    [x_out, y_out, z_out] = solution(par(1:end-l), tspanX, tspanY, tspanZ, model);
else
    x_out=[];
    y_out=[];
    z_out=[];
end



y_out=[y_pre_out,y_out];

if ~isempty(tspanYPlats)
    hb_out = solutionHB(par, tspanYPlats, model, hb_model);
else
    hb_out=[];
end
%calculate total number of cells within the BM
if ~isempty(z_out)
    if z_out(1)>0
        BM_cells=getBMTotalCount(zdata, z_out(1), refill);
    else
        BM_cells=getBMTotalCount(zdata, max(max(z_out),1), refill);
    end
else
    BM_cells=1;
end

%calculate CLL percentage on all cells
z_percentage=z_out./BM_cells;


%calculate differences 
featX = (xdata - x_out);
featY = (ydata - y_out); 
featZ = (zdata - z_percentage);
featHB=(hbdata-hb_out);

end