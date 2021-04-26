function [pb_0] = getPB_0(patient,params,model_pre)
    
    load('heavywater-all-counts-decoded.mat');
    [startDate, tspanY, ~, ~, ~] = aggBloodData(num, txt, patient,true);
    
    y_out=solution_PB_pre(params,[tspanY(1),0],model_pre);
    pb_0=y_out(end);

end