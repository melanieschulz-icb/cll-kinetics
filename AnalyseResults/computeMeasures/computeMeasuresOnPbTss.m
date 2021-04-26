function par_table = computeMeasuresOnPbTss(par_table,num, txt, chars,charDesc,cll_volume, use_complex_tss, tss_measurements, seperate_sd)
%computeMeasuresOnTwoCompartmemt
%Returns BIC/AIC/likelihood, just looking at the blood and tissue
%compartments


%params, IN:
%   num, txt, chars, charDesc
%           -> to evaluate a patients data
%   cll_volume, use_complex_tss, tss_measurements
%           -> to evaluate the tissue data taken into account
%   seperate_sd
%           -> if to compute seperate or joint sd on blood/tss compartment
%   par_table 
%           -> table for one patient, containing the model information and
%               parameters that belong to the model
%OUT:
%   par_table
%           -> added columns: BIC,AIC,likelihood based on PB, tss. compartmemts

    
    par_table=computeMeasuresOnPbTssBm(par_table,num, txt, chars,charDesc,cll_volume, use_complex_tss, tss_measurements,0,0,seperate_sd);

end

