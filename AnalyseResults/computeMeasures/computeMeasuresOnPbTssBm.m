function par_table = computeMeasuresOnPbTssBm(par_table,num, txt,cll_volume, use_complex_tss, tss_measurements,bm_measurements,refill ,seperate_sd,columnAppendix)
%computeMeasuresOnThreeCompartmemts
%Returns BIC/AIC/likelihood, just looking at the blood, tissue and bone
%marrow compartments


%params, IN:
%   num, txt, chars, charDesc
%           -> to evaluate a patients data
%   cll_volume, use_complex_tss, tss_measurements
%           -> to evaluate the tissue data taken into account
%   bm_measurement, refill
%           -> to evaluate the bm_measurements taken into account
%   seperate_sd
%           -> if to compute seperate or joint sd on blood/tss compartment
%   par_table 
%           -> table for one patient, containing the model information and
%               parameters that belong to the model
%OUT:
%   par_table
%           -> added columns: BIC,AIC,likelihood based on PB, tss. compartmemts

%     par_table.bic_TwoComp=zeros(size(par_table,1),1);
%     par_table.aicr_TwoComp=zeros(size(par_table,1),1);
%     par_table.like_TwoComp=zeros(size(par_table,1),1);

    patient=unique(par_table.ACC);
    [startDate, tspanY, ydata, ~, ~] = aggBloodData(num, txt, patient);
    [~, ~, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements);
    [~,~, tspanZ, zdata, ~]=getBMData(patient, startDate, bm_measurements);
              
    
        for i=1:size(par_table,1)
            modelId=par_table{i,'model'};
            if any(strcmp(par_table.Properties.VariableNames,'hbmodel'))
                hbmodel=par_table{i,'hbmodel'};
            else
                hbmodel=-1;
            end
            all_pars=par_table{i,getVarNames(modelId)};
            pars=all_pars(getRelevantParams(modelId));
            
            hb_pars=par_table{i,getVarNamesHB(hbmodel)};

            likelihood=likelihood_add_noise([pars,hb_pars],xdata,ydata,zdata,tspanX, tspanY,tspanZ, modelId,refill,seperate_sd,[],[],hbmodel);

            parNumb=length(pars)+1+seperate_sd*3+length(hb_pars);
            dataNumb=length([xdata, ydata,zdata]);
            bic=bicFromLikelihood(parNumb, dataNumb, likelihood);
            aicr=aicFromLikelihood(parNumb,likelihood);

            par_table(i,strcat("bic",columnAppendix))={bic};
            par_table(i,strcat("aic",columnAppendix))={aicr};
            par_table(i,strcat("likelihood",columnAppendix))={likelihood};

        end
end

