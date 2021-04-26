% Read raw data

[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');


% Define Constants

use_complex_tss=1;
cll_volume=166;

% number of tissue measurements that shall be taken into account when
% computing the measurements
tss_measurements=[3];
bm_measurements=[0];
refill=1;
seperate_sd=0;

patients=[30];
bestModels=zeros(length(patients),7);

folder="/Figures/11_13/";

AllNums=[];

counter =1;
% Get measurement data
for patient = patients
    
    [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt, chars, charDesc, patient);
    [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements);
    [tspanZAvailable, zdataAvailable, tspanZ,zdata,normZdata]=getBMData(patient, startDate, bm_measurements);

    % 
    %get all params and measures
    [num_x_y_z, txt_x_y_z,~] = xlsread(strcat(pwd,['/Figures/11_13/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
      
    
   % allNums=[num_x_y_z;num_z_y_x;num_x_z_y;num_z_x_y;num_y_z_x;num_y_x_z];
    allNums=num_x_y_z;
    AllNums=[AllNums;allNums];
    txt_x_y_z=txt_x_y_z(1,:);
    [~, model]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'model')));
    model=model(1);
    
    
    [numbOfModels,~]=size(allNums);
    measure=zeros(numbOfModels,5);
    
    
    for i = 1:numbOfModels
        modelId=allNums(i,model);
        all_pars=allNums(i,3:14);
        pars=all_pars(getRelevantParams(modelId));

        [sd_x, sd_y, sd_z,~]=compute_add_sd(pars,xdata,ydata,zdata,tspanX, tspanY, tspanZ, modelId, refill,seperate_sd,[],[],[]);
%         [sd_x, sd_y, sd_z,sd_hb]=compute_add_sd(pars,xdata,ydata,zdata,tspanX, tspanY, tspanZ, modelId, refill,seperate_sd,HB',tspanYPlats', hb_model);
        likelihood=likelihood_add_noise(pars,xdata,ydata,zdata,tspanX, tspanY,tspanZ, modelId,refill ,seperate_sd,[],[],[]);

        bic=(length(pars)+1+seperate_sd*2)*log(length([xdata, ydata, zdata]))-2*log(likelihood);
        aicr=2*(length(pars)+1+seperate_sd*2)-2*log(likelihood);
        
        measure(i,:)=[modelId, patient, bic, aicr, likelihood];
    end   
        
        bestModels(counter,1)=patient;
        bestModels(counter,2)=measure(measure(:,3)==min(measure(:,3)),1);
        bestModels(counter,3)=min(measure(:,3)); 
        bestModels(counter,4)=measure(measure(:,4)==min(measure(:,4)),1);
        bestModels(counter,5)=min(measure(:,4));
        bestModels(counter,6)=measure(measure(:,5)==max(measure(:,5)),1);
        bestModels(counter,7)=max(measure(:,5)); 
        
%         writetable(array2table(measure, 'VariableNames',["model","patient", "bic", "aicr","likelihood"]),strcat(pwd,folder,'Params_Pat',int2str(patient),'.xlsx'));
        counter=counter+1;
end 

varNames=getVarNames(1003);
writetable(array2table(bestModels, 'VariableNames',["patient","bestBic", "bic","bestAic", "aicr","bestLikelihood","likelihood"]),strcat(pwd,folder,'bestModels.xlsx'));
writetable(array2table(AllNums, 'VariableNames',["model","ACC",varNames,"sdx","sdy","sdz", "bic", "aicr","likelihood"]),strcat(pwd,folder,'params.xlsx'));




