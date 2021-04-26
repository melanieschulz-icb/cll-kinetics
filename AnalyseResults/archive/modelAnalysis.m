clear all;

patients = [1:11,13:30];

compMeasure='bic';
goodThreshold=2;
rejThreshold=10;

GoodModels=[];
BestModel=[];

RelevantParams=[];
RejectedParams=[];
names=getVarNames(3);
for patient = patients
    [num_x_y_z, txt_x_y_z,~] = xlsread(strcat(pwd,['/Figures/10_30/x_y_z/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
    [num_z_y_x, txt_z_y_x,~]=xlsread(strcat(pwd,['/Figures/10_30/z_y_x/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
    [num_x_z_y, txt_x_z_y,~]=xlsread(strcat(pwd,['/Figures/10_30/x_z_y/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
    [num_z_x_y, txt_z_x_y,~]=xlsread(strcat(pwd,['/Figures/10_30/z_x_y/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
    [num_y_z_x, txt_y_z_x,~]=xlsread(strcat(pwd,['/Figures/10_30/y_z_x/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
    [num_y_x_z, txt_y_x_z,~]=xlsread(strcat(pwd,['/Figures/10_30/y_x_z/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));

    
    
    allNums=[num_x_y_z;num_z_y_x;num_x_z_y;num_z_x_y;num_y_z_x;num_y_x_z];
    txt=txt_x_y_z(1,:);
    [~, model]=find(~cellfun(@isempty, regexp(txt, 'model')));
    model=model(1);
    [~, d_PB]=find(~cellfun(@isempty, regexp(txt, 'd_PB')));
    [~, d_LT]=find(~cellfun(@isempty, regexp(txt, 'd_LT')));
    [~, m_LT_PB]=find(~cellfun(@isempty, regexp(txt, 'm_LT_PB')));
    [~, m_PB_LT]=find(~cellfun(@isempty, regexp(txt, 'm_PB_LT')));
    [~, c_LT]=find(~cellfun(@isempty, regexp(txt, 'c_LT')));
    [~, d_BM]=find(~cellfun(@isempty, regexp(txt, 'd_BM')));
    [~, m_BM_PB]=find(~cellfun(@isempty, regexp(txt, 'm_BM_PB')));
    [~, c_PB]=find(~cellfun(@isempty, regexp(txt, 'c_PB')));
    [~, c_BM]=find(~cellfun(@isempty, regexp(txt, 'c_BM')));
    [~, m_LT_BM]=find(~cellfun(@isempty, regexp(txt, 'm_LT_BM')));
    
    [~, measure]=find(~cellfun(@isempty, regexp(txt_x_y_z, compMeasure)));
%     [~, bic]=find(~cellfun(@isempty, regexp(txt, 'bic')));
%     [~, aic]=find(~cellfun(@isempty, regexp(txt, 'aic')));
    [~, likelihood]=find(~cellfun(@isempty, regexp(txt, 'likelihood')));

%     [bestBicFirstModels,i]=min(num(:,bic));
%     [bestBicAddModels, j]=min(num2(:,bic));
    
    [bestMeasure,i]=min(allNums(:,measure));
    bestModel=allNums(i,model);
    
%     if bestBicAddModels<bestBicFirstModels
%         bestmodel=num2(j, model);
%     else
%         bestmodel=num(i,model);
%     end
    
%     bestBic=min(bestBicFirstModels, bestBicAddModels);
    
    BestModel=[BestModel;patient,bestModel,bestMeasure];
%     BestModelFirstModel=[BestModelFirstModel;patient, bestBicFirstModels];
    
    goodModels=unique(allNums(allNums(:,measure)-bestMeasure<goodThreshold))';
    rejectedModels=unique(allNums(allNums(:,measure)-bestMeasure>rejThreshold))';
    nonRejectedModels=allNums(allNums(:,measure)-bestMeasure<=rejThreshold,:);

    modelParams=zeros(length(goodModels),20);

    counter=1;
    for modelId=goodModels(:,model)
        modelParams(counter,1)=patient;
        modelParams(counter,2)=modelId;
        modelParams(counter,getRelevantParams(modelId)+2)=1;
        if isempty(find(~cellfun(@isempty, regexp(getVarNames(modelId), 'm_LT_PB')), 1))
            modelParams(counter,m_LT_PB)=0;
            modelParams(counter,16)=1;
        end
        if isempty(find(~cellfun(@isempty, regexp(getVarNames(modelId), 'm_BM_PB')), 1))
            modelParams(counter,m_BM_PB)=0;
            modelParams(counter,17)=1;
        end
        if isempty(find(~cellfun(@isempty, regexp(getVarNames(modelId), 'm_LT_BM')), 1))
            modelParams(counter,m_LT_BM)=0;
            modelParams(counter,18)=1;
        end
        modelParams(counter,19)=allNums(unique(allNums(:,model)==modelId),measure);
        modelParams(counter,20)=allNums(unique(allNums(:,model)==modelId),likelihood);

        counter=counter+1;
    end
    modelParams=modelParams(:,[1:2,find(cellfun(@isempty, regexp(names, '0')))+2,16:20]);
    GoodModels=[GoodModels; modelParams];

    nonRejectedModels(:,18:23)=nonRejectedModels(:,15:20);
    nonRejectedModels(:,15:17)=0;
    
    counter = 1;
    for modelId=nonRejectedModels(:,model)'

        if isempty(find(~cellfun(@isempty, regexp(getVarNames(modelId), 'm_LT_PB')), 1))
            nonRejectedModels(counter,15)=nonRejectedModels(counter,m_LT_PB);
            nonRejectedModels(counter,m_LT_PB)=0;

        end
        if isempty(find(~cellfun(@isempty, regexp(getVarNames(modelId), 'm_BM_PB')), 1))
            nonRejectedModels(counter,16)=nonRejectedModels(counter,m_BM_PB);
            nonRejectedModels(counter,m_BM_PB)=0;
        end
        if isempty(find(~cellfun(@isempty, regexp(getVarNames(modelId), 'm_LT_BM')), 1))
            nonRejectedModels(counter,17)=nonRejectedModels(counter,m_LT_BM);
            nonRejectedModels(counter,m_LT_BM)=0;
        end

        counter=counter+1;
    end
    
    rejectedParams=zeros(1,12);
    counter=1;
    for param = [d_PB,d_LT,m_LT_PB,c_LT,d_BM,m_BM_PB,c_BM,m_LT_BM,c_PB,15:17]
       if isempty(find(nonRejectedModels(:,param)>0, 1))
           rejectedParams(counter)=1;
       end
       counter=counter+1;
    end  
    
    RejectedParams=[RejectedParams; patient,rejectedParams];
    


    
    relevantParams=zeros(1,12);
    counter=1;
    for param = [d_PB,d_LT,m_LT_PB,c_LT,d_BM,m_BM_PB,c_BM,m_LT_BM,c_PB,15:17]
       if isempty(find(nonRejectedModels(:,param)==0, 1))
           relevantParams(counter)=1;
       end
       counter=counter+1;
    end  
    
    RelevantParams=[RelevantParams; patient,relevantParams];
%     
%     %mxz relevant?
%     if isempty(find(nonRejectedModels(:,m_LT_BM)==0 | nonRejectedModels(:,model)>=100,1))
%         relevantParams(counter)=1;
%     end
%     counter=counter+1;
%     
%     %mzx relevant?
%     if isempty(find(nonRejectedModels(:,m_BM_LT)==0 | nonRejectedModels(:,model)<100,1))
%         relevantParams(counter)=1;
%     end
%     
%     RelevantParams=[RelevantParams; patient,relevantParams];
%     complexModelRejected=min(rejectedModels)==3;
%     
%     
%     %mxz rejected?
%     if isempty(find(nonRejectedModels(:,m3)>0 & nonRejectedModels(:,model)<100,1))
%         rejectedParams(counter)=1;
%     end
%     counter=counter+1;
%     
%     %mzx rejected?
%     if isempty(find(nonRejectedModels(:,m3)>0 & nonRejectedModels(:,model)>=100,1))
%         rejectedParams(counter)=1;
%     end
%     
%     RejectedParams=[RejectedParams; patient,rejectedParams];
    
end