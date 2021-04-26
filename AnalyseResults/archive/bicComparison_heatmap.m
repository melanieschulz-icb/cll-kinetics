clear all;

patients = [1];

GoodModels=[];
BestModel=[];
UnclearModels=[];

RelevantParams=[];
RejectedParams=[];
names=getVarNames(3);

bicComparison=zeros(99,29);
duplModels=[1019,1020,1026:1029,1050:1051,1066,1067,1074:1077,1090:1093,1098,1099];

basefolder="/Figures/11_18/";

for patient = patients
    
    %get all params and bics
    [allNums, txt_x_y_z,~] = xlsread(strcat(pwd,basefolder,['patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
%     [num_z_y_x, txt_z_y_x,~]=xlsread(strcat(pwd,['/Figures/10_30/z_y_x/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
%     [num_x_z_y, txt_x_z_y,~]=xlsread(strcat(pwd,['/Figures/10_30/x_z_y/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
%     [num_z_x_y, txt_z_x_y,~]=xlsread(strcat(pwd,['/Figures/10_30/z_x_y/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
%     [num_y_z_x, txt_y_z_x,~]=xlsread(strcat(pwd,['/Figures/10_30/y_z_x/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
%     [num_y_x_z, txt_y_x_z,~]=xlsread(strcat(pwd,['/Figures/10_30/y_x_z/patient_',int2str(patient),'/Params_Pat',int2str(patient)]));

    
    
%     allNums=[num_x_y_z;num_z_y_x;num_x_z_y;num_z_x_y;num_y_z_x;num_y_x_z];
    txt=txt_x_y_z(1,:);
    [~, model]=find(~cellfun(@isempty, regexp(txt, 'model')));
    model=model(1);

    [~, bic]=find(~cellfun(@isempty, regexp(txt, 'bic')));
    [~, aic]=find(~cellfun(@isempty, regexp(txt, 'aic')));
    [~, likelihood]=find(~cellfun(@isempty, regexp(txt, 'likelihood')));
    
%     for dupl=[1019,1020,1026:1029,1050:1051,1066,1067,1074:1077,1090:1093,1098,1099]
%         allNums=allNums(allNums(:,model)~=dupl,:);
%     end
    
    %find the best model
    [bestBic,i]=min(allNums(:,bic));
    bestModel=allNums(i,model);
    
    
    
    %models with bic difference <2
    goodModels=unique(allNums(allNums(:,bic)-bestBic<=2&allNums(:,bic)-bestBic>0.1));
    bicComparison(goodModels-1002,patient)=1;
    
    %models with big difference <6
    posRejectedModels=unique(allNums(allNums(:,bic)-bestBic<=6&allNums(:,bic)-bestBic>2));
    for model = posRejectedModels'
        if ~isempty(find(goodModels==model, 1)) || ~isempty(find(bestModel==model, 1))
            UnclearModels=[UnclearModels;patient,model];
        else
            bicComparison(model-1002,patient)=2;
        end
    end
    
    %models with big difference <10
    stronglyRejectedModels=unique(allNums(allNums(:,bic)-bestBic<=10&allNums(:,bic)-bestBic>6));
    for model = stronglyRejectedModels'
        if ~isempty(find(goodModels==model, 1)) || ~isempty(find(posRejectedModels==model, 1)) || ~isempty(find(bestModel==model, 1))
            UnclearModels=[UnclearModels;patient,model];
        else
            bicComparison(model-1002,patient)=3;
        end
    end
    
    
    %models with big difference >10
    rejectedModels=unique(allNums(allNums(:,bic)-bestBic>10));
    for model = rejectedModels'
        if ~isempty(find(goodModels==model, 1)) ||~isempty(find(posRejectedModels==model, 1))||~isempty(find(stronglyRejectedModels==model, 1)) || ~isempty(find(bestModel==model, 1))
            UnclearModels=[UnclearModels;patient,model];
        else
            bicComparison(model-1002,patient)=4;
        end
    end
    

    

    %for each patient, create heatmap with all models, color depending on
    %bic difference
    f=clf;
    subplot(12,1,1);
    h=heatmap(3:92,{"Cat_1"},bicComparison(1:90,patient)');
    h.ColorLimits = [0 4];
    h.Colormap=winter(5);
    h.FontSize=6;
    plotTitle="BIC comparison, patient "+int2str(patient);
    title(plotTitle);
    subplot(12,1,2);
    h2=heatmap(93:181,{"Cat_1"},bicComparison(91:179,patient)','YLabel','');
    h2.ColorLimits = [0 4];
    h2.Colormap=winter(5);
    h2.ColorbarVisible = 'off';
    h2.FontSize=6;
    for counter=3:12
        subplot(12,1,counter);
        lines=(180+63*(counter-3)):(180+63*(counter-2));
        h=heatmap(lines+2,{"Cat_"+int2str(floor(counter/2)+mod(counter,2))},bicComparison(lines,patient)','XLabel',[]);
        h.ColorLimits = [0 4];
        h.Colormap=winter(5);
        h.ColorbarVisible = 'off';
        h.FontSize=6;

    end
    plotTitle='modelComparison_patient_'+string(patient);
    writetable(array2table(UnclearModels, 'VariableNames',["patient","model"]),[pwd,strcat(basefolder,['UnclearModels_Pat',int2str(patient),'.xlsx'])]);
    savefig(f,strcat(pwd,basefolder,plotTitle),'compact');
    close all;
    
end

%for each category, create heatmap with all patients, color depending on
%bic difference
f=clf;
h=heatmap(1:30,3:181,bicComparison(1:179,:));
plotTitle="BIC comparison, models cat_1";
title(plotTitle);
h.ColorLimits = [0 4];
h.Colormap=winter(5);
h.FontSize=6;
savefig(f,strcat(pwd,basefolder,plotTitle),'compact');

f=clf;
h=heatmap(1:30,182:309,bicComparison(180:307,:));
plotTitle="BIC comparison, models cat_2";
title(plotTitle);
h.ColorLimits = [0 4];
h.Colormap=winter(5);
h.FontSize=6;
savefig(f,strcat(pwd,basefolder,plotTitle),'compact');

f=clf;
h=heatmap(1:30,310:437,bicComparison(308:435,:));
plotTitle="BIC comparison, models cat_3";
title(plotTitle);
h.ColorLimits = [0 4];
h.Colormap=winter(5);
h.FontSize=6;
savefig(f,strcat(pwd,basefolder,plotTitle),'compact');

f=clf;
h=heatmap(1:30,438:565,bicComparison(436:563,:));
plotTitle="BIC comparison, models cat_4";
title(plotTitle);
h.ColorLimits = [0 4];
h.Colormap=winter(5);
h.FontSize=6;
savefig(f,strcat(pwd,basefolder,plotTitle),'compact');

f=clf;
h=heatmap(1:30,566:693,bicComparison(564:691,:));
plotTitle="BIC comparison, models cat_5";
title(plotTitle);
h.ColorLimits = [0 4];
h.Colormap=winter(5);
h.FontSize=6;
savefig(f,strcat(pwd,basefolder,plotTitle),'compact');

f=clf;
h=heatmap(1:30,694:821,bicComparison(692:819,:));
plotTitle="BIC comparison, models cat_6";
title(plotTitle);
h.ColorLimits = [0 4];
h.Colormap=winter(5);
h.FontSize=6;
savefig(f,strcat(pwd,basefolder,plotTitle),'compact');


