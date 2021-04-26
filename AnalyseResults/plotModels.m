% Read raw data

[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');


% Define Constants

tspan = [-100:0.002:100,101:1:700];
tspanx=[0:0.002:100,101:1:700];
%tspanz=[0:0.002:100,101:1:700];
tspanz=[];
numb_compartments=2;
use_complex_tss=1;
cll_volume=166;
tss_measurements=3;
counter=0;
patients=[1:11,13:30];

refill=1;

compMeasure='bic';
goodThreshold=2;
posThreshold=6;
strongThreshold=10;

%If you want to plot a particular model
modelToPlot=[0];
modelPre=[1];

%basefolders=["/Figures/02_03/","/Figures/11_18/","/Figures/02_12/"];
%basefolders=["/Figures/02_24/"];
basefolders=["/Figures/03_22/"];
duplModelDefiner=[1003,2003];
%duplModelDefiner=[3,1003];
%bestModels={376,1047};
bestModels={1};
outputFolder="/Figures/03_22/plots_0/";
titlename="";
sheets=[1];

colors=['m','c','r','g','b','w','k','y'];
colors=[[0.00,0.45,0.74];[1.00,0.41,0.16];[0.47,0.67,0.19];...
    [1.00,1.00,0.07];[0.49,0.18,0.56];[0.50,0.50,0.50]];
bestColor=[0.64,0.08,0.18];
          
% 
% Get measuremen data
for patient = patients
    counter =0;
    
    [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt,patient,true);
    [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements);
    [tspanZ, zdata, ~]=getBMData(patient, startDate, 0);
    
    ydata_relative=turnToRelativeCounts(ydata,patient);
    y_outlier_relative=turnToRelativeCounts(y_outlier,patient);
    
    count=1;
    goodModelsAll=[];
    parTables{1,2}=[];
    
    for basefolder=basefolders
    % 
    %get all params and measures
        [num_x_y_z, txt_x_y_z,~] = xlsread(strcat(pwd,basefolder,['patient_',int2str(patient),'/Params_Pat',int2str(patient)]),sheets(count));
        parTable=readtable(strcat(pwd,basefolder,['patient_',int2str(patient),'/Params_Pat',int2str(patient)]),"Sheet",sheets(count));
        parTables(1,count)={parTable};

        allNums=num_x_y_z;
        txt_x_y_z=txt_x_y_z(1,:);
        [~, model]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'model')));
        model=model(1);

        [~, measure]=find(~cellfun(@isempty, regexp(txt_x_y_z, compMeasure)));
        [~, z0]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'BM_0')));
        [~, aic]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'aic')));
        [~, likelihood]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'likelihood')));

        measure=measure(1);

        %find the best model
        [bestMeasure,i]=min(allNums(:,measure));
        bestModel=allNums(i,model);

        firstModel=allNums(allNums(:,model)==3,:);
        [~,i]=unique(firstModel(:,model));
        firstModel=firstModel(i,:);



        %models with measure difference <goodThreshold
        goodModels=allNums(allNums(:,measure)-bestMeasure<=goodThreshold,:);
        [~,i]=unique(goodModels(:,model));
        goodModels=goodModels(i,:);

        duplModels=getDuplicateModels(duplModelDefiner(count));
        for dupl = duplModels.Variables'
            goodModels=goodModels(goodModels(:,model)~=dupl,:);
        end
        
        for bestModel=bestModels{1,count}
            goodModels=[goodModels;allNums(allNums(:,model)==bestModel,:)];
        end
        [~,i]=unique(goodModels(:,model));
        goodModels=goodModels(i,:);
        

        goodModelsAll=[goodModelsAll;[goodModels(:,model),count*ones(length(goodModels(:,model)),1)]];
        count=count+1;
    end    
    %models with big difference <posThreshold
%     posRejectedModels=allNums(allNums(:,measure)-bestMeasure<=posThreshold&allNums(:,measure)-bestMeasure>goodThreshold,:);
%     [~,i]=unique(posRejectedModels(:,model));
%     posRejectedModels=posRejectedModels(i,:);
%     
%     
%     %models with big difference <strongThreshold
%     stronglyRejectedModels=allNums(allNums(:,measure)-bestMeasure<=strongThreshold&allNums(:,measure)-bestMeasure>posThreshold,:);
%     [~,i]=unique(stronglyRejectedModels(:,model));
%     stronglyRejectedModels=stronglyRejectedModels(i,:);
%     
%     %models with big difference >strongThreshold
%     rejectedModels=allNums(allNums(:,measure)-bestMeasure>strongThreshold,:);
%     [~,i]=unique(rejectedModels(:,model));
%     rejectedModels=rejectedModels(i,:);
    %
    goodModels=goodModelsAll;
    if ~isempty(modelToPlot)
        goodModels=[modelToPlot,1];
    end
%    goodModels=[[2006,1];[2016,1]];
    % Create Data and Plot
    f=figure;
    counter=counter+1;
    subplot(numb_compartments,1,counter);
    
    [goodModelsCount,~]=size(goodModels);
    
    X_Solutions=zeros(goodModelsCount,length(tspanx));
    Y_Solutions=zeros(goodModelsCount,length(tspan));
    Y_Solutions_relative=zeros(goodModelsCount,length(tspan));
    Z_Solutions=zeros(goodModelsCount, length(tspanz));
    bics=zeros(goodModelsCount,1);
    
    
    for i = 1:goodModelsCount
        modelNr=goodModels(i,1);
        folderCounter=goodModels(i,2);
        parTable=parTables{1,folderCounter};
        preModel=parTable(parTable.model_pre==modelPre,:);
        modelPars=preModel(preModel.model==modelNr,:);
        nam=getVarNames_PB_pre(modelPre);
        nam(nam=="groth")="g_all_PB";

        bestRun=modelPars(modelPars.bic==min(modelPars.bic),:);
        prePars=bestRun{:,nam};
        bestRun.PB_0=getPB_0(patient,prePars,modelPre);
        bic=bestRun{:,'bic'};
        all_pars=bestRun{:,getVarNames(modelNr)};
%        pars=zeros(1,12);
%        pars(getRelevantParams(modelNr))=all_pars;
        %pars=all_pars;
        
        
        [x, y, z]=solution(all_pars(getRelevantParams(modelNr)),tspanx, tspan(tspan>0), tspanz,modelNr);
        X_Solutions(i,:)=x;
        y_pre=solution_PB_pre([getPB_0(patient,prePars,modelPre),prePars(2)],tspan(tspan<=0),modelPre+3);
        Y_Solutions(i,:)=[y_pre,y];
        Y_Solutions_relative(i,:)=turnToRelativeCounts([y_pre,y],patient);
        Z_Solutions(i,:)=z;
        bics(i,:)=bic;
    end   
        %outliner not taken into account


    %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
    
    for i=1:goodModelsCount
        if find([bestModels{:}] == goodModels(i,model))>0
            linestyle='-';
            color=bestColor;
        else
            linestyle='--';
            color=colors(mod(i-1,length(colors))+1,:);
        end
        plot(tspan, Y_Solutions(i,:),'color',color,'MarkerSize',10,'LineWidth',1.5,'LineStyle',linestyle); hold on;
    end

    
    plot(tspanY, ydata,  '.k','MarkerSize',12); hold on;
    plot(ty_outlier, y_outlier,  '.r','MarkerSize',12); 
    %if counter==1
        xlabel('days after therapy start', 'FontWeight','bold','FontSize',13)
        ylabel('TLC, blood','FontWeight','bold','FontSize',14)
     if counter==1
%        legend([string(goodModels(:,model))+[", bic: "]+string(bics);"Total lymphocyte counts (est.)"],'FontSize',12)
         legend(["Total lymphocyte counts, blood - model fit";"Total lymphocyte counts (est.)"],'FontSize',12)
     end

    %end
    
    

    xlim([-150 700])
    ylim([0, max(ydata)+1e11]);
    if patient==14
        set(gca, 'ytick',0:5e11:round(max(ydata)+1e11)) 
    elseif patient==16 || patient==21 || patient==23 || patient ==25
        set(gca, 'ytick',0:2e11:round(max(ydata)+1e11)) 
    else
        set(gca, 'ytick',0:1e11:round(max(ydata)+1e11))
    end
    set(gca,'xtick',0:100:700);
    
    plotTitle=strcat(titlename,' Patient '+string(patient));
    title(plotTitle, 'Fontsize',16);
    
    counter=counter+1;
    subplot(numb_compartments,1,counter)
    
    
    for i=1:goodModelsCount
        if find([bestModels{:}] == goodModels(i,model))>0
            linestyle='-';
            color=bestColor;
        else
            linestyle='--';
            color=colors(mod(i-1,length(colors))+1,:);
        end
        plot(tspanx, X_Solutions(i,:),'color',color,'MarkerSize',10,'LineWidth',1.5,'LineStyle',linestyle); hold on;
    end
    %tissue measurements taken into account
    plot(tspanXAvailable(1:min(tss_measurements,length(tspanXAvailable))), XdataAvailable(1:min(tss_measurements,length(tspanXAvailable))), '.k','MarkerSize',12); hold on;
    %tissue measurements not taken into account
    plot(tspanXAvailable(min(tss_measurements,length(tspanXAvailable))+1:end), XdataAvailable(min(tss_measurements,length(tspanXAvailable))+1:end), '.r','MarkerSize',12); hold on;
    %if counter==2
        xlabel('days after therapy start', 'FontWeight','bold','FontSize',13);
        ylabel('TLC, lymph tissue','FontWeight','bold','FontSize',14)

  
    %legend([string(goodModels(:,model))+[", bic: "]+string(bics);"Total lymphocyte counts (est.)"],'FontSize',12)
    
    legend(["Lymphocyte counts, lymph tissue - model fit";"Total lymphocyte counts, lymph tissue (est.)"],'FontSize',12)
        %end
    xlim([-150 700])
    ylim([0, max(xdata)+5e11]);
    set(gca,'xtick',-150:100:700);
    set(gca, 'ytick',0:1e12:round(max(xdata)+5e11)) 

    
    if ~isempty(tspanz)
    
        counter=counter+1;
        subplot(3,1,counter);


        for i=1:goodModelsCount
            if find([bestModels{:}] == goodModels(i,model))>0
                linestyle='-';
                color=bestColor;
            else
                linestyle='--';
                color=colors(mod(i-1,length(colors))+1,:);
            end
           %yyaxis left;
    %         plot(ax,tspanz, Z_Solutions(i,:),'Color',colors(i),'MarkerSize',10); hold on;
    % 
    %         yyaxis(ax(1),'right');
            BM_cells=getBMTotalCount(zdata,Z_Solutions(i,1), refill);
            plot(tspanz, Z_Solutions(i,:)./BM_cells(1),'Color',color,'MarkerSize',10,'LineWidth',1.5,'LineStyle',linestyle); hold on;
    %         yyaxis right;
    %         plot(tspanz, Z_Solutions(i,:),'Color',color,'MarkerSize',10); hold on;

        end
        %yyaxis left;
        plot(tspanZ, zdata,  '.k','MarkerSize',12); hold on;


        xlabel('days after therapy start', 'FontWeight','bold','FontSize',13);
        ylabel('CLL cells (fraction)','FontWeight','bold','FontSize',14)

    %    legend([string(goodModels(:,model))+[", bic: "]+string(bics);"CLL cells (fraction)"],'FontSize',12)
        legend(["Fraction of CLL cells, BM - model fit";"CLL cells (fraction)"],'FontSize',12)

    end
    f=gcf;
    savefig(f,strcat(pwd,outputFolder,plotTitle,".fig"),'compact');
    close;
    
    f=figure;
    
    for i=1:goodModelsCount
        if find([bestModels{:}] == goodModels(i,model))>0
            linestyle='-';
            color=bestColor;
        else
            linestyle='--';
            color=colors(mod(i-1,length(colors))+1,:);
        end
        plot(tspan, Y_Solutions_relative(i,:),'color',color,'MarkerSize',10,'LineWidth',1.5,'LineStyle',linestyle); hold on;
    end

    
    plot(tspanY, ydata_relative,  '.k','MarkerSize',12); hold on;
    plot(ty_outlier, y_outlier_relative,  '.r','MarkerSize',12); 
    %if counter==1
        xlabel('days after therapy start', 'FontWeight','bold','FontSize',13)
        ylabel('ALC per nanooliter','FontWeight','bold','FontSize',14)

    legend([string(goodModels(:,model))+", bic: "+string(bics);"RLC"],'FontSize',12)


    %end
    
    

    xlim([0 700])
    ylim([0, max(Y_Solutions_relative(i,:))+2]);
    if patient==14
        set(gca, 'ytick',0:10:round(max(ydata_relative)+10)) 
    elseif patient==16 || patient==21 || patient==23 || patient ==25
        set(gca, 'ytick',0:20:round(max(ydata_relative)+10)) 
    else
        set(gca, 'ytick',0:10:round(max(ydata_relative)+10))
    end
    set(gca,'xtick',0:100:700);
    
    plotTitle=strcat(titlename,' Patient '+string(patient));
    title(plotTitle, 'Fontsize',16);
    close;
    
    
end 