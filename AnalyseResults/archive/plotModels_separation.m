% Read raw data

[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');


% Define Constants

tspan = 0:1:500;
tspanx=0:1:400;
tspanz=0:1:700;
use_complex_tss=1;
cll_volume=166;
tss_measurements=3;
counter=0;

compMeasure='bic';
goodThreshold=2;
posThreshold=6;
strongThreshold=10;


basefolder="/Figures/12_09ui/";

colors=['m','c','r','g','b','w','k','y'];
          
% 
% Get measuremen data
for patient = [1:11,13:30]
    counter =0;
    
    [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt, chars, charDesc, patient);
    [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements);
    [tspanZ, zdata, ~]=getBMData(patient, startDate);
    
    % 
    %get all params and measures
    [allNums, txt_x_y_z,~] = xlsread(strcat(pwd,basefolder,['patient_',int2str(patient),'/Params_Pat',int2str(patient)])); 
    

    txt_x_y_z=txt_x_y_z(1,:);
    [~, model]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'model')));
    model=model(1);

    [~, measure]=find(~cellfun(@isempty, regexp(txt_x_y_z, compMeasure)));
    [~, z0]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'BM_0')));
    [~, aic]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'aic')));
    [~, likelihood]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'likelihood')));

    
    %find the best model
    [bestMeasure,i]=min(allNums(:,measure));
    bestModel=allNums(i,model);
    
    
    
    %models with measure difference <goodThreshold
    goodModels=allNums(allNums(:,measure)-bestMeasure<=goodThreshold,:);
    [~,i]=unique(goodModels(:,model));
    goodModels=goodModels(i,:);

    
    %models with big difference <posThreshold
    posRejectedModels=allNums(allNums(:,measure)-bestMeasure<=posThreshold&allNums(:,measure)-bestMeasure>goodThreshold,:);
    [~,i]=unique(posRejectedModels(:,model));
    posRejectedModels=posRejectedModels(i,:);
    
    
    %models with big difference <strongThreshold
    stronglyRejectedModels=allNums(allNums(:,measure)-bestMeasure<=strongThreshold&allNums(:,measure)-bestMeasure>posThreshold,:);
    [~,i]=unique(stronglyRejectedModels(:,model));
    stronglyRejectedModels=stronglyRejectedModels(i,:);
    
    %models with big difference >strongThreshold
    rejectedModels=allNums(allNums(:,measure)-bestMeasure>strongThreshold,:);
    [~,i]=unique(rejectedModels(:,model));
    rejectedModels=rejectedModels(i,:);

    % Create Data and Plot
    f=figure;
    counter=counter+1;
    subplot(1,3,counter);
    
    [goodModelsCount,~]=size(goodModels);
    
    X_Solutions=zeros(goodModelsCount,length(tspanx));
    Y_Solutions=zeros(goodModelsCount,length(tspan));
    Z_Solutions=zeros(goodModelsCount, length(tspanz));
    
    
    for i = 1:goodModelsCount
        all_pars=goodModels(i,3:14);
        pars=all_pars(getRelevantParams(goodModels(i,model)));
        [x, y, z]=solution(pars,tspanx, tspan, tspanz,goodModels(i,model));
        X_Solutions(i,:)=x;
        Y_Solutions(i,:)=y;
        Z_Solutions(i,:)=z;
    end   
        %outliner not taken into account


    %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
    
    for i=1:goodModelsCount
        plot(tspan, Y_Solutions(i,:),colors(mod(i-1,length(colors))+1),'MarkerSize',10); hold on;
    end

    
    plot(tspanY, ydata, '.k'); hold on;
    plot(ty_outlier, y_outlier, '.r'); 
    %if counter==1
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16)
        ylabel('ALC in blood','FontWeight','bold','FontSize',16)
     if counter==1
        legend([string(goodModels(:,model));"Estimated ALC"],'FontSize',12)
     end

    %end
    
    
    plotTitle='goodModelsPatient'+string(patient);
    title(plotTitle, 'Fontsize',18);
    xlim([0 500])
    ylim([0, max(ydata)+1e11]);
    if patient==14
        set(gca, 'ytick',0:5e11:round(max(ydata)+1e11)) 
    elseif patient==16 || patient==21 || patient==23 || patient ==25
        set(gca, 'ytick',0:2e11:round(max(ydata)+1e11)) 
    else
        set(gca, 'ytick',0:1e11:round(max(ydata)+1e11))
    end
    set(gca,'xtick',0:100:500);
    
    counter=counter+1;
    subplot(1,3,counter)
    
    
    for i=1:goodModelsCount
        plot(tspanx, X_Solutions(i,:),colors(mod(i-1,length(colors))+1),'MarkerSize',10); hold on;
    end
    %tissue measurements taken into account
    plot(tspanXAvailable(1:min(tss_measurements,length(tspanXAvailable))), XdataAvailable(1:min(tss_measurements,length(tspanXAvailable))), '.k','MarkerSize',10); hold on;
    %tissue measurements not taken into account
    plot(tspanXAvailable(min(tss_measurements,length(tspanXAvailable))+1:end), XdataAvailable(min(tss_measurements,length(tspanXAvailable))+1:end), '.r'); hold on;
    %if counter==2
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
        ylabel('ALC in tissue','FontWeight','bold','FontSize',16)
  
    legend([string(goodModels(:,model));"Estimated ALC"],'FontSize',12)
        %end
    ylim([0, max(xdata)+5e11]);
    set(gca,'xtick',0:100:500);
    set(gca, 'ytick',0:1e12:round(max(xdata)+5e11)) 

    
        
    counter=counter+1;
    subplot(1,3,counter);
    
     
    for i=1:goodModelsCount
%         yyaxis(ax,'left');
%         plot(ax,tspanz, Z_Solutions(i,:),'Color',colors(i),'MarkerSize',10); hold on;
% 
%         yyaxis(ax(1),'right');
        %BM_cells=getBMTotalCount(zdata, z_out(1), refill);
        plot(tspanz, Z_Solutions(i,:)./Z_Solutions(i,1),'Color',colors(mod(i-1,length(colors))+1),'MarkerSize',10); hold on;

    end
    
    plot(tspanZ, zdata/zdata(1), '*k'); hold on;
     
     
    xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
    ylabel('ALC in BM','FontWeight','bold','FontSize',16)
     
    legend([string(goodModels(:,model));"Percentage of CLL cells"],'FontSize',12)
    
    
    savefig(f,strcat(pwd,basefolder,plotTitle),'compact');
end 