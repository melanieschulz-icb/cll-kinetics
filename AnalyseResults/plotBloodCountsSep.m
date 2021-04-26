% Read raw data

[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');

[numPlat, txtPlat, ~] = xlsread('Heavywater-all counts-decoded.xlsx','Platelets');

[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');




% Define Constants

tspan = 0:1:500;
tspanx=0:1:400;
tspanz=0:1:700;
use_complex_tss=1;
cll_volume=166;
tss_measurements=3;
bm_measurements=5;
counter=0;

compMeasure='bic';
goodThreshold=2;
posThreshold=6;
strongThreshold=10;


basefolder="/Figures/11_18/";

colors=['m','c','r','g','b','w','k','y'];
counter=0;

          
% [9,10,14,15,18,21,28:30];
%[1,2,4,5,7,11,13,16,17,20,22:26]
% Get measuremen data
for patient = [1:2]
    counter =0;
    
    
    [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt, chars, charDesc, patient);
    [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements,false);
    [tspanZ, zdata, ~]=getBMData(patient, startDate,bm_measurements);
    
    [~,tspanYPlats,Platelets,PlateletsScaled, HB,HBScaled, Neutroph,NeutrophScaled, RBC, RBCScaled] = getBloodCounts(numPlat, txtPlat, chars, charDesc,patient);


    
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
    
    for i = 1:goodModelsCount
        % Create Data and Plot for each good model
        counter =0;
        g=figure;
        counter=counter+1;
        subplot(3,3,counter);


        plot(tspan, Y_Solutions(i,:),'c','MarkerSize',10); hold on;

        plot(tspanY, ydata, '.k'); hold on;
        plot(ty_outlier, y_outlier, '.r'); 
        %if counter==1
            xlabel('Time(days)', 'FontWeight','bold','FontSize',16)
            ylabel('ALC in blood','FontWeight','bold','FontSize',16)
         if counter==1
            legend('Best Fit','Estimated ALC','FontSize',12)
         end

        %end



        title('Patient '+ string(patient)+'_Model_'+string(goodModels(i,model)), 'Fontsize',18);
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
        subplot(3,3,counter)

        plot(tspanx, X_Solutions(i,:),'c'); hold on;
    %tissue measurements taken into account
    plot(tspanXAvailable(1:min(tss_measurements,length(tspanXAvailable))), XdataAvailable(1:min(tss_measurements,length(tspanXAvailable))), '.k','MarkerSize',10); hold on;
    %tissue measurements not taken into account
    plot(tspanXAvailable(min(tss_measurements,length(tspanXAvailable))+1:end), XdataAvailable(min(tss_measurements,length(tspanXAvailable))+1:end), '.r'); hold on;
        %if counter==2
            xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
            ylabel('ALC in tissue','FontWeight','bold','FontSize',16)

        legend('Best Fit','estimated ALC','FontSize',12)
            %end
        ylim([0, max(xdata)+5e11]);
        set(gca,'xtick',0:100:500);
        set(gca, 'ytick',0:1e12:round(max(xdata)+5e11)) 


        counter=counter+1;
        subplot(3,3,counter);

        plot(tspanz, Z_Solutions(i,:)/Z_Solutions(i,1), 'c'); hold on;
%         yyaxis right;
%         normal_cells=((1-zdata(1))/zdata(1))*params(startrow,10);
%         percentage=100*z_opt_four./(z_opt_four+normal_cells);
% 
%         plot(tspanz, percentage,'--m',xlabel('days'),ylabel('CLL cell percentage'));hold on;
        plot(tspanZ, zdata/zdata(1), '*k'); hold on;
        %plot(tspanZ, zdatanorm*params(startrow,10), '*k'); hold on;

        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
%         ylabel('ALC in BM','FontWeight','bold','FontSize',16)

        legend('Best Fit','Percentage of CLL cells','FontSize',12)

        counter=counter+1;
        subplot(3,3,counter);


        plot(tspanYPlats',Platelets'-min(Platelets),'*k'); hold on;
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
        ylabel('Platelets, shifted','FontWeight','bold','FontSize',16)
        yyaxis right;
        plot(tspanz, (max(Z_Solutions(i,:))-Z_Solutions(i,:)),'-c'); hold on;
        %plot(tspanZ, max(z_opt_four)-params(startrow,10)*zdatanorm, '*r'); hold on;
        ylabel('negative RLC in BM, shifted','FontWeight','bold','FontSize',12);
        %legend('Platelets development','Reversed BM development - fits','BM percentage','FontSize',10)
        legend('Platelets development','Reversed BM development - fits','FontSize',10)

        counter=counter+1;
        subplot(3,3,counter);


        plot(tspanYPlats',HB'-min(HB),'*k'); hold on;
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
        ylabel('HB, shifted','FontWeight','bold','FontSize',16)
        yyaxis right;
       plot(tspanz, (max(Z_Solutions(i,:))-Z_Solutions(i,:)),'-c'); hold on;
       % plot(tspanZ, max(z_opt_four)-params(startrow,10)*zdatanorm, '*r'); hold on;
    %    ylabel('negative RLC in BM, shifted','FontWeight','bold','FontSize',12);
        %legend('HB development','Reversed BM development - fits','BM percentage','FontSize',10)
     %   legend('HB development','Reversed BM development - fits','FontSize',10)

        counter=counter+1;
        subplot(3,3,counter);


        plot(tspanYPlats',RBC'-min(RBC),'*k'); hold on;
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
        ylabel('RBC, shifted','FontWeight','bold','FontSize',16)
        yyaxis right;
        plot(tspanz, (max(Z_Solutions(i,:))-Z_Solutions(i,:)),'-c'); hold on;
        %plot(tspanZ, max(z_opt_four)-params(startrow,10)*zdatanorm, '*r'); hold on;
        ylabel('negative RLC in BM, shifted','FontWeight','bold','FontSize',12);
        %legend('RBC development','Reversed BM development - fits','BM percentage','FontSize',10)
        legend('RBC development','Reversed BM development - fits','FontSize',10)
        
        
        counter=counter+1;
        subplot(3,3,counter);


        plot(tspanYPlats',Platelets'-min(Platelets),'*k'); hold on;
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
        ylabel('Platelets, shifted','FontWeight','bold','FontSize',16)
        yyaxis right;
        plot(tspan, (max(Y_Solutions(i,:))-Y_Solutions(i,:)),'-c'); hold on;
        %plot(tspanZ, max(z_opt_four)-params(startrow,10)*zdatanorm, '*r'); hold on;
        ylabel('negative ALC in PB, shifted','FontWeight','bold','FontSize',12);
        %legend('Platelets development','Reversed BM development - fits','BM percentage','FontSize',10)
        legend('Platelets development','Reversed PB development - fits','FontSize',10)

        counter=counter+1;
        subplot(3,3,counter);


        plot(tspanYPlats',HB'-min(HB),'*k'); hold on;
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
        ylabel('HB, shifted','FontWeight','bold','FontSize',16)
        yyaxis right;
       plot(tspan, (max(Y_Solutions(i,:))-Y_Solutions(i,:)),'-c'); hold on;
       % plot(tspanZ, max(z_opt_four)-params(startrow,10)*zdatanorm, '*r'); hold on;
    %    ylabel('negative RLC in BM, shifted','FontWeight','bold','FontSize',12);
        %legend('HB development','Reversed BM development - fits','BM percentage','FontSize',10)
     %   legend('HB development','Reversed BM development - fits','FontSize',10)

        counter=counter+1;
        subplot(3,3,counter);


        plot(tspanYPlats',RBC'-min(RBC),'*k'); hold on;
        xlabel('Time(days)', 'FontWeight','bold','FontSize',16);
        ylabel('RBC, shifted','FontWeight','bold','FontSize',16)
        yyaxis right;
        plot(tspan, (max(Y_Solutions(i,:))-Y_Solutions(i,:)),'-c'); hold on;
        %plot(tspanZ, max(z_opt_four)-params(startrow,10)*zdatanorm, '*r'); hold on;
        ylabel('negative ALC in PB, shifted','FontWeight','bold','FontSize',12);
        %legend('RBC development','Reversed BM development - fits','BM percentage','FontSize',10)
        legend('RBC development','Reversed PB development - fits','FontSize',10)
        
        savefig(g,strcat(pwd,'/Analysis/plots/refill/Patient '+ string(patient)+'_Model_'+string(goodModels(i,model))),'compact');

        
        f=figure;
    
        plot(tspanYPlats', Neutroph','*k'); hold on;
        ylabel('Neutroph');
        yyaxis right;
        plot(tspanz, (max(Z_Solutions(i,:))-Z_Solutions(i,:)),'-c'); hold on;

        savefig(f,strcat(pwd,'/Analysis/plots/refill/Patient '+ string(patient)+'_Model_'+string(goodModels(i,model))+'Neutr.'),'compact');
    end  
    
    
end  
%f=gca;
% exportgraphics(f,[pwd,strcat('/Figures/ACC',string(1),'.png')],'Resolution',300)
%saveas(f, 'a'+string(patient), 'jpg')