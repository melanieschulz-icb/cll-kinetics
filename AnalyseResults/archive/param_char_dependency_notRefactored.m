%mutated=[3,6,9,10,14,15,18,21,28:30]; %+30
%unmutated=[1,2,4,5,7,11,13,16,17,19,20,22:26];
%unmutated=[1,2,4,5,7,11,13,16,19,20,22,23,25,26];

mutated=[9,10,14,15,18,21,28:30]; %+30 %3,6,
unmutated=[1,2,4,5,7,11,13,16,17,20,22:26];%19
nr=[8,27];



[params, txt,~] = xlsread(strcat(pwd,'/Figures/10_30/modelAnalysis'));
%[params, txt,~] = xlsread(strcat(pwd,'/Analysis/Dependencies/Params_best'));
[chars, txtChars,~]=xlsread(strcat(pwd,'/Heavy Water patients characterists_IPI'));

%folder=strcat(pwd,'/Figures/09_11/plots/');
%folder=strcat(pwd,'/Analysis/Dependencies/Plots/');


%paramNames=["d2","d1", "m", "alpha", "Tissue", "Blood", "redist", "died","model"];
paramNames=["paramCluster - 5","paramCluster - 6"];
charNames=["SEX_decoded","RAI","B2M","MS_decoded","F1_decoded","ZAPn","IPI"];

for paramName = paramNames
    for charName = charNames
        paramMut=[];
        paramUnmut=[];
        paramNr=[];

        charMut=[];
        charUnmut=[];
        charNr=[];
    %extract the column of the param that is to be analysed
        [idx, idy]=find(~cellfun(@isempty, regexp(txt, paramName)));
        [idxY, idyX]=find(~cellfun(@isempty, regexp(txtChars, charName)));
        idy=idy(1);
        idyX=idyX(1);
        for patient = mutated
            startrow=find(params(:,1)==patient,1);

            paramMut=[paramMut,params(startrow,idy)];
            charMut=[charMut,chars(patient,idyX)];
        end


        for patient = unmutated

            startrow=find(params(:,1)==patient,1);

            paramUnmut=[paramUnmut,params(startrow,idy)];
            charUnmut=[charUnmut,chars(patient,idyX)];
        end

        for patient = nr

            startrow=find(params(:,1)==patient,1);
            paramNr=[paramNr,params(startrow,idy)];
            charNr=[charNr,chars(patient,idyX)];
        end
        indMut=paramMut>0;
        indUnmut=paramUnmut>0;
        indNr=paramNr>0;
        
        paramMut=paramMut(indMut);
        paramUnmut=paramUnmut(indUnmut);
        paramNr=paramNr(indNr);
        
        charMut=charMut(indMut);
        charUnmut=charUnmut(indUnmut);
        charNr=charNr(indNr);
        


        sz = 60;

        figure;

        plotTitle=sgtitle('Param Analysis - '+charName+'-->'+paramName);

        s=scatter(charMut, paramMut,80,'filled','s','MarkerFaceColor',[0.6350 0.0780 0.1840]); hold on;
        scatter(charUnmut, paramUnmut,80,'filled','^','MarkerFaceColor',[0 0.4470 0.7410]); hold on;
        scatter(charNr, paramNr,80,'*','MarkerEdgeColor',[0.1290 0.6940 0.1250]); hold on;
        set(gca, 'yscale','log');
        %ylim([0.05 10]);
        xlabel(charName);
        ylabel(paramName);

        fig = gca;
        fig.YColor='black';


        f=gcf;
        

    
        [h, icons]=legend('MUTATED', 'UNMUTATED', 'NR','Fontsize',18);
        
%        savefig(f,strcat(folder,plotTitle.String),'compact')
    end
end