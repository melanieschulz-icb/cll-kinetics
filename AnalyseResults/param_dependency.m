% Get Params


mutated=[3,6,9,10,14,15,18,21,28:30]; %+30
unmutated=[1,2,4,5,7,11,13,16,17,19,20,22:26];
nr=[8,27];

folder=strcat(pwd,'/Figures/09_07/Plots/');

[params, txt,~] = xlsread(strcat(pwd,'/Figures/09_07/Params_best'));

paramNames=["d2","d1", "m", "alpha"];
ParamNamesX=["d2","d1", "m", "alpha"];

for paramName = paramNames
    for paramNameX = ParamNamesX
        paramMut=[];
        paramUnmut=[];
        paramNr=[];

        paramMutX=[];
        paramUnmutX=[];
        paramNrX=[];
    %extract the column of the param that is to be analysed
        [idx, idy]=find(~cellfun(@isempty, regexp(txt, paramName)));
        [idxY, idyX]=find(~cellfun(@isempty, regexp(txt, paramNameX)));
        idy=idy(1);
        idyX=idyX(1);
        for patient = mutated
            if patient < 12
                startrow=patient;
            else
                startrow=patient-1;
            end

            paramMut=[paramMut,100*params(startrow,idy)];
            paramMutX=[paramMutX,100*params(startrow,idyX)];
        end


        for patient = unmutated

            if patient < 12
                startrow=patient;
            else
                startrow=patient-1;
            end

            paramUnmut=[paramUnmut,100*params(startrow,idy)];
            paramUnmutX=[paramUnmutX,100*params(startrow,idyX)];
        end

        for patient = nr

            if patient < 12
                startrow=patient;
            else
                startrow=patient-1;
            end

            paramNr=[paramNr,100*params(startrow,idy)];
            paramNrX=[paramNrX,100*params(startrow,idyX)];
        end
        indMut=paramMut>0;
        indUnmut=paramUnmut>0;
        indNr=paramNr>0;
        
        paramMut=paramMut(indMut);
        paramUnmut=paramUnmut(indUnmut);
        paramNr=paramNr(indNr);
        
        paramMutX=paramMutX(indMut);
        paramUnmutX=paramUnmutX(indUnmut);
        paramNrX=paramNrX(indNr);
        


        sz = 60;

        figure;

        plotTitle=sgtitle('Param Analysis  - '+paramNameX+'-->'+paramName);
        

        s=scatter(paramMutX, paramMut,80,'filled','s','MarkerFaceColor',[0.6350 0.0780 0.1840]); hold on;
        scatter(paramUnmutX, paramUnmut,80,'filled','^','MarkerFaceColor',[0 0.4470 0.7410]); hold on;
        scatter(paramNrX, paramNr,80,'*','MarkerEdgeColor',[0.1290 0.6940 0.1250]); hold on;
        set(gca, 'yscale','lin');
        %ylim([0.05 10]);
        xlabel(paramNameX);
        ylabel(paramName);

        fig = gca;
        fig.YColor='black';

        f=gcf;
        

    
        [h, icons]=legend('MUTATED', 'UNMUTATED', 'NR','Fontsize',18);
        
        savefig(f,strcat(folder,plotTitle.String),'compact')
    end
end