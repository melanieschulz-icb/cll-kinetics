% Get Params


mutated=[3,6,9,10,14,15,18,21,28:30]; %+30
unmutated=[1,2,4,5,7,11,13,16,17,19,20,22:26];
nr=[8,27];

folder=strcat(pwd,'/Figures/09_07/plots/');

[params, txt,~] = xlsread(strcat(pwd,'/Figures/09_07/Params_best'));

paramNames=["d2","d1", "m", "alpha", "Tissue", "Blood", "redist", "died"];

for paramName = paramNames
    paramMut=[];
    paramUnmut=[];
    paramNr=[];
%extract the column of the param that is to be analysed
    [idx, idy]=find(~cellfun(@isempty, regexp(txt, paramName)));

    for patient = mutated
        if patient < 12
            startrow=patient;
        else
            startrow=patient-1;
        end

        paramMut=[paramMut,100*params(startrow,idy)];
    end


    for patient = unmutated

        if patient < 12
            startrow=patient;
        else
            startrow=patient-1;
        end

        paramUnmut=[paramUnmut,100*params(startrow,idy)];
    end

    for patient = nr

        if patient < 12
            startrow=patient;
        else
            startrow=patient-1;
        end

        paramNr=[paramNr,100*params(startrow,idy)];
    end

    paramMut=paramMut(paramMut>0);
    paramUnmut=paramUnmut(paramUnmut>0);
    paramNr=paramNr(paramNr>0);

    cll_duration_mut=100./paramMut;
    cll_duration_unmut=100./paramUnmut;
    cll_duration_nr=100./paramNr;

    sz = 60;

    figure;
    
    plotTitle=sgtitle('Param Analysis - '+paramName);
    
    yyaxis left
    s=scatter(ones(1,length(paramMut)), paramMut,80,'filled','s','MarkerFaceColor',[0.6350 0.0780 0.1840]); hold on;
    scatter(ones(1,length(paramUnmut)), paramUnmut,80,'filled','^','MarkerFaceColor',[0 0.4470 0.7410]); hold on;
    scatter(ones(1,length(paramNr)), paramNr,80,'*','MarkerEdgeColor',[0.1290 0.6940 0.1250]); hold on;
    set(gca, 'yscale','log');


    fig = gca;
    fig.YColor='black';

    ylabel('% CLL cells dying per day');
    set(gca,'ytick',[0.1,1,10,100,1000]);
    yyaxis right
    scatter(2*ones(1,length(paramMut)), cll_duration_mut,80,'filled','s','MarkerFaceColor',[0.6350 0.0780 0.1840]); hold on;
    scatter(2*ones(1,length(paramUnmut)),cll_duration_unmut,80,'filled','^','MarkerFaceColor',[0 0.4470 0.7410]); hold on;
    scatter(2*ones(1,length(paramNr)), cll_duration_nr,80,'*','MarkerEdgeColor',[0.1290 0.6940 0.1250]); hold on;
    set(gca, 'yscale','log');

    set(gca,'xtick',1:2, 'xtickLabel',{paramName,paramName+'^-^1'},'ytick',[0.1,1,10,100,1000]);
    fig.YColor='black';

    ylabel('Average life span of tissue cells');

    xlim([0 3]);
    [h, icons]=legend('MUTATED', 'UNMUTATED', 'NR','Fontsize',18);   
    
    f=gcf;
    savefig(f,strcat(folder,plotTitle.String),'compact')
end