% Get Params


mutated=[9,10,14,15,18,21,28:30]; %+30 %3,6,
unmutated=[1,2,4,5,7,11,13,16,17,20,22:26];%19

%unmutated=[1,2,4,5,7,11,13,16,19,20,22,23,25,26];
nr=[8,27];

folder=strcat(pwd,'/Analysis/Dependencies/Plots/');


[params, txt,~] = xlsread(strcat(pwd,'/Analysis/Dependencies/Params_best'));

paramNames=["d2","d1", "m1", "d3", "m2","m3"];

for paramName = paramNames
    paramMut=[];
    paramUnmut=[];
    paramNr=[];
%extract the column of the param that is to be analysed
    [idx, idy]=find(~cellfun(@isempty, regexp(txt, paramName)));

    for patient = mutated
        startrow=find(params(:,1)==patient);

        paramMut=[paramMut,100*params(startrow,idy)];
    end


    for patient = unmutated

        startrow=find(params(:,1)==patient);

        paramUnmut=[paramUnmut,100*params(startrow,idy)];
    end

    for patient = nr

        startrow=find(params(:,1)==patient);

        paramNr=[paramNr,100*params(startrow,idy)];
    end

    paramMut=paramMut(paramMut>0);
    paramUnmut=paramUnmut(paramUnmut>0);
    paramNr=paramNr(paramNr>0);


    sz = 60;

    figure;

    plotTitle=sgtitle('Param Analysis - '+paramName+'_v2');
    
    s=scatter(ones(1,length(paramMut)), paramMut,80,'filled','s','MarkerFaceColor',[0.6350 0.0780 0.1840]); hold on;
    scatter(2*ones(1,length(paramUnmut)), paramUnmut,80,'filled','^','MarkerFaceColor',[0 0.4470 0.7410]); hold on;
    scatter(3*ones(1,length(paramNr)), paramNr,80,'*','MarkerEdgeColor',[0.1290 0.6940 0.1250]); hold on;
    set(gca, 'yscale','log');
    %ylim([0.05 10]);

    fig = gca;
    fig.YColor='black';

    

    xlim([0 4]);
    [h, icons]=legend('MUTATED', 'UNMUTATED', 'NR','Fontsize',18);
    
    f=gcf;
    
    savefig(f,strcat(folder,plotTitle.String),'compact')
end