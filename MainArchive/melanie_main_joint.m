%% melanie_main_joint (Tss, PB, BM, HB)
% can be used to fit the Burger data
%-> use lsq with normalization
close all; clearvars;

%% Choose parameter

%models = [3:51,52,55,58,59,61,63,65]; %:66;
%models = [1001:1006];
%models=182:309;
%models=[3:821];
%models=2003:2101;
models=[1003];
duplModels=[getDuplicateModels(1003);getDuplicateModels(2003)];

for dupl = duplModels.Variables'
    models=models(models~=dupl);
end

hb_model=[];

%patients you want to fit; 12 will be excluded as ibrutinib dose was paused
%(available patients: 1-30)
patients=[1];
runns=[1:2];

%assumed volume of a CLL cell (in fl)
cll_volume=166;

%tissue measurements to be taken into account (max 3 available, Burger
%takes only one)
tss_measurements=3;

%bm measurements to be taken into account 
bm_measurements=0;

%if one should use the complex computation of the tissue cll burden
use_complex_tss=0;

%if one assumes that vanishing CLL cells in BM are instantly replaced by
%healthy cells
refill=1;

%if one should normaize when using the least square optimizer
normalize_lls=0;

%assume multiplicative noise
mult_noise=0;

%if optimizing for sd
opt_sd=1;

%if computing/optimizing joint or seperate sd for each compartment
seperate_sd=0;

%where to save output
folderbase='/Figures/02_17_3tss_smpl/';

%folder where previous runs have been completed, to determine best models
basefolder="/Figures/02_03/";
%measurement for comparison
compMeasure="negLik";
%threshold for comparison
threshold=0;

%save settings
params=[tss_measurements, use_complex_tss,normalize_lls, mult_noise, opt_sd];

n_starts=25;

%% Read the PB data
load('heavywater-all-counts-decoded.mat');
load('heavywater-all-counts-decoded_platelets.mat');
load('heavy-water-patients-characterists.mat');

%%Get VarNames
namesModelOne=getVarNames(1003);
ltPb=find(namesModelOne=="m_LT_PB");
bmPb=find(namesModelOne=="m_BM_PB");
ltBm=find(namesModelOne=="m_LT_BM");
%% Prepare Variables for storing the opt. results


%% Perform fitting and plot the results

Zpredictions=[];
for patient=patients
    Xmulti=[];
    mkdir([folderbase(2:end),'patient_',int2str(patient),'/']);
    folder=[folderbase,'patient_',int2str(patient),'/'];
    mkdir([folder(2:end),'/plots/']);
    plotFolder=[folder,'/plots/'];
    writetable(array2table(params, 'VariableNames',{'tss_measurements', 'use_complex_tissue','normalize_lls', 'mult_noise', 'opt_sd'}),[pwd,strcat(folder,'Settings.txt')]);
    save(strcat(pwd,folderbase,'/params'),'cll_volume','tss_measurements','bm_measurements','use_complex_tss',...
        'refill');
    
    %% Set data
    
    %use available mess points as time span
    %use available mess values as true data

    [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt, patient);
    [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements);
    [tspanZAvailable,ZdataAvailable, tspanZ, zdata, ~]=getBMData(patient, startDate, bm_measurements);
                
    [~,tspanYPlats,Platelets,PlateletsScaled, HB,HBScaled, Neutroph,NeutrophScaled, RBC, RBCScaled] = getBloodCounts(numPlat, txtPlat, chars, charDesc,patient, false);
    
%    models=getBestModels(patient, basefolder, compMeasure, threshold);
    for model=models
        for runn=runns
            if patient==12
                'Skipping patient_nr 12 as therapy was paused in between';
            else
                'Fitting Patient_nr '+ string(patient);
               clearvars -except num txt chars charDesc patient use_constraint patients cll_volume Xmulti Xllh Errormulti Errorllh runn runns tss_measurements x normalize_lls tss_measurements use_complex_tss weight mult_noise opt_sd folder folderbase plotFolder model models...
                   seperate_sd refill params numPlat txtPlat hb_model bm_measurements...
                   namesModelOne ltPb bmPb ltBm basefolder compMeasure threshold Zpredictions...
                   startDate tspanY ydata ty_outlier y_outlier...
                   tspanXAvailable XdataAvailable tspanX xdata...
                   tspanZAvailable ZdataAvailable tspanZ zdata...
                   tspanYPlats Platelets PlateletsScaled  HB HBScaled  Neutroph NeutrophScaled  RBC  RBCScaled...
                   n_starts
              
               close all; 
                %% optimization using integral




                
                %% Set Bounds

                [lb_p, ub_p,lb_sd, ub_sd, par_init_p, sdXYZ] = getBounds(xdata(1), ydata(1), zdata(1),model,patient);  
                [lb_hb, ub_hb,~, ~, par_init_hb, ~] = getBoundsHB(HB(1), hb_model);  
                



                %% test log pars

                lb=log10([lb_p,lb_hb]);
                ub=log10([ub_p,ub_hb]);
                par_init=log10([par_init_p,par_init_hb]);


%                 lb=[log10([lb_p,lb_hb(1)]),-log10(-lb_hb(2)),log10(lb_hb(3:4))];
%                 ub=[log10([ub_p,ub_hb(1)]),-log10(-ub_hb(2)),log10(ub_hb(3:4))];
%                 par_init=[log10([par_init_p,par_init_hb(1)]),-log10(-par_init_hb(2)),log10(par_init_hb(3:4))];

                %% Get latin hypercube sampled starts

                [par_starts,~] = lhsdesign_modified(n_starts,lb,ub);

                startpoints = CustomStartPointSet(par_starts(:,:));

                %% define objective function.
                %if normalize_lls=true, use the appropriate objective function
                if normalize_lls==true
                    if mult_noise==true
                            lsqfun=@(par)diff_log_norm_summed(par,xdata,ydata, tspanX, tspanY);
                    else
                        lsqfun = @(par)diff_integral_norm_summed(10.^par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ);
                    end
                else
                    if mult_noise==true
                        if opt_sd==false
                            lsqfun=@(par)diff_log_summed(par,xdata,ydata,tspanX, tspanY);
                        else
                            lsqfun=@(par)diff_log_summed_sd(10.^par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill, seperate_sd);
                        end
                    else
                        if opt_sd==false
                            lsqfun = @(par)diff_integral_summed(par,xdata,ydata,tspanX, tspanY);
                        else
%                             lsqfun = @(par)diff_integral_summed_sd(10.^par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,seperate_sd,HB',tspanYPlats', hb_model);
                              lsqfun = @(par)diff_integral_summed_sd(10.^par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,seperate_sd,[],[],[]);
                        end
                    end 
                end

                %% create OptimProblem

                oldOptions = optimoptions('fmincon');
                options = optimoptions(oldOptions,'MaxFunctionEvaluations',15000,'MaxIterations',6000,'OptimalityTolerance',1e-8);

                names=getVarNames(model);
                names=names(getRelevantParams(model));
                
               aineq=zeros(1,length(names));
               aineq(1,names=="BM_0")=-1;
               aineq(1,names=="BM_healthy")=1;
     %          aineq(2,length(names)+1)=1;
     %          aineq(2,length(names)+2)=-1;
                
               problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
               'lb',lb,'ub',ub,'Aineq',aineq,'bineq',0);

%                  problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
%                  'lb',lb,'ub',ub, 'options', options);

                %% run the optimization problem 
                ms = MultiStart('PlotFcns',{@gsplotbestf, @errorlandscape},'FunctionTolerance',1e-8);
%                ms = MultiStart('FunctionTolerance',1e-8);
                [xmulti,errormulti,exitflag,output,solutions] = run(ms,problem,startpoints);
                xmulti=10.^xmulti;
                if mult_noise==1
                    [sd_x, sd_y, sd_z]=compute_mult_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd);
                    likelihood=likelihood_mult_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,seperate_sd);

                else
  %                  [sd_x, sd_y, sd_z, sd_hb]=compute_add_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd, HB',tspanYPlats', hb_model);
                    [sd_x, sd_y, sd_z, sd_hb]=compute_add_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd, [],[], hb_model);

 %                   likelihood=likelihood_add_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill ,seperate_sd,HB',tspanYPlats', hb_model);
                   likelihood=likelihood_add_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill ,seperate_sd,[],[], hb_model);
  
                    likelihood_twoComp=likelihood_add_noise(xmulti,xdata,ydata,[],tspanX, tspanY,[], model,refill ,seperate_sd,[], [], hb_model);
 %                   likelihood=likelihood_twoComp;
 
   
                end
 %               bic=(length(xmulti)+1+seperate_sd*3)*log(length([xdata, ydata, zdata,HB']))-2*log(likelihood);
                bic=(length(xmulti)+1+seperate_sd*3)*log(length([xdata, ydata, zdata,[]]))-2*log(likelihood);

                bic_twoComp=(length(xmulti)+1+seperate_sd*3)*log(length([xdata, ydata]))-2*log(likelihood_twoComp);
%                bic=bic_twoComp;
                aicr=2*(length(xmulti)+1+seperate_sd*3)-2*log(likelihood);
                aicr_twoComp=2*(length(xmulti)+1+seperate_sd*3)-2*log(likelihood_twoComp);
%                aicr=aicr_twoComp;



                %% save params and values

                storeParams=zeros(1,19+length(getRelevantHBParams(hb_model)));
                storeParams(getRelevantParams(model))=xmulti(1:end-length(getRelevantHBParams(hb_model)));
                names=getVarNames(model);
                if isempty(find(names=="m_LT_PB", 1))
                    storeParams(13)=storeParams(ltPb);
                    storeParams(ltPb)=0;
                end
                if isempty(find(names=="m_BM_PB", 1))
                    storeParams(14)=storeParams(bmPb);
                    storeParams(bmPb)=0;
                end
                if isempty(find(names=="m_LT_BM", 1))
                    storeParams(15)=storeParams(ltBm);
                    storeParams(ltBm)=0;
                end
                    
                    
                storeParams(16:19)=[sd_x, sd_y, sd_z, sd_hb];
                storeParams(20:end)=xmulti(end-length(getRelevantHBParams(hb_model))+1:end);
%                 Xmulti=[Xmulti;[model,patient,storeParams,bic, aicr,likelihood]];
                Xmulti=[Xmulti;[model,hb_model,patient,storeParams,errormulti,bic, aicr,likelihood, bic_twoComp, aicr_twoComp, likelihood_twoComp]];


                %% plot data

                tspan = 0:1:600;
                tspanx=0:1:500;
                tspanz=0:1:700;
                [x_opt, y_opt, z_opt] = solution(xmulti(1:length(getRelevantParams(model))), tspanx, tspan,tspanz, model);

                figure;

%                 subplot(2,1,1)

                %tissue measurements taken into account
                plot(tspanXAvailable(1:min(tss_measurements,length(tspanXAvailable))), XdataAvailable(1:min(tss_measurements,length(tspanXAvailable))), '*k'); hold on;
                %tissue measurements not taken into account
                if length(tspanXAvailable)>tss_measurements
                    plot(tspanXAvailable(tss_measurements+1:end), XdataAvailable(tss_measurements+1:end), 'og'); hold on;
                end

                plot(tspanx, x_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); 
                legend('measurements used for fitting','additional measurements','solution lsq');

                plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'_tissue_run_'+string(runn),'Interpreter','none');
                plotTitle.FontSize=16;

                %save figure as PDF
                f=gcf;
                exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)

                figure;
%                 subplot(2,1,1)
                plot(tspanY, ydata, '*k'); hold on;
                %outliner not taken into account


                %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
                plot(tspan, y_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); hold on;
                plot(ty_outlier, y_outlier, 'og'); 
                legend('measurements','solution lsq')

                plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'_PB_run_'+string(runn),'Interpreter','none');
                %plotTitle.Color='red';
                plotTitle.FontSize=16;
% 
% 
                %save figure as PDF
                f=gcf;
                exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)

%                 figure;
%                 yyaxis left;
% %                 plot(tspanz, z_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); hold on;
%                 plot(tspanz, z_opt, '--m'); hold on;
%                 
%                 yyaxis right;
%                 if refill
%                     bm_cells=(z_opt(1)/ZdataAvailable(1)).*ones(1,length(tspanz));
%                 else
%                     normal_cells=((1-ZdataAvailable(1))/ZdataAvailable(1))*z_opt(1);
%                     bm_cells=z_opt+normal_cells;
%                 end
% %                 plot(tspanz, 100*z_opt./bm_cells,'--r',xlabel('days'),ylabel('CLL cell percentage'));
%                 plot(tspanz, 100*z_opt./bm_cells,'-r');
%                 plot(tspanZAvailable(1:min(bm_measurements,length(tspanZAvailable))), 100*ZdataAvailable(1:min(bm_measurements,length(tspanZAvailable))),'*k'); hold on;
%                                 %bm measurements not taken into account
%                 if length(tspanZAvailable)>bm_measurements
%                     plot(tspanZAvailable(bm_measurements+1:end), 100*ZdataAvailable(bm_measurements+1:end), 'og'); hold on;
%                 end
% 
%                 legend('estimated ALC','fitted CLL percentage','measurements');
% 
%                 
% 
% 
%                 plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'_BM_run_'+string(runn),'Interpreter','none');
%                 %plotTitle.Color='red';
%                 plotTitle.FontSize=16;
% 
%                 %save figure as PDF
%                 f=gcf;
%                 exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)
                
%                  tspanPlats = 0:1:600;
%                 
%                 [hb_opt] = solutionHB(xmulti,tspanPlats,model,hb_model);
% 
%                  figure;
% %                 subplot(2,1,1)
%                 plot(tspanYPlats, HB, '*k'); hold on;
%                 %outliner not taken into account
% 
% 
%                 %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
%                 plot(tspanPlats, hb_opt, '--m',xlabel('days'),ylabel('HB value')); hold on;
% 
%                 legend('measurements','solution lsq')
% 
%                 plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'_HB_run_'+string(runn),'Interpreter','none');
%                 %plotTitle.Color='red';
%                 plotTitle.FontSize=16;
% 
% 
%                 %save figure as PDF
%                 f=gcf;
%                 exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)
% 
% % predictions:
%                 [~, ~, z_opt_pred] = solution(xmulti(1:length(getRelevantParams(model))), [], [],tspanZAvailable(1:4), model);
%                 bm_cells=(z_opt_pred(1)/ZdataAvailable(1)).*ones(1,4);
%                 Zpredictions=[Zpredictions;patient,z_opt_pred./bm_cells];
                %% Profile Likelihood Analysis

                if mult_noise==1 && opt_sd==1
                    original_lsqfun=@(par)diff_log_summed_sd(par,xdata,ydata,zdata,tspanX, tspanY, tspanZ);
                    [parval, parres] = Profiles_lsq(original_lsqfun, xmulti, 10*n_starts, [10.^lb,1e-8,1e-8, 1e-8], [10.^ub,0.2,1,1]);
                elseif mult_noise==0 && opt_sd==1
                    original_lsqfun=@(par)diff_integral_summed_sd_given(par,[sd_x,sd_y,sd_z,sd_hb],xdata,ydata,zdata,tspanX, tspanY,tspanZ,model,refill,[],[],[]);
                    [parval, parres] = Profiles_lsq(original_lsqfun, [xmulti,sd_x,sd_y,sd_z,sd_hb], 1000*n_starts, [10.^lb,10e4,10e4, 10e4,10e4], [10.^ub,10e14,10e14, 10e14,10e14]);
                else
                    [parval, parres] = Profiles_lsq(lsqfun, log10(xmulti), 1000*n_starts, lb, ub);
                end
                
    
                figure;
                plotTitle=sgtitle('parval_ls_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
                %plotTitle.Color='red';
                plotTitle.FontSize=16;
                modelNames=getVarNames(model);
                modelNames=modelNames(getRelevantParams(model));
                for ind=1:length(getRelevantParams(model))
                    subplot(3,5,ind)
                    plot(parval(ind,:),parres(ind,:));xlabel('parameter value');ylabel('sum of squared differences');title(modelNames(ind));
                end

    
                if opt_sd==true
                   subplot(3,5,12)
                   plot(parval(end-3,:),parres(end-3,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdx');
                   subplot(3,5,13)
                   plot(parval(end-2,:),parres(end-2,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdy');
                   subplot(3,5,14)
                   plot(parval(end-1,:),parres(end-1,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdz');
                   subplot(3,5,15)
                   plot(parval(end,:),parres(end,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdhb');
    
                end
    % 
    % 
    %             %save figure as PDF
    %             f=gcf;
    %             exportgraphics(f,[pwd,strcat(folder,plotTitle.String,'.png')],'Resolution',300)
    %             savefig(f,[pwd,strcat(folder,plotTitle.String)],'compact')
            end 
        end
   
    
 %% Store Parameters and errors
    end
    if opt_sd==false
        writetable(array2table(Xmulti, 'VariableNames',["model","ACC",varNames]),[pwd,strcat(folder,'Params.xlsx')]);
    else
        hb_model_store=hb_model;
        if isempty(hb_model_store)
            hb_model_store=5;
        end
        writetable(array2table(Xmulti, 'VariableNames',["model","ACC",namesModelOne,"m_PB_LT","m_PB_BM",...
            "m_BM_LT","sd_x","sd_y","sd_z","sd_hb",getVarNamesHB(hb_model),...
            "negLik","bic", "aicr","likelihood", "bic_TwoComp", "aicr_TwoComp", "likelihood_TwoComp"]),[pwd,strcat(folder,['Params_Pat',int2str(patient),'.xlsx'])]);
        save(['folder/Params_Pat',int2str(patient)],array2table(Xmulti, 'VariableNames',["model","ACC",namesModelOne,"m_PB_LT","m_PB_BM",...
            "m_BM_LT","sd_x","sd_y","sd_z","sd_hb",getVarNamesHB(hb_model),...
            "negLik","bic", "aicr","likelihood", "bic_TwoComp", "aicr_TwoComp", "likelihood_TwoComp"]));
%        writetable(array2table(Zpredictions,'VariableNames',["patient","z1","z2","z3","z4"]),[pwd,strcat(folder,['modeling_01','.xlsx'])])
    end
    
end







