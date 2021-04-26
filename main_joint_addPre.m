function main_joint_addPre(num, txt, chars, charDesc, patients, cll_volume, runns, tss_measurements, normalize_lls, use_complex_tss, mult_noise, opt_sd,folderbase,models,...
                       seperate_sd, refill, use_constraint,bounds,params, numPlat, txtPlat, hb_model,pre_therapy_models,combine_pre_post, bm_measurements,...
                       namesModelOne, ltPb, bmPb, ltBm,n_starts)
    for patient=patients
        Xmulti=[];
        mkdir(strcat(convertStringsToChars(folderbase),'patient_',int2str(patient),'/'));
        folder=strcat("/",folderbase,'patient_',int2str(patient),'/');
        mkdir(strcat(convertStringsToChars(folder.extractAfter(1)),'/plots/'));
        plotFolder=strcat(folder,'/plots/');
        writetable(array2table(params, 'VariableNames',{'tss_measurements', 'use_complex_tissue','normalize_lls', 'mult_noise', 'opt_sd'}),strcat(pwd,folder,'Settings.txt'));
        save(strcat(pwd,"/",folderbase,'/params'),'cll_volume','tss_measurements','bm_measurements','use_complex_tss',...
            'combine_pre_post','refill','n_starts');

        %% Set data

        %use available mess points as time span
        %use available mess values as true data

        [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt, patient,combine_pre_post);

         
        [tspanXAvailable, XdataAvailable, tspanX, xdata]= aggTissueData(patient, cll_volume, use_complex_tss, startDate, tss_measurements);
        [~, XdataAvailableRaw, ~, ~]= aggTissueData(patient, cll_volume, false, startDate, tss_measurements);
        [tspanZAvailable,ZdataAvailable, tspanZ, zdata, ~]=getBMData(patient, startDate, bm_measurements);

        [~,tspanYPlats,Platelets,PlateletsScaled, HB,HBScaled, Neutroph,NeutrophScaled, RBC, RBCScaled] = getBloodCounts(numPlat, txtPlat, chars, charDesc,patient, false);

    %    models=getBestModels(patient, basefolder, compMeasure, threshold);
%         load("compTable.mat");
%         patComp=compTable(compTable.ACC==patient,:);
%         models=unique(patComp(patComp.rating<3,"model").Variables');
%         models=[models,1003,1004,1005,1006,1016,2003,2004,2006,2016];
%         models=unique(models);
        
        if ~combine_pre_post
            pre_therapy_models=-1;
        end
        for model=models
            for model_pb_pre=pre_therapy_models

                for runn=runns
                    if patient==12
                        'Skipping patient_nr 12 as therapy was paused in between';
                    else
                        'Fitting Patient_nr '+ string(patient);
                       clearvars -except num txt chars charDesc patient use_constraint bounds patients cll_volume Xmulti Xllh Errormulti Errorllh runn runns tss_measurements x normalize_lls tss_measurements use_complex_tss weight mult_noise opt_sd folder folderbase plotFolder model models...
                           seperate_sd refill use_constraint params numPlat txtPlat hb_model bm_measurements...
                           pre_therapy_models model_pb_pre combine_pre_post...
                           namesModelOne ltPb bmPb ltBm basefolder compMeasure threshold Zpredictions...
                           startDate tspanY ydata ty_outlier y_outlier...
                           tspanXAvailable XdataAvailable tspanX xdata XdataAvailableRaw...
                           tspanZAvailable ZdataAvailable tspanZ zdata...
                           tspanYPlats Platelets PlateletsScaled  HB HBScaled  Neutroph NeutrophScaled  RBC  RBCScaled...
                           n_starts

                       close all; 
                        %% optimization using integral





                        %% Set Bounds

                        [lb_p, ub_p,~, ~, par_init_p, ~] = getBounds(XdataAvailable(1),XdataAvailableRaw(end),ydata(tspanY==0),ZdataAvailable(1), model,patient,combine_pre_post,bounds);  
                        [lb_hb, ub_hb,~, ~, par_init_hb, ~] = getBoundsHB(HB(1), hb_model); 
                        if model_pb_pre<=3 || model_pb_pre==7
                            pb_0=ydata(1);
                        elseif model_pb_pre<=6
                            pb_0=ydata(tspanY==0);
                        else
                            error("Specify if pb_o is first or last datapoint");
                        end

                        [lb_pre, ub_pre,~, ~, par_init_pre, ~] = getBounds_PB_pre(pb_0, model_pb_pre);


                        if patient==patients(1) && runn==1
                            save(strcat(pwd,"/",folderbase,"boundaries.mat"));
                        end

                        %% test log pars

                        lb=log10([lb_p,lb_hb,lb_pre]);
                        ub=log10([ub_p,ub_hb,ub_pre]);
                        par_init=log10([par_init_p,par_init_hb,par_init_pre]);


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
                                      lsqfun = @(par)diff_integral_summed_sd(10.^par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,seperate_sd,[],[],[],model_pb_pre);
                                end
                            end 
                        end

                        %% create OptimProblem

                        oldOptions = optimoptions('fmincon');
                        options = optimoptions(oldOptions,'MaxFunctionEvaluations',15000,'MaxIterations',6000,'OptimalityTolerance',1e-8);

                        names=getVarNames(model);
                        names=names(getRelevantParams_comb(model,combine_pre_post));
                        
                        if use_constraint
                            aineq=zeros(1,length(names)+length(getRelevantParams_PB_pre(model_pb_pre)));
                            aineq(1,names=="d_LT")=-1;
                            aineq(1,names=="d_PB")=1;
                            problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
                        'lb',lb,'ub',ub,'Aineq',aineq,'bineq',[0],'options',options);
                        else
%                        aineq=zeros(3,length(names)+length(getRelevantParams_PB_pre(model_pb_pre)));
%                        aineq(1,names=="BM_0")=-1;
%                        aineq(1,names=="BM_healthy")=1;
% %                        aineq(2,names=="PB_0")=-1;
% %                        aineq(2,names=="Blood_healthy")=1;
%                        aineq(3,names=="LT_0")=-1;
%                        aineq(3,names=="LT_healthy")=1;
             %          aineq(2,length(names)+1)=1;
             %          aineq(2,length(names)+2)=-1;

%                        problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
%                        'lb',lb,'ub',ub,'Aineq',aineq,'bineq',[0;0;0],'options',options);

                         problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
                         'lb',lb,'ub',ub, 'options', options);

                        end
%                          problemPar = problem('particleswarm','objective',lsqfun,...
%                          'lb',lb,'ub',ub, 'options', options);
                        %% run the optimization problem 
                        ms = MultiStart('PlotFcns',{@gsplotbestf, @errorlandscape},'FunctionTolerance',1e-8);
        %                ms = MultiStart('FunctionTolerance',1e-8);
                        [xmulti,errormulti,exitflag,output,solutions] = run(ms,problem,startpoints);
        %                [xpsw,fvalpsw,eflagpsw,outputpsw] = particleswarm(lsqfun,12,lb,ub);
        %                options = optimoptions('patternsearch','MaxIterations',50000);
        %                [xps,errormulti,eflagps,outputps] = patternsearch(lsqfun,par_init,[],[],[],[],lb,ub,[],options);
                        xmulti=10.^xmulti;
        %                xmulti=10.^xps;
                        if mult_noise==1
                            [sd_x, sd_y, sd_z]=compute_mult_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd);
                            likelihood=likelihood_mult_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,seperate_sd);

                        else
          %                  [sd_x, sd_y, sd_z, sd_hb]=compute_add_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd, HB',tspanYPlats', hb_model);
                            [sd_x, sd_y, sd_z, sd_hb]=compute_add_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd, [],[], hb_model,model_pb_pre);

         %                   likelihood=likelihood_add_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill ,seperate_sd,HB',tspanYPlats', hb_model);
                           likelihood=likelihood_add_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill ,seperate_sd,[],[], hb_model,model_pb_pre);

                            likelihood_twoComp=likelihood_add_noise(xmulti,xdata,ydata,[],tspanX, tspanY,[], model,refill ,seperate_sd,[], [], hb_model,model_pb_pre);
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

                        storeParams=zeros(1,19+length(getRelevantHBParams(hb_model))+2);
                        storeParams(getRelevantParams_comb(model,combine_pre_post))=xmulti(1:end-length(getRelevantHBParams(hb_model))-length(getRelevantParams_PB_pre(model_pb_pre)));
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
                        storeParams(20:end-2+length(getRelevantParams_PB_pre(model_pb_pre)))=xmulti(end-length(getRelevantHBParams(hb_model))-length(getRelevantParams_PB_pre(model_pb_pre))+1:end);
        %                 Xmulti=[Xmulti;[model,patient,storeParams,bic, aicr,likelihood]];
                        Xmulti=[Xmulti;[model,hb_model,model_pb_pre,patient,storeParams,errormulti,bic, aicr,likelihood, bic_twoComp, aicr_twoComp, likelihood_twoComp]];


                        %% plot data
%                         if combine_pre_post    
%                             tspan_pre=-100:1:0;
%                         else
%                             tspan_pre=[];
%                         end
%                         tspan = 0:1:600;
%                         tspanx=0:1:500;
%                         tspanz=0:1:700;
%                         y_opt_pre=solution_PB_pre(xmulti(end-1:end),tspan_pre,model_pb_pre);
%                         par=xmulti(1:length(getRelevantParams_comb(model,combine_pre_post)));
%                         if ~isempty(tspan_pre)
%                             par_nums=getRelevantParams(model);
%                             
%                             ind_pb0=find(par_nums==5,1);
%                             par=[par(1:ind_pb0-1),y_opt_pre(end),par(ind_pb0:end)];
%                         end
%                         [x_opt, y_opt, z_opt] = solution(par, tspanx, tspan,tspanz, model);
%                         y_opt=[y_opt_pre,y_opt];
%                         figure;
%     
%         %                 subplot(2,1,1)
%     
%                         %tissue measurements taken into account
%                         plot(tspanXAvailable(1:min(tss_measurements,length(tspanXAvailable))), XdataAvailable(1:min(tss_measurements,length(tspanXAvailable))), '*k'); hold on;
%                         %tissue measurements not taken into account
%                         if length(tspanXAvailable)>tss_measurements
%                             plot(tspanXAvailable(tss_measurements+1:end), XdataAvailable(tss_measurements+1:end), 'og'); hold on;
%                         end
%     
%                         plot(tspanx, x_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); 
%                         legend('measurements used for fitting','additional measurements','solution lsq');
%     
%                         plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'pre_model'+string(model_pb_pre)+'_tissue_run_'+string(runn),'Interpreter','none');
%                         plotTitle.FontSize=16;
%     
%                         %save figure as PDF
%                         f=gcf;
%                         exportgraphics(f,strcat(pwd,plotFolder,plotTitle.String,'.png'),'Resolution',300);
%     
%                         figure;
%         %                 subplot(2,1,1)
%                         plot(tspanY, ydata, '*k'); hold on;
%                         %outliner not taken into account
%     
%     
%                         %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
%                         plot([tspan_pre,tspan], y_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); hold on;
%                         plot(ty_outlier, y_outlier, 'og'); 
%                         legend('measurements','solution lsq')
%     
%                         plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'pre_model'+string(model_pb_pre)+'_PB_run_'+string(runn),'Interpreter','none');
%                         %plotTitle.Color='red';
%                         plotTitle.FontSize=16;
%         % 
%         % 
%                         %save figure as PDF
%                         f=gcf;
%                         exportgraphics(f,strcat(pwd,plotFolder,plotTitle.String,'.png'),'Resolution',300);
%     
%                         figure;
%                         yyaxis left;
%         %                 plot(tspanz, z_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); hold on;
%                         plot(tspanz, z_opt, '--m'); hold on;
%                         
%                         yyaxis right;
%                         if refill
%                             bm_cells=(z_opt(1)/ZdataAvailable(1)).*ones(1,length(tspanz));
%                         else
%                             normal_cells=((1-ZdataAvailable(1))/ZdataAvailable(1))*z_opt(1);
%                             bm_cells=z_opt+normal_cells;
%                         end
%         %                 plot(tspanz, 100*z_opt./bm_cells,'--r',xlabel('days'),ylabel('CLL cell percentage'));
%                         plot(tspanz, 100*z_opt./bm_cells,'-r');
%                         plot(tspanZAvailable(1:min(bm_measurements,length(tspanZAvailable))), 100*ZdataAvailable(1:min(bm_measurements,length(tspanZAvailable))),'*k'); hold on;
%                                         %bm measurements not taken into account
%                         if length(tspanZAvailable)>bm_measurements
%                             plot(tspanZAvailable(bm_measurements+1:end), 100*ZdataAvailable(bm_measurements+1:end), 'og'); hold on;
%                         end
%         
%                         legend('estimated ALC','fitted CLL percentage','measurements');
%         
%                         
%         
%         
%                         plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'_mod_pre_'+string(model_pb_pre)+'_BM_run_'+string(runn),'Interpreter','none');
%                         %plotTitle.Color='red';
%                         plotTitle.FontSize=16;
        
                        %save figure as PDF
%                         f=gcf;
%                         exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)

%                          tspanPlats = 0:1:600;
%                         
%                         [hb_opt] = solutionHB(xmulti,tspanPlats,model,hb_model);
%         
%                          figure;
%         %                 subplot(2,1,1)
%                         plot(tspanYPlats, HB, '*k'); hold on;
%                         %outliner not taken into account
%         
%         
%                         %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
%                         plot(tspanPlats, hb_opt, '--m',xlabel('days'),ylabel('HB value')); hold on;
%         
%                         legend('measurements','solution lsq')
%         
%                         plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'_HB_run_'+string(runn),'Interpreter','none');
%                         %plotTitle.Color='red';
%                         plotTitle.FontSize=16;
%         
%         
%                         %save figure as PDF
%                         f=gcf;
%                         exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)
        % 
        % % predictions:
        %                 [~, ~, z_opt_pred] = solution(xmulti(1:length(getRelevantParams(model))), [], [],tspanZAvailable(1:4), model);
        %                 bm_cells=(z_opt_pred(1)/ZdataAvailable(1)).*ones(1,4);
        %                 Zpredictions=[Zpredictions;patient,z_opt_pred./bm_cells];
%                         %% Profile Likelihood Analysis
%     
%                         if mult_noise==1 && opt_sd==1
%                             original_lsqfun=@(par)diff_log_summed_sd(par,xdata,ydata,zdata,tspanX, tspanY, tspanZ);
%                             [parval, parres] = Profiles_lsq(original_lsqfun, xmulti, 10*n_starts, [10.^lb,1e-8,1e-8, 1e-8], [10.^ub,0.2,1,1]);
%                         elseif mult_noise==0 && opt_sd==1
%                             original_lsqfun=@(par)diff_integral_summed_sd_given(par,[sd_x,sd_y,sd_z,sd_hb],xdata,ydata,zdata,tspanX, tspanY,tspanZ,model,refill,[],[],[]);
%                             [parval, parres] = Profiles_lsq(original_lsqfun, [xmulti,sd_x,sd_y,sd_z,sd_hb], 1000*n_starts, [10.^lb,10e4,10e4, 10e4,10e4], [10.^ub,10e14,10e14, 10e14,10e14]);
%                         else
%                             [parval, parres] = Profiles_lsq(lsqfun, log10(xmulti), 1000*n_starts, lb, ub);
%                         end
%     
%     
%                         figure;
%                         plotTitle=sgtitle('parval_ls_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
%                         %plotTitle.Color='red';
%                         plotTitle.FontSize=16;
%                         modelNames=getVarNames(model);
%                         modelNames=modelNames(getRelevantParams(model));
%                         for ind=1:length(getRelevantParams(model))
%                             subplot(3,5,ind)
%                             plot(parval(ind,:),parres(ind,:));xlabel('parameter value');ylabel('sum of squared differences');title(modelNames(ind));
%                         end
%     
%     
%                         if opt_sd==true
%                            subplot(3,5,12)
%                            plot(parval(end-3,:),parres(end-3,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdx');
%                            subplot(3,5,13)
%                            plot(parval(end-2,:),parres(end-2,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdy');
%                            subplot(3,5,14)
%                            plot(parval(end-1,:),parres(end-1,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdz');
%                            subplot(3,5,15)
%                            plot(parval(end,:),parres(end,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdhb');
%     
%                         end
%             % 
%             % 
%                         %save figure as PDF
%                         f=gcf;
%                         exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)
%                         savefig(f,strcat(pwd,plotFolder,plotTitle.String),'compact');
                    end 
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
            namesModelOne
            tab=array2table(Xmulti, 'VariableNames',convertStringsToChars(["model","model_pre","ACC",namesModelOne,"m_PB_LT","m_PB_BM",...
                "m_BM_LT","sd_x","sd_y","sd_z","sd_hb",getVarNamesHB(hb_model),"PB_pre_0","g",...
                "negLik","bic", "aicr","likelihood", "bic_TwoComp", "aicr_TwoComp", "likelihood_TwoComp"]));
            writetable(tab,strcat(pwd,folder,['Params_Pat',int2str(patient),'.xlsx']));
            save(strcat(pwd,folder,'/Params_Pat',int2str(patient)),'tab');
    %        writetable(array2table(Zpredictions,'VariableNames',["patient","z1","z2","z3","z4"]),[pwd,strcat(folder,['modeling_01','.xlsx'])])
        end

    end
end