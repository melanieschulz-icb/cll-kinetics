function main_joint_run(varFolder)
% %% melanie_main_joint (Tss, PB, BM, HB)
%% for run on the server, variables are read out of the varFile saved in the specified folder
    % % can be used to fit the Burger data
    % %-> use lsq with normalization
    % close all; clearvars;
    % 


    runParameters=readtable(strcat(pwd,'/',varFolder,'/runParameters.xlsx'));
    
    %models that should be fitted
    models=runParameters.models';
    models=models(~isnan(models));
    
    duplModels=[getDuplicateModels(1003);getDuplicateModels(2003)];

    for dupl = duplModels.Variables'
        models=models(models~=dupl);
    end

     %HB model that should be fitted -> leave empty if not required
    hb_model=runParameters.hb_model';
    hb_model=hb_model(~isnan(hb_model));
    
    %model that should be used for fitting preTherapy data -> enter
    %anything if not required    
    models_pb_pre=runParameters.models_pb_pre';
    models_pb_pre=models_pb_pre(~isnan(models_pb_pre));
    
    %if pre and post therapy should be combined
    % 1: add pre therapy
    % 0: only therapy   
    combine_pre_post=runParameters.combine_pre_post';
    combine_pre_post=combine_pre_post(~isnan(combine_pre_post));

    %patients you want to fit; 12 will be excluded as ibrutinib dose was paused
    %(available patients: 1-30)
    patients=runParameters.patients';
    patients=patients(~isnan(patients));

    %add if you want to run multiple times. Format: 1:n if n runs required
    runns=runParameters.runns';
    runns=runns(~isnan(runns));

    %assumed volume of a CLL cell (in fl)
    cll_volume=runParameters.cll_volume';
    cll_volume=cll_volume(~isnan(cll_volume));


    %tissue measurements to be taken into account (max 3 available, Burger
    %takes only one)
    tss_measurements_run=runParameters.tss_measurements_run';
    tss_measurements_run=tss_measurements_run(~isnan(tss_measurements_run));

    %bm measurements to be taken into account (max 4-5 available) 
    bm_measurements=runParameters.bm_measurements';
    bm_measurements=bm_measurements(~isnan(bm_measurements));


    %if one should use the complex computation of the tissue cll burden
    use_complex_tss_run=runParameters.use_complex_tss_run';
    use_complex_tss_run=use_complex_tss_run(~isnan(use_complex_tss_run));

    %if one assumes that vanishing CLL cells in BM are instantly replaced by
    %healthy cells
    refill=runParameters.refill';
    refill=refill(~isnan(refill));

    %if one should normaize when using the least square optimizer
    normalize_lls=runParameters.normalize_lls';
    normalize_lls=normalize_lls(~isnan(normalize_lls));

    %assume multiplicative noise
    mult_noise=runParameters.mult_noise';
    mult_noise=mult_noise(~isnan(mult_noise));

    %if optimizing for sd
    opt_sd=runParameters.opt_sd';
    opt_sd=opt_sd(~isnan(opt_sd));

    %if computing/optimizing joint or seperate sd for each compartment
    seperate_sd=runParameters.seperate_sd';
    seperate_sd=seperate_sd(~isnan(seperate_sd));

    %where to save output
    folderbase_run=runParameters.folderbase_run;

        %;
    %folder where previous runs have been completed, to determine best models
%     basefolder=runParameters.basefolder';
%     basefolder=basefolder(~isnan(basefolder));

    %measurement for comparison
%     compMeasure="negLik";
    %threshold for comparison
%     threshold=0;

    %if constraint d_LT>d_PB should be used
    use_constraint=runParameters.use_constraint';
    use_constraint=use_constraint(~isnan(use_constraint));
    
    %which bounds should be used (see getBounds.m for details)
    bounds=runParameters.bounds';
    bounds=bounds(~isnan(bounds));





    %% Read the PB data
    load('heavywater-all-counts-decoded.mat');
    load('heavywater-all-counts-decoded_platelets.mat');
    load('heavy-water-patients-characterists.mat');

    %%Get VarNames
    namesModelOne=getVarNames(1003);
    ltPb=find(namesModelOne=="m_LT_PB");
    bmPb=find(namesModelOne=="m_BM_PB");
    ltBm=find(namesModelOne=="m_LT_BM");

    n_starts=runParameters.n_starts';
    n_starts=n_starts(~isnan(n_starts));



    for i=1:length(tss_measurements_run)
        %save settings
        params=[tss_measurements_run(i), use_complex_tss_run(i),normalize_lls, mult_noise, opt_sd];
        main_joint_addPre(num, txt, chars, charDesc, patients, cll_volume, runns, tss_measurements_run(i), normalize_lls, use_complex_tss_run(i), mult_noise, opt_sd,string(folderbase_run(i)),models,...
                           seperate_sd, refill, use_constraint,bounds,params, numPlat, txtPlat, hb_model, models_pb_pre,combine_pre_post,bm_measurements,...
                           namesModelOne, ltPb, bmPb, ltBm,n_starts);
    end


end
