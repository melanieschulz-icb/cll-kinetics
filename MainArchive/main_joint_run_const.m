% %% melanie_main_joint (Tss, PB, BM, HB)
% % can be used to fit the Burger data
% %-> use lsq with normalization
% close all; clearvars;
% 
% %% Choose parameter
% 
%models = [3:51,52,55,58,59,61,63,65]; %:66;
%models = [1001:1006];
%models=182:309;
%models=[3:821];
%models=2003:2101;
models=[0.5];
duplModels=[getDuplicateModels(1003);getDuplicateModels(2003)];

for dupl = duplModels.Variables'
    models=models(models~=dupl);
end

hb_model=[];
models_pb_pre=[1];
combine_pre_post=0;

%patients you want to fit; 12 will be excluded as ibrutinib dose was paused
%(available patients: 1-30)
patients=[1:11,13:30];
runns=[1];

%assumed volume of a CLL cell (in fl)
cll_volume=166;

%tissue measurements to be taken into account (max 3 available, Burger
%takes only one)
tss_measurements_run=[0,1,2,2,3,3];
%tss_measurements_run=[2];

%bm measurements to be taken into account 
bm_measurements=0;

%if one should use the complex computation of the tissue cll burden
use_complex_tss_run=[0,1,0,1,0,1];

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
folderbase_run=["Figures/04_07_0tss_constraint/";"Figures/04_07_1tss_constraint/";"Figures/04_07_2tss_s_constraint/";"Figures/04_07_2tss_c_constraint/";"Figures/04_07_3tss_s_constraint/";"Figures/04_07_3tss_c_constraint/"];
%;
%folder where previous runs have been completed, to determine best models
basefolder="/Figures/02_03/";
%measurement for comparison
compMeasure="negLik";
%threshold for comparison
threshold=0;

%if constraint d_LT>d_PB should be used
use_constraint=1;




%% Read the PB data
load('heavywater-all-counts-decoded.mat');
load('heavywater-all-counts-decoded_platelets.mat');
load('heavy-water-patients-characterists.mat');

%%Get VarNames
namesModelOne=getVarNames(1003);
ltPb=find(namesModelOne=="m_LT_PB");
bmPb=find(namesModelOne=="m_BM_PB");
ltBm=find(namesModelOne=="m_LT_BM");

n_starts=40;

for i=1:length(folderbase_run)
    %save settings
    params=[tss_measurements_run(i), use_complex_tss_run(i),normalize_lls, mult_noise, opt_sd];
    main_joint_addPre(num, txt, chars, charDesc, patients, cll_volume, runns, tss_measurements_run(i), normalize_lls, use_complex_tss_run(i), mult_noise, opt_sd,folderbase_run(i),models,...
                       seperate_sd, refill, use_constraint,params, numPlat, txtPlat, hb_model, models_pb_pre,combine_pre_post,bm_measurements,...
                       namesModelOne, ltPb, bmPb, ltBm,n_starts);
end



% close all; clearvars;

%% Choose parameter

% models = [3,4,367,102,316,376,1003:1101]; %:66;
% models = [1001:1006];
% models=182:309;
% models=[3:821];
% models=[0,1];
% % models=[3,4,376,102,316,376,1003,1004,1006,1016,1041,1047,2003,2004,2006,2016];
% duplModels=[getDuplicateModels(1003);getDuplicateModels(2003)];
% 
% for dupl = duplModels.Variables'
%     models=models(models~=dupl);
% end
% 
% hb_model=[];
% 
% % patients you want to fit; 12 will be excluded as ibrutinib dose was paused
% % (available patients: 1-30)
% patients=[1:11,13:30];
% runns=[1:2];
% 
% %assumed volume of a CLL cell (in fl)
% cll_volume=166;
% 
% %tissue measurements to be taken into account (max 3 available, Burger
% %takes only one)
% tss_measurements_run=[1,2,2,3,3];
% 
% 
% %bm measurements to be taken into account 
% bm_measurements=0;
% 
% %if one should use the complex computation of the tissue cll burden
% use_complex_tss_run=[0,0,1,0,1];
% 
% %if one assumes that vanishing CLL cells in BM are instantly replaced by
% %healthy cells
% refill=1;
% 
% %if one should normaize when using the least square optimizer
% normalize_lls=0;
% 
% %assume multiplicative noise
% mult_noise=0;
% 
% %if optimizing for sd
% opt_sd=1;
% 
% %if computing/optimizing joint or seperate sd for each compartment
% seperate_sd=0;
% 
% %where to save output
% folderbase_run=["Figures/02_21_1ss_2comp/","Figures/02_21_2ss_smpl_2comp/","Figures/02_21_2ss_cmpl_2comp/",...
%     "Figures/02_21_3ss_smpl_2comp/","Figures/02_21_3ss_cmpl_2comp/"];
% 
% %folder where previous runs have been completed, to determine best models
% basefolder="/Figures/02_20_4_2c/";
% %measurement for comparison
% compMeasure="negLik";
% threshold=0;
% 
% 
% 
% 
% % Read the PB data
% load('heavywater-all-counts-decoded.mat');
% load('heavywater-all-counts-decoded_platelets.mat');
% load('heavy-water-patients-characterists.mat');
% 
% %Get VarNames
% namesModelOne=getVarNames(1003);
% ltPb=find(namesModelOne=="m_LT_PB");
% bmPb=find(namesModelOne=="m_BM_PB");
% ltBm=find(namesModelOne=="m_LT_BM");
% 
% n_starts=25;
% 
% for i=1:length(folderbase_run)
%     %save settings
%     params=[tss_measurements_run(i), use_complex_tss_run(i),normalize_lls, mult_noise, opt_sd];
%     main_joint(num, txt, chars, charDesc, patients, cll_volume, runns, tss_measurements_run(i), normalize_lls, use_complex_tss_run(i), mult_noise, opt_sd,folderbase_run(i),models,...
%                        seperate_sd, refill, params, numPlat, txtPlat, hb_model, bm_measurements,...
%                        namesModelOne, ltPb, bmPb, ltBm,n_starts);
% end
% 
% % %% 2 compartment models
% % 
% % %% melanie_main_joint (Tss, PB, BM, HB)
% % % can be used to fit the Burger data
% % %-> use lsq with normalization
% % close all; clearvars;
% % 
% %% Choose parameter
% 
% %models = [3:51,52,55,58,59,61,63,65]; %:66;
% %models = [1001:1006];
% %models=182:309;
% %models=[3:821];
% %models=2003:2101;
% % models=[0,1];
% % duplModels=[getDuplicateModels(1003);getDuplicateModels(2003)];
% % 
% % for dupl = duplModels.Variables'
% %     models=models(models~=dupl);
% % end
% % 
% % hb_model=[];
% % 
% % %patients you want to fit; 12 will be excluded as ibrutinib dose was paused
% % %(available patients: 1-30)
% % patients=[1:11,13:30];
% % runns=[1:2];
% % 
% % %assumed volume of a CLL cell (in fl)
% % cll_volume=166;
% % 
% % %tissue measurements to be taken into account (max 3 available, Burger
% % %takes only one)
% % tss_measurements_run=[3];
% % 
% % %bm measurements to be taken into account 
% % bm_measurements=0;
% % 
% % %if one should use the complex computation of the tissue cll burden
% % use_complex_tss_run=[1];
% % 
% % %if one assumes that vanishing CLL cells in BM are instantly replaced by
% % %healthy cells
% % refill=1;
% % 
% % %if one should normaize when using the least square optimizer
% % normalize_lls=0;
% % 
% % %assume multiplicative noise
% % mult_noise=0;
% % 
% % %if optimizing for sd
% % opt_sd=1;
% % 
% % %if computing/optimizing joint or seperate sd for each compartment
% % seperate_sd=0;
% % 
% % %where to save output
% % folderbase_run=["Figures/02_20_5_2c"];
% % 
% % %folder where previous runs have been completed, to determine best models
% % basefolder="/Figures/02_03/";
% % %measurement for comparison
% % compMeasure="negLik";
% % %threshold for comparison
% % threshold=0;
% % 
% % 
% % 
% % 
% % %% Read the PB data
% % load('heavywater-all-counts-decoded.mat');
% % load('heavywater-all-counts-decoded_platelets.mat');
% % load('heavy-water-patients-characterists.mat');
% % 
% % %%Get VarNames
% % namesModelOne=getVarNames(1003);
% % ltPb=find(namesModelOne=="m_LT_PB");
% % bmPb=find(namesModelOne=="m_BM_PB");
% % ltBm=find(namesModelOne=="m_LT_BM");
% % 
% % n_starts=25;
% % 
% % for i=1:5
% %     %save settings
% %     params=[tss_measurements_run(i), use_complex_tss_run(i),normalize_lls, mult_noise, opt_sd];
% %     main_joint(num, txt, chars, charDesc, patients, cll_volume, runns, tss_measurements_run(i), normalize_lls, use_complex_tss_run(i), mult_noise, opt_sd,folderbase_run(i),models,...
% %                        seperate_sd, refill, params, numPlat, txtPlat, hb_model, bm_measurements,...
% %                        namesModelOne, ltPb, bmPb, ltBm,n_starts);
% % end