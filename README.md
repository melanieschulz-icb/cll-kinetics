# cll-kinetics

Project: Modeling CLL-kinetics under Ibrutinib, Data from Jan Burger
Usage: 
- enter your Parameters into runParameters.xlsx (use template)
- save the parameter file in a folder 
- run: main_joint_run(foldername)

Folders:
- AnalyseResults: functions for model/parameter analysis after fitting the data
- Data: data
- Figures: often used as outputFolder -> empty, no code
- Helpfunctions: helpfunctions needed for fitting (data aggregation, technical functions, boundaries etc.) and model comparison
- likelihood: implementation of the likelihood functions to compute likelihood/bic, etc.
- MainArchive: archive for deprecated main functions
- MixedEffectModeling: pre- and postprocessing for mixed effect modeling, to be combined with Monolix
- ObjectiveFunctions: implementation of the objective functions for model fitting
- Solutions: implementation of model solutions
- TissueData: tissue data
