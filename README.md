# lipoprotein_kinetics

The computational model for hepatic lipoprotein kinetics.
The MATLAB codes and data files published in this repository are intended to accompany the paper submitted to the Journal of Clinical Investigations:Insights.
The computational model is intended to analyze stable isotopic enrichment and biochemical data collected from bariatric surgery patients.
Provided codes and data files are created for population average data.
Each script is explained in detail with given inline remarks. 
The readers and users should refer to the published manuscript and accompanying appendix.


Published data files and scripts are as follows;

master_run.m : The file can be used to run the model.

k_matrix.m : The function that creates the coefficient matrix for the hepatic, plasma and tracer injection modules.

kdigest.m : The function that creates the coefficient matrix for gastrointestinal module.

ode.m : The function that defines the system of ordinary differential equations to represent the systems dynamics.

simul.m : The function that generates the simulations and model outputs. The function simulates, preprandial stage, meal consumption stage and postprandial stage 

and outputs 3 MATLAB solution structures.

sols.m : The function that concatenates the 3 solution structures to generate a continuous solution.

inp.mat : The data file that holds injected labelled leucine (inl) and glycerol (ing) during the kinetics study. First row is pre-surgery and second row is post surgery.

par_table.mat : Table of estimated parameter values. Columns are as follows;

	name : parameter name
	
	lb: lower bound
	
	ub: upper bound
	
	adj: 1 for estimated parameters, 0 for fixed ones
	
	est1: pre-surgery estimate
	
	std1: pre-surgery parameter standard deviations
	
	r1: ratio of std to estimate (std1/est1)
	
	est2: post-surgery estimate
	
	std2: post-surgery parameter standard deviations
	
	r2: ratio of std to estimate (std2/est2)

mata.mat : Cell structure that holds population averages for biochemical concentration and enrichment data. The first cell is for pre-surgery and the second cell is for post-surgery data. Each cell contains 15 variables and for each variable the first row is measurement time points (hour), second row is population average and the third row is population standard deviation. The variables and their units are as follows;

	vldl1_tg_pool: Plasma VLDL1 TG concentration (mg/L)
	
	vldl2_tg_pool: Plasma VLDL2 TG concentration (mg/L)
	
	vldl1_apob_pool: Plasma VLDL1 apoB concentration (mg/L)
	
	vldl2_apob_pool: Plasma VLDL2 apoB concentration (mg/L)	
	
	vldl1_glyc_enr: Plasma VLDL1 TG glycerol enrichment
	
	vldl2_glyc_enr: Plasma VLDL2 TG glycerol enrichment
	
	vldl1_leu_enr: Plasma VLDL1 apoB leucine enrichment
	
	vldl2_leu_enr: Plasma VLDL2 apoB leucine enrichment
	
	pl_gly_enr: Plasma glycerol enrichment
	
	pl_leu_enr: Plasma leucine enrichment
	
	PL_TG : Total plasma TG concentration (mg/L)
	
	CM_TG : Plasma chylomicron TG concentration (mg/L)
	
	insulin : Plasma insulin concentration (pool/L)
	
	glucose : Plasma glucose concentration (mmol/L)
	
	hmusig: Estimated alpha, beta and delta values used in f function in ode.m these values are estimated from
	plasma insulin data by fiting a Gaussian curve.
