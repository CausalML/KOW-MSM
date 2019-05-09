R

The data generation was done by using the code in the R files:

* Over L:
	Linear 		data_gene_3Xs_211117_L_L_S.R
	Nonlinear 	data_gene_3Xs_211117_NL_NL_S.R
* Over N:
	Linear 		data_gene_3Xs_211117_L_L_nn_S.R
	Nonlinear	data_gene_3Xs_211117_NL_NL_nn_S.R


The simulations were done by using:

* Over Sample Size
	Linear - Correct: 				- KOW:	simu_3Xs_291117_LL_linPoly1_nn_prova.R simu_3Xs_291117_LL_linPoly1_nn_prova_b.R
														- CBPS: simu_3Xs_291117_LL_linPoly1_nn_prova_CBPS.R
	Linear - Over: 						- KOW: 	simu_3Xs_291117_LL_linPoly2_nn_prova_b.R,
														- CBPS:	simu_3Xs_291117_LL_linPoly2_nn_prova.R simu_3Xs_291117_LL_linPoly2_nn_prova_CBPS
	NonLinear - Misspecified: - KOW:	simu_3Xs_291117_NLNL_linPoly1_nn_prova_b.R simu_3Xs_291117_NLNL_linPoly1_nn_prova.R
														- CBPS: simu_3Xs_291117_NLNL_linPoly1_nn_prova_CBPS.R
	NonLinear - Correct: 			- KOW: 	simu_3Xs_291117_NLNL_linPoly2_nn_prova.R
														- CBPS: simu_3Xs_291117_NLNL_linPoly2_nn_prova_CBPS.R


* Over Lambda
	Linear - Correct: simu_3Xs_291117_LL_linPoly1_LL_prova
	Linear - Over: simu_3Xs_291117_LL_linPoly2_LL_prova
	NonLinear - Misspecified: simu_3Xs_291117_NLNL_linPoly1_NLNL_prova
	NonLinear - Correct: simu_3Xs_291117_NLNL_linPoly2_NLNL_prova
