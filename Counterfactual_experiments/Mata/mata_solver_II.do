
mata 
mata set matastrict off

		


/**********************************
* Solve
**********************************/
	// BASELINE
	p=pass()
	
	// COUNTERFACTUAL CGE
	//pCFL=pass_CFL_CGE()

	// COUNTERFACTUAL GE
	pge=pass_CFL_ge()
	

	

function results(struct pass scalar p,  struct pass_CFL_ge pge)
{
/*
real vector S 
real vector v 
real vector xp_avw 
real vector check1 
*/

	stata("use tempfile_3.dta, clear")

	stata("qui keep if year == 2013")


/**********************************************
*** Load and define Variables
**********************************************/		
	year = st_data(.,"year") 			
		displayas("result")
		printf("%f\n", max(year))	

		Dx=st_data(., "ibn.ex_id")					/*... Dummy design matrix for importer (N² x N)*/
		Dm_all=st_data(., "ibn.im_id")				/*... Dummy design matrix for exporter (N² x N-1)*/
		Dm=Dm_all[.,1..cols(Dx)-1]	
			D=(Dx,Dm)								/*... Dummy design matrix for both ex and im (N² x 2N-1)*/
			D_all=(Dx,Dm_all)								/*... Dummy design matrix for both ex and im (N² x 2N-1)*/
		type = st_data(.,("ibn.type")) 					/*... load indicator matrix for types*/
		
	kappa_i=st_data(., "kappa_i") 					/*... load production shares (repeated)*/
	theta_j=st_data(., "theta_j") 					/*... load expenditure shares (repeated)*/
		kappa_i=(Dx'*kappa_i)/cols(Dx) 						/*... Production shares (κ_i) (1 x N-1)*/
		theta_j_all=(Dm_all'*theta_j)/cols(Dx) 		/*... Expenditure shares without reference (θ_j) (1 x N)*/
		theta_j=(Dm'*theta_j)/cols(Dx) 						/*... Expenditure shares (θ_j) (1 x N)*/
	kapthe=(kappa_i\theta_j) 						/*... collect in one vector (2N-1 x 1) yi=yj .. balanced trade*/
		theta_j = theta_j\0							/*... add value for reference country to complete vector*/

	// Reload trade cost vectors for CS
	 t_BLN = st_data(.,"t_BLN") 					/*... load trade cost vector BLN*/
	 t_CFL = st_data(.,"t_CFL") 					/*... load trade cost vector CFL*/
	
	// Solve Conditional Model to get starting values
	stata("qui glm x ibn.ex_id ib($N_j).im_id , nocon offset(t_BLN) family(poisson) vce(robust)")
	start = st_matrix("e(b)")'			
		start = start[1..2*cols(Dx)-1,1]
		beta0 = start[1..cols(Dx),1]
		
	stata("qui glm x ibn.ex_id ib($N_j).im_id , nocon offset(t_CFL) family(poisson) vce(robust)")
	start_CFL = st_matrix("e(b)")'			
		start_CFL = start_CFL[1..2*cols(Dx)-1,1]
		beta0_CFL = start_CFL[1..cols(Dx),1]

displayas("text")
printf("done step: (re)load data into mata for year %f\n", max(year))	
	
	
	
	
/*********************************************
* Solver 
**********************************************/
	
 	// BASELINE	
	
	xp_avw = exp(t_BLN + D*start)

	printf("done step: solve baseline \n")	
	

	// COUNTERFACTUAL GE
	
	pge.t_CFL = t_CFL
	pge.kapthe = kapthe
	pge.D = D
	pge.Dx = Dx
	pge.kappa_i = kappa_i
	pge.theta_j = theta_j
	pge.beta0 = beta0_CFL
	pge.sigma = $sigma


	S = solvenl_init()
	solvenl_init_argument(S,1,pge)
	solvenl_init_evaluator(S, &avwge())
	solvenl_init_type(S, "zero")
	solvenl_init_technique(S, "broydenpowell")
	solvenl_init_numeq(S, cols(D))
	solvenl_init_damping(S, 0.05)
	solvenl_init_startingvals(S, start_CFL)
	solvenl_init_conv_maxiter(S, 1000)
	solvenl_init_conv_iterchng(S, 1e-3)
	solvenl_init_conv_nearzero(S, 1e-3)
	solvenl_init_iter_log(S, "off")
	C= solvenl_result_converged(S)
		vge = solvenl_solve(S)
		xp_avw_CFL_ge = exp(t_CFL + D*vge)
		//st_store(., "xp_avw_CFL_ge", (xp_avw_CFL_ge))
	
		
displayas("text")
printf("done step: solve counterfactual GE \n")	
	
		
/***********************
Results
************************/
/***** GE trade impacts *****/		
    /*Calculate (weighted) means of percentage changes by country and type*/
	ex = st_data(.,("ibn.ex_id")) 
	d_xp_CFL_ge = (xp_avw_CFL_ge :/ xp_avw) - J(rows(xp_avw), 1, 1) // observed EU to non EU
	d_xp_ex =  (ex'* d_xp_CFL_ge)/ cols(ex) 	
 	d_xp_type =  (type'*xp_avw_CFL_ge) :/ (type'*xp_avw) - J(cols(type), 1, 1)


/***** Welfare impacts (approximiation folowing Arkolakis et al. 2012) *****/
	/* Domestic trade and expenditure shares*/
	X_jj = diagonal(rowshape(xp_avw, cols(ex))')
	X_jj_GE = diagonal(rowshape(xp_avw_CFL_ge, cols(ex))')
		kappa_i_GE= st_matrix("kappa_i_GE")
		theta_j_GE = theta_j_all:* (kappa_i_GE :/ pge.kappa_i)	
			stata("collapse (sum) type if ex==im, by(ex year)")
			stata("qui replace type = 4 if type ==1")
			stata("qui replace type = 5 if type ==2")
			stata("qui replace type = 8 if type ==3")
			type_w = st_data(.,("ibn.type")) 		

	/* Calculate percentage change BLN to CFL */
	W_BLN = (X_jj:/theta_j_all)  :^(1/(1-$sigma))
	W_CFL = (X_jj_GE:/theta_j_GE)  :^(1/(1-$sigma))
	d_W_ex = (W_CFL :/ W_BLN) - J(cols(ex), 1, 1)

	/* weighted welfare impact by groups */		
	d_W_type =  (type_w'*W_CFL) :/ (type_w'*W_BLN)  - J(cols(type_w), 1, 1)


st_matrix("d_xp_ex_"+strofreal(max(year)), d_xp_ex)
st_matrix("d_xp_type_"+strofreal(max(year)), d_xp_type)
st_matrix("d_W_ex_"+strofreal(max(year)), d_W_ex)	
st_matrix("d_W_type_"+strofreal(max(year)), d_W_type)	
st_matrix("C_"+strofreal(max(year)), C)	

displayas("text")
printf("done step: store change in trade and welfare  \n")	


}

end 



