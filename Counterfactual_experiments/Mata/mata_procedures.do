mata 

/***************************
* Programm solver AvW
***************************/

 // BASELINE
	struct pass{
	matrix t_BLN
	vector kapthe
	matrix D


	}

	void function avw(real colvector s, real colvector values, struct pass scalar p)
	{
	values = p.D'*exp(p.t_BLN + p.D*s) - p.kapthe
	}
	//(values, p.D'*exp(p.t_BLN), p.D'*exp(p.D*s) ,  p.kapthe)

/* // COUNTERFACTUAL
	struct pass_CFL_CGE{
		matrix t_CFL
		vector kapthe
		matrix D

	}
 
	void function avwCFL(real colvector s, real colvector values, struct pass_CFL_CGE scalar pCFL)
	{
	values = pCFL.D'*exp(pCFL.t_CFL + pCFL.D*s) - pCFL.kapthe
	}
*/




/***************************
* Programm solver GE AvW
***************************/
	// COUNTERFACTUAL

	struct pass_CFL_ge {
		real matrix t_CFL
		real vector kapthe
		real matrix D
		real matrix Dx
		real matrix kappa_i
		real matrix theta_j
		real scalar sigma	
		real vector beta0

	}

	void function avwge(real colvector s, real colvector values, struct pass_CFL_ge scalar pge)
	{
	real colvector z 
	real colvector beta
	real vector d_p
	real colvector kappa_i_GE
	real vector theta_j_GE
	real vector kaptheGE
	
	z = pge.t_CFL + pge.D*s		// solve for new beta_i & gamma_j
	beta = s[1..cols(pge.Dx),1]					    			// extract beta_i (N)
	d_p = exp( (beta - pge.beta0)*(1 / (1-pge.sigma)) ) 	// calculate change of factory gate prices (delta_price)
	kappa_i_GE = pge.kappa_i :* d_p		// calculate new kappa_j
	kappa_i_GE = kappa_i_GE / sum(kappa_i_GE)
	theta_j_GE = pge.theta_j :* (kappa_i_GE :/ pge.kappa_i)
	// equivalent thmge = pge.thm :* d_p								// recalculate theta_j shares
	
	kaptheGE = (kappa_i_GE \ theta_j_GE) 	// "row-join" kappa_i & theta_j using inverted dividend sign \operator, makes a new row
	kaptheGE = kaptheGE[1..rows(kaptheGE)-1,1]				// subtract last importer for reference
	values = pge.D'*exp(z) - kaptheGE						// combine values [N^2]
	
	st_matrix("kappa_i_GE", kappa_i_GE)						// store
	st_matrix("theta_j_GE", theta_j_GE)						// store
	//z
	//s
	//beta // (40 x 1)
	// d_p //. (40 x 1)
	// kappa_i_GE .
	//theta_j_GE .
	//kaptheGE .
	//pge.D*values 
	//exp(pge.t_CFL + pge.D*values) 
		 }
		 
		 

		 
end

