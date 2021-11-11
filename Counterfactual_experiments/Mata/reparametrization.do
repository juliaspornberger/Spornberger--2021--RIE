		/**************************************************	
		Use drawn parameters 
		**************************************************/
		estimates restore first_step
			matrix fe = parameters_initial[1,1..$FE]
			matrix c = parameters_initial[1,$FE+$K+1]
				matrix colname c =_cons
			matrix b_draw = fe, draw_parameters, c		
		
		erepost b = b_draw, rename
		
		/****************************
		recover pair fixed effects
		****************************/
		predict z_delta_BLN, xb								// linear prediction of regressors without pair fe
		gen double ez_delta_BLN = exp(z_delta_BLN)
			egen double ez_delta_BLN_bar = mean(ez_delta_BLN), by(pair_id)
			egen double x_bar = mean(x), by(pair_id)
			gen double emu_ij = (x_bar)/(ez_delta_BLN_bar)

		
	/**************************************************	
	* Reparameterization (Oberhofer, Pfaffermayr. 2017)
			mu_ij_tilde= mu_ij - mu_jj - mu_iN_tilde
			beta_it_tilde= beta_it - mu_iN + mu_NN 
			gamma_jt_tilde= gamma_jt + mu_jj		
	**************************************************/
		* Generate mu_jj, mu_iN, mu_NN
		gen double mu_ij = log(emu_ij)
		
		gen double temp0 = mu_ij if ex==im 				
		gen double mu_jj = (temp0) 
			qui replace mu_jj = 0 if mu_jj == .
		egen double temp_mu_NN = mean(temp0) if im_id == $N_j & ex_id == $N_j 
		gen double mu_NN = (temp_mu_NN) 
				qui replace mu_NN = 0 if mu_NN==.
			
		gen double temp1 = mu_ij if im_id==$N_j										
		egen double mu_iN=mean(temp1) if im_id==$N_j, by(ex year) 
			qui replace mu_iN = 0 if mu_iN == .
		* Genereta mu_iN_tilde (such that mu_ij_tilde = 0 for i=j and j=N)
		gen mu_iN_tilde = - mu_iN +mu_NN
			qui replace mu_iN_tilde = 0 if ex_id==im_id & ex!=im
			qui replace mu_iN_tilde = 0 if mu_iN_tilde == .
		drop temp*
		
		* Subtract from estimated fixed effects: mu_ij_tilde = mu_ij - mu_jj + mu_iN_tilde
		gen double mu_ij_tilde = mu_ij - mu_jj + mu_iN_tilde
	
	
		* Add subtracted FE to keep equivalence: z_delta_tilde = z_delta + mu_ii + mu_iN - 1 (mu_NN twice	
		gen double z_delta_tilde_BLN = z_delta_BLN + mu_jj  + mu_iN - mu_NN 
		
		
