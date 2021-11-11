	/**************************************************
	* Save trade costs mata
	**************************************************/
mata	
	 st_addvar("double", ("t_BLN", "t_CFL")) 			// add variable in stata

	 mu_tilde =  st_data(.,"mu_ij_tilde")
	// BASELINE
			 W = st_data(.,"$indepI")
			 t_BLN =  (W*b) + mu_tilde
			 st_store(., "t_BLN", (t_BLN))
end
	// COUNTERFACTUAL (effects over time to be canncelled out (CFLII / BLN) )
			mata W_CFL = st_data(.,"$indepI")
			mata t_CFL =  (W_CFL*b) + mu_tilde
			mata st_store(., "t_CFL", (t_CFL))
