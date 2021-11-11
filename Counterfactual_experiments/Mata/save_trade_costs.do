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

	* COUNTERFACTUAL
			* Omit EU integration effects
			forvalues year= 2004(3)2013{
				replace EUnewx`year' = 0 if type == 4 
				replace EUnewxEU15x`year' = 0 	if type == 6 
			}
			forvalues year= 1998(3)2013{
				replace EU15x`year' = 0 if type == 5 
			}
			
			* Simulate EU and other deep agreements as average RTA
			replace schengen_noEU = 0 if schengen_noEU~=9 
			forvalues i= 3(3)6{
					replace EEA_noEU = 0 
					replace EEA_noEU_lag`i' = 0 
				replace EU_CHE_noEU = 0 
				replace EU_CHE_noEU_lag`i' = 0 
					replace EU_TUR_noEU = 0 
					replace EU_TUR_noEU_lag`i' = 0 
			}	
			
			replace RTA = 1 if EU == 1  & type~=5
			forvalues i = 3(3)6 {
				replace RTA_lag`i' = 1 if EU_lag`i' == 1  & type~=5 
			}
			*/
			
			
			mata W_CFL = st_data(.,"$indepI")
			mata t_CFL =  (W_CFL*b) + mu_tilde
			mata st_store(., "t_CFL", (t_CFL))
			
/* check direct & conditional effects
mata	ex = st_data(.,("ibn.im_id")) 
mata	type = st_data(.,("ibn.type")) 
mata	year = st_data(.,("ibn.time")) 

mata  ex'*((W-W_CFL)*b)
mata  type'*((W-W_CFL)*b)
mata year'*((W-W_CFL)*b)

/*desmat  ex_id.year im_id.year=ind($N_j), full	
			rename _x_# _FE_#, renumber 
	drop _FE_43	_FE_44 _FE_127 _FE_128 _FE_225 _FE_226 _FE_344 _FE_345 _FE_428 _FE_429 _FE_526 _FE_527
*/
estimates restore first_step
predict xp_CFL, mu


gen bln = 0
gen cfl = 0

forvalues i = 1995(3)2013{
glm x  ibn.ex_id ib($N_j).im_id if year==`i', nocon offset(t_BLN) family(poisson) vce(robust) irls
predict bln`i'  if year==`i', mu
	replace bln = bln + bln`i' if  bln`i'~=.
glm x  ibn.ex_id ib($N_j).im_id if year==`i', nocon offset(t_CFL) family(poisson) vce(robust) irls
predict cfl`i' if year==`i', mu
	replace cfl = cfl + cfl`i' if  cfl`i'~=.
}


*glm x  _FE*, nocon offset(t_BLN) family(poisson) vce(robust) irls
*predict bln, mu

*glm x _FE* , nocon offset(t_CFL) family(poisson) vce(robust) irls
*predict cfl, mu

gen double diff_conditional= ((bln-cfl)/cfl)*100
gen double diff_direct= ((xp_BLN-xp_CFL)/xp_CFL)*100


  table ex  year  , c(m diff_direct m diff_conditional)  format(%9.2f)
  table type year [aweight=cfl]   , c( m diff_conditional)  format(%9.2f)
*/