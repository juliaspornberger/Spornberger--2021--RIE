	/**************************************************
	* Second step
	**************************************************/
	
	save tempfile_2.dta, replace
		sort   pair_id year
		
			
		collapse (first) ex im ex_id im_id  RTA schengen* EU* EEA* (mean) emu_ij (max) _x_*   ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type, by(pair_id)	
	
	
	glm emu_ij _x* $indepII  , family(poisson) irls nocons vce(robust)

			* Predictions
			predict emup_BLN, mu			
				* Define CFL experiment I: EU Integration 
				* To obtain the starting value of EU trade flows, set border effect within the EU equal to average border effect
				replace EUxBORDER_ij = 0 
				predict emup_CFL, mu
					gen double ed_mu = (emup_CFL / emup_BLN)
					gen d_mu = log(ed_mu)
		save temp0.dta, replace
		
	use tempfile_2.dta, replace	
	merge m:1 pair_id using temp0.dta
	erase temp0.dta
	drop _merge
	
		mata d_mu = st_data(., "d_mu")
