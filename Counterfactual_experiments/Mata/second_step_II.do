	/**********************		
	* - Second step: estimate time-invariant parameters
	**********************/	
		* Get EU starting values for expansion countries
		replace type = 4 if ex!=im & old_new_id==1
		replace type = 6 if ex!=im & old_new_id==3 
		replace type = 1 if ex==im & old_new_id==1

		***********************
		*Define CFLII: Border effect halves within _ years if the trend continues (Barro & Sala-i-Martin (2004))
		***********************
		/* Implied rate of reduction per year
		global year = 2013-1995
		sca ror_15=  exp(_b[BORDERx2013]/$year)-1
		sca ror_new=   exp(_b[EUnewx2013]/$year)-1
		sca ror_new15=   exp(_b[EUnewxEU15x2013]/$year)-1
				* Border within EU decreases by _% per year
				di ror_15*100
				di ror_new*100
				di ror_new15*100
				* Halving the remaining border barrier will be achieved in _ years
				sca year_ror_15=  log(2)/ror_15
				sca year_ror_new= log(2)/ror_new
				sca year_ror_new15= log(2)/ror_new15
					di year_ror_15
					di year_ror_new
					di year_ror_new15
		*/		

		save tempfile_2.dta, replace

		collapse (last) ex im ex_id im_id  RTA schengen* EU*  _x_*   ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type (mean) emu_ij, by(pair_id)	
			
		glm emu_ij _x* $indepII  , family(poisson) irls nocons vce(robust)
		estimates store second_step
		* Predictions 
			predict emup_BLN, mu			
			matrix second = e(b)	
				*CFLII: Border effect declines with exp(- alpha*t). Half the EU's border effect with different rate of reductions. ln(2) or (ror_15*year_ror_15)
				replace BORDER = ln(2) if type == 4	| type ==5 | type==6				
				*replace EUxBORDER_ij = ln(2) if type == 4	| type ==5 | type==6	
				
									
											
				predict emup_CFL, mu
					gen double ed_mu = (emup_CFL / emup_BLN)
					gen d_mu = log(ed_mu)
		save temp0.dta, replace
		
	use tempfile_2.dta, replace	
	merge m:1 pair_id using temp0.dta
	erase temp0.dta
	drop _merge
	
		mata d_mu = st_data(., "d_mu")
