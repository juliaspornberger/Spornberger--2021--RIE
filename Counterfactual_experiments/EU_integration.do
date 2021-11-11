args seed Breps
*local Breps  3
*local seed 1234

clear all
*matrix drop _all
*timer clear
clear mata
capture set matsize 10000 
capture log close 
capture set more off
program drop _all
sca drop _all



********************************************************************************
*********************EU-Border Effect*******************************************
********************************************************************************
qui{
* Chosen arguments
global Breps = `Breps'
global seed = `seed'
set seed $seed

global sigma = 6
sca number=0

*Load data
use Data/WIOD_long_manu.dta, clear


/*******************
* 3-year interval
*******************/
	keep if year==1995 | year==1998 | year==2001 | year==2004  | year==2007 | year==2010 | year==2013

* Define id's
egen ex_id=group(ex)
egen im_id=group(im)
egen ex_time_id=group(ex year)
egen im_time_id=group(im year)
egen pair_id=group(ex im) 
egen time=group(year)
egen id=group(ex im year)
xtset pair_id time

}

/*******************
* Define variables for baseline and counterfactual scenario
*******************/
* Define variables for baseline and counterfactual scenario
	* First step	
xi i.schengen_noEU
		rename _Ischengen* schengen_noEU*
		
	global indepI = "RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013"
	
	global K: word count $indepI
	di $K
			

	* Second step (CONTROLS explaining endogeneity of RTA (EU): lDISTii lDISTij noCONTIG noCOLONY noCOMLANG)
	global indepII = "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij"
	global K_II: word count $indepII
	di $K_II

qui{		
*******************
* Dummies for conditional estimator
*******************
	*Choose reference country
		qui tab im
		global N_j=r(r)

			* Origin-time & destination-time FE [(N-1)(T-1) + N(T-1) + N^2 ]
			desmat  ex_id.year=ind($N_j) im_id.year=ind($N_j)
			gen _x_499=0
				replace _x_499=1 if ex_id==$N_j & year==1998		
			gen _x_500=0
				replace _x_500=1 if ex_id==$N_j & year==2001
			gen _x_501=0
				replace _x_501=1 if ex_id==$N_j & year==2004
			gen _x_502=0
				replace _x_502=1 if ex_id==$N_j & year==2007
			gen _x_503=0
				replace _x_503=1 if ex_id==$N_j & year==2010
			gen _x_504=0
				replace _x_504=1 if ex_id==$N_j & year==2013
			* Drop those reference dummy of CHE(7), HRV(19), NOR(33) in 2001 (do not have one yet)
			drop _x_37 _x_108 _x_191 _x_286 _x_357 _x_440	
			rename _x_# _I_#, renumber 
		unab vars : _I*
		global FE  `: word count `vars''
		di $FE
		
		* FE for solver [(N-1)T + NT + (N-1)^2 ]
		desmat  ex_id.year im_id.year=ind($N_j), full	
					rename _x_# _FE_#, renumber 
			drop _FE_43	_FE_44 _FE_127 _FE_128 _FE_225 _FE_226 _FE_344 _FE_345 _FE_428 _FE_429 _FE_526 _FE_527
				
		* (N-1) pair FE for 2) step
		desmat  ex_id im_id=ind($N_j) , full  
				global exporter_fe = "$term1"
				global importer_fe = "$term2"
			rename _x_# _x_#, renumber 
		unab vars : _x*
		global FE_II  `: word count `vars''
		di $FE_II
				
				

			
		

********************************************************************************
* Step 1) estimate time-variant and time-invariant parameters using the two step procedure
********************************************************************************

* Exclude EU from RTA
replace RTA = RTA_noEU
forvalues i= 1(1)12 {
	replace RTA_lag`i' =RTA_lag`i'_noEU
}	



	/**********************		
	 - First step: estimate time-variant parameters
	**********************/		
	qui ppmlhdfe x _I_* $indepI , absorb(i.im_id#i.ex_id) cluster(pair) nocons d(ln_pairFE)
		estimates store first_step
		matrix parameters_initial = e(b)
		predict xp_BLN, mu
		
		** Bias Correction **
		qui ppmlhdfe x  $indepI, absorb(im_id#year ex_id#year i.im_id#i.ex_id) cluster(im_id#ex_id) nocons d
		predict lambda
        matrix beta = e(b)
		
       ppml_fe_bias x  $indepI, i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) 
		
			* Store parameters b (without FE and constant)	
			matrix b_initial = e(b)	
			matrix V_initial = e(V)
		
	

*********************
*********************
*********************
* Start Program
*********************
*********************
*********************
save tempfile_1.dta, replace


program define bootprogram, eclass

	clear mata 
	preserve	
		/**************************************************
		* Random draw from estimated parameter distribution 
		**************************************************/
		drawnorm $indepI , means(b_initial) cov(V_initial) n(1)  clear
	
		mkmat _all, matrix(draw_parameters)
		
		use tempfile_1.dta,clear 	
		/**************************************************
		* Reparametrization
		**************************************************/
		qui do Counterfactual_experiments/Mata/reparametrization.do	
		/* switch parameters with randomly drawn parameters*/
			mata b = st_matrix("e(b)")'						
			mata b = b[$FE+1..$FE+$K,.]	

		display as result "done step: reparametrization"	


		/**************************************************
		* Save trade costs in mata
		**************************************************/
		qui do Counterfactual_experiments/Mata/save_trade_costs.do
		display as result "done step: save trade costs"	
		
		save tempfile_3.dta, replace
		
		/**************
		* Load mata procedures
		**************/				
		qui do Counterfactual_experiments/Mata/mata_procedures.do

		/***********************
		* Solve Baseline and Counterfactual model and calculate trade changes (nolevel)
		************************/
		display as result "start step: solve "	
			 do Counterfactual_experiments/Mata/mata_solver.do
			mata: results(p, pge)
	
		forvalues year = 1995(3)2013{
			matrix GE_type_`year' = (d_xp_type_`year')'
			matrix GE_country_`year' = (d_xp_ex_`year')'
			matrix W_type_`year' = (d_W_type_`year')'
			matrix W_country_`year' = (d_W_ex_`year')'
		}
		
		mat nl_GE_country = (GE_country_1995 , GE_country_1998 , GE_country_2001 , GE_country_2004 , GE_country_2007 , GE_country_2010 , GE_country_2013)
		mat nl_GE_type = (GE_type_1995 , GE_type_1998 , GE_type_2001 , GE_type_2004 , GE_type_2007 , GE_type_2010 , GE_type_2013)
		mat nl_W_country = (W_country_1995 , W_country_1998 , W_country_2001 , W_country_2004 , W_country_2007 , W_country_2010 , W_country_2013)
		mat nl_W_type = (W_type_1995 , W_type_1998 , W_type_2001 , W_type_2004 , W_type_2007 , W_type_2010 , W_type_2013)
		mat nl_C = (C_1995 , C_1998 , C_2001 , C_2004 , C_2007 , C_2010 , C_2013)

		mata mata drop results()

		/**************************************************			
		* Second step: estimate time variant parameters
		**************************************************/
		use tempfile_3.dta, replace
		qui do Counterfactual_experiments/Mata/second_step.do
			mata t_CFL = t_CFL + d_mu
			mata st_store(., "t_CFL", (t_CFL))
		display as result "done step: second step"	
			unab vars : _x*
			global FE_II  `: word count `vars''
			*di $FE_II
			matrix second = e(b)
			matrix second = second[.,$FE_II+1..$FE_II+$K_II]	
			matrix draw_parameters = (draw_parameters, second)

		save tempfile_3.dta, replace
		
		/**************
		* Load mata procedures
		**************/		
		*clear mata
		*qui do Counterfactual_experiments/Mata/mata_procedures.do
	
		/***********************
		* Solve Baseline and Counterfactual model and calculate trade changes (with level)
		************************/
		display as result "start step: solve "	
			 do Counterfactual_experiments/Mata/mata_solver.do
			mata: results(p, pge)

		/***********************
		* Store results by type and country
		***********************/
		* Combine GE results for each year
		forvalues year = 1995(3)2013{
			matrix GE_type_`year' = (d_xp_type_`year')'
			matrix GE_country_`year' = (d_xp_ex_`year')'
			matrix W_type_`year' = (d_W_type_`year')'
			matrix W_country_`year' = (d_W_ex_`year')'
		}
		

		mat GE_country = (GE_country_1995 , GE_country_1998 , GE_country_2001 , GE_country_2004 , GE_country_2007 , GE_country_2010 , GE_country_2013)
		mat GE_type = (GE_type_1995 , GE_type_1998 , GE_type_2001 , GE_type_2004 , GE_type_2007 , GE_type_2010 , GE_type_2013)
		mat W_country = (W_country_1995 , W_country_1998 , W_country_2001 , W_country_2004 , W_country_2007 , W_country_2010 , W_country_2013)
		mat W_type = (W_type_1995 , W_type_1998 , W_type_2001 , W_type_2004 , W_type_2007 , W_type_2010 , W_type_2013)
		mat C = (C_1995 , C_1998 , C_2001 , C_2004 , C_2007 , C_2010 , C_2013)
	
		
		* Rename columns GE country
		local colnames
		local nl_colnames
		forvalues year = 1995(3)2013 {
			forvalues ex = 1/43 {
				local colnames `colnames' " xp_c`=`ex'`year ''" 
				local nl_colnames `nl_colnames' " nl_xp_c`=`ex'`year ''" 
			}
		}
		matrix colnames GE_country = `colnames'
		matrix colnames nl_GE_country = `nl_colnames'
			* Rename columns W country		
			local colnames
			local nl_colnames
			forvalues year = 1995(3)2013 {
				forvalues ex = 1/43 {
					local colnames `colnames' " w_c`=`ex'`year ''" 
					local nl_colnames `nl_colnames' " nl_w_c`=`ex'`year ''" 
				}
			}
			matrix colnames W_country = `colnames'
			matrix colnames nl_W_country = `nl_colnames'
		* Rename columns GE type
		local colnames
		local nl_colnames
		forvalues year = 1995(3)2013 {
			forvalues type = 1/8 {
				local colnames `colnames' " xp_t`=`type'`year ''" 
				local nl_colnames `nl_colnames' " nl_xp_t`=`type'`year ''" 
			}
		}
		matrix colnames GE_type = `colnames'
		matrix colnames nl_GE_type = `nl_colnames'
			* Rename columns W type		
			local colnames
			local nl_colnames
			forvalues year = 1995(3)2013 {
				forvalues type = 1/3 {
					local colnames `colnames' " w_t`=`type'`year ''" 
					local nl_colnames `nl_colnames' " nl_w_t`=`type'`year ''" 
				}
			}
			matrix colnames W_type = `colnames'
			matrix colnames nl_W_type = `nl_colnames'
			
			* colnames C
			local C_colnames
			local nl_C_colnames
			forvalues year = 1995(3)2013 {
					local C_colnames `C_colnames' " c`=`year ''" 
					local nl_C_colnames `nl_C_colnames' " nl_c`=`year ''" 
			}
			matrix colnames C = `C_colnames'
			matrix colnames nl_C = `nl_C_colnames'


		* Store coefficient estimates
		matrix results = (draw_parameters, GE_country, GE_type, W_country, W_type, nl_GE_country, nl_GE_type, nl_W_country, nl_W_type, nl_C, C)
		
		* Store in ereturn b (drops all other values in e)
		ereturn post results

		di as red "done run" " number= " `seed'  "
*********************
*********************
*********************
* end Program
*********************
*********************

	restore
end


********************************************************************************
}


simulate _b,  reps($Breps): bootprogram

save Estimates/EU_integration/EU_integration.dta, replace


xxxx

save Estimates/EU_integration/EU_integration.dta, replace
 *save Estimates/EU_integration/EU_integration_paper.dta, replace
 
erase tempfile_1.dta
erase tempfile_2.dta
erase tempfile_3.dta

********************************************************************************
********************************************************************************
* Organize estimates
********************************************************************************
qui{
****************************************************
* CFL experiment I: estimates by type (nolevel)
****************************************************
*use Estimates/EU_integration/EU_integration.dta, clear
use Estimates/EU_integration/EU_integration_paper.dta, clear

missings dropobs
foreach v of varlist _b_c* {
	drop if `v'==.
}
keep _b_nl*t* 

*** Rename variables ***
forvalues year=1995(3)2015{
	forvalues type = 1(1)8{
		rename _b_nl_xp_t`type'`year' xp_`year'_`type'
	}
	forvalues type = 1(1)3{
		rename _b_nl_w_t`type'`year' w_`year'_`type'
	}
}

foreach v of varlist  w* {
	qui sum `v',d
}

*** Reshape ***
gen sim_id= _n
	reshape long xp_1995_  xp_1998_ xp_2001_ xp_2004_ xp_2007_ xp_2010_ xp_2013_ w_1995_  w_1998_ w_2001_ w_2004_ w_2007_ w_2010_ w_2013_, i(sim_id) j(type)
	rename xp_*_ xp_*
	rename w_*_ w_*
gen id = _n
	reshape long xp_ w_ , i(id) j(year)
	drop id sim_id
	rename xp_ xp
	rename w_ w
*graph box w if type==2 , over(year)
*graph box xp if type==5 , over(year)

replace xp =xp+1
replace w = w+1
gen double lnxp=(xp)


*** Store (geometric) mean estimates by country-year and year means into Matrix ***
*matrix year_mean = J(1,7,.)
sca ct = 0
forvalues year = 1995(3)2013{	
	* Generate matrices
	matrix xp_mean_`year' = J(1,8,.)
	matrix xp_CI_l_`year' = J(1,8,.)
	matrix xp_CI_u_`year' = J(1,8,.)
	matrix w_mean_`year' = J(1,8,.)
	matrix w_CI_l_`year' = J(1,8,.)
	matrix w_CI_u_`year' = J(1,8,.)
	* Create yearly mean
			*qui ameans b_ if year == `year' & b_ !=0
			*sca ct= ct+1
			*local time = ct
			*di `time'
			*matrix year_mean[1, `time'] = r(mean_g) 
	* Groupwise trade
	forvalues type=1(1)8{
		*sum lnxp if type==`type' & year == `year'
		*matrix xp_mean_`year'[1, `type'] = exp(r(mean))
		qui ameans xp if type==`type' & year == `year'
		matrix xp_mean_`year'[1, `type'] = r(mean_g)
			qui centile xp if type==`type' & year == `year', centile(2.5 97.5)
			matrix xp_CI_l_`year'[1, `type'] = r(c_1)
			matrix xp_CI_u_`year'[1, `type'] = r(c_2)
	}
		* Groupwise welfare
	forvalues type=1(1)3{
		qui ameans w if type==`type' & year == `year'
		matrix w_mean_`year'[1, `type'] = r(mean)
			qui centile w if type==`type' & year == `year', centile(2.5 97.5)
			matrix w_CI_l_`year'[1, `type'] = r(c_1)
			matrix w_CI_u_`year'[1, `type'] = r(c_2)
	}

}



*** Combine estimates for all years ***
mat xp=(xp_mean_1995\ xp_mean_1998 \ xp_mean_2001 \ xp_mean_2004 \ xp_mean_2007 \xp_mean_2010 \ xp_mean_2013  )'
mat xp_CI_l=(xp_CI_l_1995\ xp_CI_l_1998 \ xp_CI_l_2001 \ xp_CI_l_2004 \ xp_CI_l_2007 \xp_CI_l_2010 \ xp_CI_l_2013  )'
mat xp_CI_u=(xp_CI_u_1995\ xp_CI_u_1998 \ xp_CI_u_2001 \ xp_CI_u_2004 \ xp_CI_u_2007 \xp_CI_u_2010 \ xp_CI_u_2013  )'
mat w=(w_mean_1995\ w_mean_1998 \ w_mean_2001 \ w_mean_2004 \ w_mean_2007 \w_mean_2010 \ w_mean_2013  )'
mat w_CI_l=(w_CI_l_1995\ w_CI_l_1998 \ w_CI_l_2001 \ w_CI_l_2004 \ w_CI_l_2007 \w_CI_l_2010 \ w_CI_l_2013  )'
mat w_CI_u=(w_CI_u_1995\ w_CI_u_1998 \ w_CI_u_2001 \ w_CI_u_2004 \ w_CI_u_2007 \w_CI_u_2010 \ w_CI_u_2013  )'


mat colnames xp_CI_l = 1995 1998 2001 2004 2007 2010 2013 
mat colnames xp_CI_u = 1995 1998 2001 2004 2007 2010 2013 
mat colnames xp = 1995 1998 2001 2004 2007 2010 2013 
mat rownames xp_CI_l = 1 2 3 4 5 6 7 8
mat rownames xp_CI_u = 1 2 3 4 5 6 7 8
mat rownames xp = 1 2 3 4 5 6 7 8
mat colnames w_CI_l = 1995 1998 2001 2004 2007 2010 2013 
mat colnames w_CI_u = 1995 1998 2001 2004 2007 2010 2013 
mat colnames w = 1995 1998 2001 2004 2007 2010 2013 
mat rownames w_CI_l = 1 2 3 
mat rownames w_CI_u = 1 2 3 
mat rownames w = 1 2 3 

clear
svmat xp, names(matcol)
svmat xp_CI_l, names(matcol)
svmat xp_CI_u, names(matcol)
svmat w, names(matcol)
svmat w_CI_l, names(matcol)
svmat w_CI_u, names(matcol)

* Percentage changes
foreach v of varlist  _all {
	replace `v' = (`v'-1)*100
}

gen type = _n
reshape long xp xp_CI_l xp_CI_u w w_CI_l w_CI_u , i(type) j(year)

gen type_w = type
	replace type_w = 4 if type_w==1
	replace type_w = 5 if type_w==2
	replace type_w = 8 if type_w==3

	*Define labels	
		label define type_label 1 "EUnew (domestic)" 2 "EU15 (domestic)" ///
								3 "ROW (domestic)" 4 "EUnew" ///
								5 "EU15" 6 "EUnew-EU15, EU15-EUnew" ///
								7 "EU-ROW, ROW-EU" 8 "ROW", replace 
		label values type type_label
		label values type_w type_label

		
*save Estimates/EU_integration/EU_integration_nolevel_type.dta, replace
 save Estimates/EU_integration/EU_integration_nolevel_type_paper.dta, replace

di as red "done step: save results by groups of trade flows (no level)"

****************************************************
* CFL experiment I: estimates by type (level)
****************************************************
*use Estimates/EU_integration/EU_integration.dta, clear
use Estimates/EU_integration/EU_integration_paper.dta, clear

missings dropobs
foreach v of varlist _b_c* {
	drop if `v'==.
}
keep _b_*t* 
drop *nl*


*** Rename variables ***
forvalues year=1995(3)2015{
	forvalues type = 1(1)8{
		rename _b_xp_t`type'`year' xp_`year'_`type'
	}
	forvalues type = 1(1)3{
		rename _b_w_t`type'`year' w_`year'_`type'
	}
}


*** Reshape ***
gen sim_id= _n
	reshape long xp_1995_  xp_1998_ xp_2001_ xp_2004_ xp_2007_ xp_2010_ xp_2013_ w_1995_  w_1998_ w_2001_ w_2004_ w_2007_ w_2010_ w_2013_, i(sim_id) j(type)
	rename xp_*_ xp_*
	rename w_*_ w_*
gen id = _n
	reshape long xp_ w_ , i(id) j(year)
	drop id sim_id
	rename xp_ xp
	rename w_ w
	
*** Store (geometric) mean estimates by country-year and year means into Matrix ***
replace xp =xp+1
replace w = w+1	

*matrix year_mean = J(1,7,.)
sca ct = 0
forvalues year = 1995(3)2013{	
	* Generate matrices
	matrix xp_mean_`year' = J(1,8,.)
	matrix xp_CI_l_`year' = J(1,8,.)
	matrix xp_CI_u_`year' = J(1,8,.)
	matrix w_mean_`year' = J(1,8,.)
	matrix w_CI_l_`year' = J(1,8,.)
	matrix w_CI_u_`year' = J(1,8,.)
	* Create yearly mean
			*qui ameans b_ if year == `year' & b_ !=0
			*sca ct= ct+1
			*local time = ct
			*di `time'
			*matrix year_mean[1, `time'] = r(mean_g) 
	* Groupwise trade
	forvalues type=1(1)8{
		qui ameans xp if type==`type' & year == `year'
		matrix xp_mean_`year'[1, `type'] = r(mean_g)
			qui centile xp if type==`type' & year == `year', centile(2.5 97.5)
			matrix xp_CI_l_`year'[1, `type'] = r(c_1)
			matrix xp_CI_u_`year'[1, `type'] = r(c_2)
	}
		* Groupwise welfare
	forvalues type=1(1)3{
		qui ameans w if type==`type' & year == `year'
		matrix w_mean_`year'[1, `type'] = r(mean_g)
			qui centile w if type==`type' & year == `year', centile(2.5 97.5)
			matrix w_CI_l_`year'[1, `type'] = r(c_1)
			matrix w_CI_u_`year'[1, `type'] = r(c_2)
	}

}



*** Combine estimates for all years ***
mat xp=(xp_mean_1995\ xp_mean_1998 \ xp_mean_2001 \ xp_mean_2004 \ xp_mean_2007 \xp_mean_2010 \ xp_mean_2013  )'
mat xp_CI_l=(xp_CI_l_1995\ xp_CI_l_1998 \ xp_CI_l_2001 \ xp_CI_l_2004 \ xp_CI_l_2007 \xp_CI_l_2010 \ xp_CI_l_2013  )'
mat xp_CI_u=(xp_CI_u_1995\ xp_CI_u_1998 \ xp_CI_u_2001 \ xp_CI_u_2004 \ xp_CI_u_2007 \xp_CI_u_2010 \ xp_CI_u_2013  )'
mat w=(w_mean_1995\ w_mean_1998 \ w_mean_2001 \ w_mean_2004 \ w_mean_2007 \w_mean_2010 \ w_mean_2013  )'
mat w_CI_l=(w_CI_l_1995\ w_CI_l_1998 \ w_CI_l_2001 \ w_CI_l_2004 \ w_CI_l_2007 \w_CI_l_2010 \ w_CI_l_2013  )'
mat w_CI_u=(w_CI_u_1995\ w_CI_u_1998 \ w_CI_u_2001 \ w_CI_u_2004 \ w_CI_u_2007 \w_CI_u_2010 \ w_CI_u_2013  )'


mat colnames xp_CI_l = 1995 1998 2001 2004 2007 2010 2013 
mat colnames xp_CI_u = 1995 1998 2001 2004 2007 2010 2013 
mat colnames xp = 1995 1998 2001 2004 2007 2010 2013 
mat rownames xp_CI_l = 1 2 3 4 5 6 7 8
mat rownames xp_CI_u = 1 2 3 4 5 6 7 8
mat rownames xp = 1 2 3 4 5 6 7 8
mat colnames w_CI_l = 1995 1998 2001 2004 2007 2010 2013 
mat colnames w_CI_u = 1995 1998 2001 2004 2007 2010 2013 
mat colnames w = 1995 1998 2001 2004 2007 2010 2013 
mat rownames w_CI_l = 1 2 3 
mat rownames w_CI_u = 1 2 3 
mat rownames w = 1 2 3 

clear
svmat xp, names(matcol)
svmat xp_CI_l, names(matcol)
svmat xp_CI_u, names(matcol)
svmat w, names(matcol)
svmat w_CI_l, names(matcol)
svmat w_CI_u, names(matcol)
	* Percentage changes
	foreach v of varlist  _all {
		replace `v' = (`v'-1)*100
	}
	
gen type = _n
reshape long xp xp_CI_l xp_CI_u w w_CI_l w_CI_u , i(type) j(year)

gen type_w = type
	replace type_w = 4 if type_w==1
	replace type_w = 5 if type_w==2
	replace type_w = 8 if type_w==3

	*Define labels	
		label define type_label 1 "EUnew (domestic)" 2 "EU15 (domestic)" ///
								3 "ROW (domestic)" 4 "EUnew" ///
								5 "EU15" 6 "EUnew-EU15, EU15-EUnew" ///
								7 "EU-ROW, ROW-EU" 8 "ROW", replace 
		label values type type_label
		label values type_w type_label

		
*save Estimates/EU_integration/EU_integration_level_type.dta, replace
 save Estimates/EU_integration/EU_integration_level_type_paper.dta, replace

di as red "done step: save results by groups of trade flows (level)"

*/
****************************************************
* CFL experiment I: estimates by country (nolevel)
****************************************************
*use Estimates/EU_integration/EU_integration.dta, clear
use Estimates/EU_integration/EU_integration_paper.dta, clear

missings dropobs, force
gen id= 0
foreach v of varlist _b_c*  {
	drop if `v'==.
}

keep *nl*
*** Rename variables ***
forvalues year=1995(3)2015{
	forvalues country = 1(1)43{
	rename _b_nl_xp_c`country'`year' xp`year'_`country'
	rename _b_nl_w_c`country'`year' w`year'_`country'
	}
}

keep xp* w*
*** Reshape ***
gen sim_id= _n
	reshape long xp1995_ xp1998_ xp2001_ xp2004_ xp2007_ xp2010_ xp2013_ w1995_ w1998_ w2001_ w2004_ w2007_ w2010_ w2013_ , i(sim_id) j(country)
	rename xp*_ xp*
	rename w*_ w*
gen id = _n
	reshape long xp w , i(id) j(year)
	drop id sim_id
*** Store mean estimates by country-year and (geometric) year means into Matrix ***
*matrix trade_year_mean = J(1,7,.)
replace xp =xp+1
replace w = w+1	


sca ct = 0
forvalues year = 1995(3)2013{	
	matrix xp_mean_`year' = J(1,43,.)
	matrix w_mean_`year' = J(1,43,.)
	matrix xp_CI_l_`year' = J(1,43,.)
	matrix xp_CI_u_`year' = J(1,43,.)
			*qui ameans trade if year == `year' & trade !=0
			*sca ct= ct+1
			*local time = ct
			*di `time'
			*matrix trade_year_mean[1, `time'] = r(mean_g)	
	forvalues country=1(1)43{
		qui ameans xp if country==`country' & year == `year'
		matrix xp_mean_`year'[1, `country'] = r(mean_g)
			qui centile xp if country==`country' & year == `year', centile(2.5 97.5)
			matrix xp_CI_l_`year'[1, `country'] = r(c_1)
			matrix xp_CI_u_`year'[1, `country'] = r(c_2)		
				qui ameans w if country==`country' & year == `year'
				matrix w_mean_`year'[1, `country'] = r(mean_g)
	}
}



*** Combine estimates for all years ***
mat xp=(xp_mean_1995\ xp_mean_1998 \ xp_mean_2001 \ xp_mean_2004 \ xp_mean_2007 \xp_mean_2010 \ xp_mean_2013  )'
	mat xp_CI_l=(xp_CI_l_1995\ xp_CI_l_1998 \ xp_CI_l_2001 \ xp_CI_l_2004 \ xp_CI_l_2007 \xp_CI_l_2010 \ xp_CI_l_2013  )'
	mat xp_CI_u=(xp_CI_u_1995\ xp_CI_u_1998 \ xp_CI_u_2001 \ xp_CI_u_2004 \ xp_CI_u_2007 \xp_CI_u_2010 \ xp_CI_u_2013  )'

mat w=(w_mean_1995\ w_mean_1998 \ w_mean_2001 \ w_mean_2004 \ w_mean_2007 \w_mean_2010 \ w_mean_2013  )'


mat colnames xp = 1995 1998 2001 2004 2007 2010 2013 
	mat colnames xp_CI_l = 1995 1998 2001 2004 2007 2010 2013 
	mat colnames xp_CI_u = 1995 1998 2001 2004 2007 2010 2013 
mat colnames w = 1995 1998 2001 2004 2007 2010 2013

	mat rownames xp = Australia Austria Belgium Bulgaria Brazil Canada Switzerland China Cyprus CzechRepublic Germany Denmark Spain Estonia Finland France GreatBritain Greece Croatia Hungary India Indonesia Ireland Italy Japan SouthKorea Lithuania Luxembourg Latvia Mexico Malta Netherlands Norway Poland Portugal Romania Russia Slovakia Slovenia Sweden Turkey Taiwan USA  
		mat rownames xp_CI_l = Australia Austria Belgium Bulgaria Brazil Canada Switzerland China Cyprus CzechRepublic Germany Denmark Spain Estonia Finland France GreatBritain Greece Croatia Hungary India Indonesia Ireland Italy Japan SouthKorea Lithuania Luxembourg Latvia Mexico Malta Netherlands Norway Poland Portugal Romania Russia Slovakia Slovenia Sweden Turkey Taiwan USA  
		mat rownames xp_CI_u = Australia Austria Belgium Bulgaria Brazil Canada Switzerland China Cyprus CzechRepublic Germany Denmark Spain Estonia Finland France GreatBritain Greece Croatia Hungary India Indonesia Ireland Italy Japan SouthKorea Lithuania Luxembourg Latvia Mexico Malta Netherlands Norway Poland Portugal Romania Russia Slovakia Slovenia Sweden Turkey Taiwan USA  

	mat rownames w = Australia Austria Belgium Bulgaria Brazil Canada Switzerland China Cyprus CzechRepublic Germany Denmark Spain Estonia Finland France GreatBritain Greece Croatia Hungary India Indonesia Ireland Italy Japan SouthKorea Lithuania Luxembourg Latvia Mexico Malta Netherlands Norway Poland Portugal Romania Russia Slovakia Slovenia Sweden Turkey Taiwan USA  

clear
svmat xp, names(matcol)
svmat xp_CI_l, names(matcol)
svmat xp_CI_u, names(matcol)
svmat w, names(matcol)

	* Percentage changes
	foreach v of varlist  _all {
		replace `v' = (`v'-1)*100
	}

gen country_id = _n

label define country_label 1 "Australia" 2 "Austria" 3 "Belgium" 4 "Bulgaria" 5 "Brazil" 6 "Canada" 7 "Switzerland" 8 "China" 9 "Cyprus" 10 "CzechRepublic" 11 "Germany" 12 "Denmark" 13 "Spain" 14 "Estonia" 15 "Finland" 16 "France" 17 "GreatBritain" 18 "Greece" 19 "Croatia" 20 "Hungary" 21 "India" 22 "Indonesia" 23 "Ireland" 24 "Italy" 25 "Japan" 26 "SouthKorea" 27 "Lithuania" 28 "Luxembourg" 29 "Latvia" 30 "Mexico" 31 "Malta" 32 "Netherlands" 33 "Norway" 34 "Poland" 35 "Portugal" 36 "Romania" 37 "Russia" 38 "Slovakia" 39 "Slovenia" 40 "Sweden" 41 "Turkey" 42 "Taiwan" 43 "USA" 44 "Mean" , modify
label values country country_label


*save Estimates/EU_integration/EU_integration_nolevel_country.dta, replace
save Estimates/EU_integration/EU_integration_nolevel_country_paper.dta, replace
di as red "done step: save results countrywise (no level)"



****************************************************
* CFL experiment I: estimates by country (level)
****************************************************
*use Estimates/EU_integration/EU_integration.dta, clear
use Estimates/EU_integration/EU_integration_paper.dta, clear

missings dropobs
foreach v of varlist _b_c* {
	drop if `v'==.
}

drop *nl*

*** Rename variables ***
forvalues year=1995(3)2015{
	forvalues country = 1(1)43{
	rename _b_xp_c`country'`year' xp`year'_`country'
	rename _b_w_c`country'`year' w`year'_`country'
	}
}

keep xp* w*
*** Reshape ***
gen sim_id= _n
	reshape long xp1995_ xp1998_ xp2001_ xp2004_ xp2007_ xp2010_ xp2013_ w1995_ w1998_ w2001_ w2004_ w2007_ w2010_ w2013_ , i(sim_id) j(country)
	rename xp*_ xp*
	rename w*_ w*
gen id = _n
	reshape long xp w , i(id) j(year)
	drop id sim_id
*** Store mean estimates by country-year and (geometric) year means into Matrix ***
replace xp =xp+1
replace w = w+1	

sca ct = 0
forvalues year = 1995(3)2013{	
	matrix xp_mean_`year' = J(1,43,.)
	matrix w_mean_`year' = J(1,43,.)
	forvalues country=1(1)43{
		qui ameans xp if country==`country' & year == `year'
		matrix xp_mean_`year'[1, `country'] = r(mean_g)
			qui ameans w if country==`country' & year == `year'
			matrix w_mean_`year'[1, `country'] = r(mean_g)
	}
}



*** Combine estimates for all years ***
mat xp=(xp_mean_1995\ xp_mean_1998 \ xp_mean_2001 \ xp_mean_2004 \ xp_mean_2007 \xp_mean_2010 \ xp_mean_2013  )'
mat w=(w_mean_1995\ w_mean_1998 \ w_mean_2001 \ w_mean_2004 \ w_mean_2007 \w_mean_2010 \ w_mean_2013  )'


mat colnames xp = 1995 1998 2001 2004 2007 2010 2013 
mat colnames w = 1995 1998 2001 2004 2007 2010 2013

	mat rownames xp = Australia Austria Belgium Bulgaria Brazil Canada Switzerland China Cyprus CzechRepublic Germany Denmark Spain Estonia Finland France GreatBritain Greece Croatia Hungary India Indonesia Ireland Italy Japan SouthKorea Lithuania Luxembourg Latvia Mexico Malta Netherlands Norway Poland Portugal Romania Russia Slovakia Slovenia Sweden Turkey Taiwan USA  
	mat rownames w = Australia Austria Belgium Bulgaria Brazil Canada Switzerland China Cyprus CzechRepublic Germany Denmark Spain Estonia Finland France GreatBritain Greece Croatia Hungary India Indonesia Ireland Italy Japan SouthKorea Lithuania Luxembourg Latvia Mexico Malta Netherlands Norway Poland Portugal Romania Russia Slovakia Slovenia Sweden Turkey Taiwan USA  

clear
svmat xp, names(matcol)
svmat w, names(matcol)
	* Percentage changes
	foreach v of varlist  _all {
		replace `v' = (`v'-1)*100
	}

gen country_id = _n

label define country_label 1 "Australia" 2 "Austria" 3 "Belgium" 4 "Bulgaria" 5 "Brazil" 6 "Canada" 7 "Switzerland" 8 "China" 9 "Cyprus" 10 "CzechRepublic" 11 "Germany" 12 "Denmark" 13 "Spain" 14 "Estonia" 15 "Finland" 16 "France" 17 "GreatBritain" 18 "Greece" 19 "Croatia" 20 "Hungary" 21 "India" 22 "Indonesia" 23 "Ireland" 24 "Italy" 25 "Japan" 26 "SouthKorea" 27 "Lithuania" 28 "Luxembourg" 29 "Latvia" 30 "Mexico" 31 "Malta" 32 "Netherlands" 33 "Norway" 34 "Poland" 35 "Portugal" 36 "Romania" 37 "Russia" 38 "Slovakia" 39 "Slovenia" 40 "Sweden" 41 "Turkey" 42 "Taiwan" 43 "USA" 44 "Mean" , modify
label values country country_label


*save Estimates/EU_integration/EU_integration_level_country.dta, replace
save Estimates/EU_integration/EU_integration_level_country_paper.dta, replace

di as red "done step: save results countrywise (level)"


****************************************************
* CFL experiment I: coefficient estimates
****************************************************
*use Estimates/EU_integration/EU_integration.dta, clear
use Estimates/EU_integration/EU_integration_paper.dta, clear

	**************************************************
	* Generate matrices with coefs and CI for first & second step
	**************************************************
	missings dropobs, force
	foreach v of varlist _b_c* {
		drop if `v'==.
	}
	rename _b_* *
	drop xp* w* nl* c*

	* Gen empty matrizes
	matrix coefs_I = J(1,$K,.)	
	matrix CI_l_I = J(1,$K,.)
	matrix CI_u_I = J(1,$K,.)
		matrix coefs_II = J(1,$K_II,.)	
		matrix CI_l_II = J(1,$K_II,.)
		matrix CI_u_II = J(1,$K_II,.)
		
	* Store values for mean, percentile CI(2.5) & CI(9.75)	
	local i 0	
	foreach var in $indepI {
	local ++ i
		qui sum `var'
		matrix coefs_I[1, `i'] = r(mean)
			 centile `var', centile(2.5 97.5)
			matrix CI_l_I[1, `i'] = r(c_1)
			matrix CI_u_I[1, `i'] = r(c_2)
   }			
	matrix CI_I = (CI_l_I\CI_u_I)	
	
		local i 0	
		foreach var in $indepII {
		local ++ i
			qui sum `var'
			matrix coefs_II[1, `i'] = r(mean)
				centile `var', centile(2.5 97.5)
				matrix CI_l_II[1, `i'] = r(c_1)
				matrix CI_u_II[1, `i'] = r(c_2)
		}			
	matrix CI_II = (CI_l_II\CI_u_II)	
	
	
	* Rename cols as variables
	matrix colnames coefs_I = RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013
	matrix colnames CI_I = RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013
	
	matrix colnames coefs_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij		
	matrix colnames CI_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij
* Check data
clear
mat li coefs_I
mat li coefs_II

mat coefs_I_tsfe = coefs_I'
mat CI_I_tsfe = CI_I'
mat coefs_II_tsfe = coefs_II'
mat CI_II_tsfe = CI_II'

clear
svmat coefs_I_tsfe	
svmat CI_I_tsfe			
svmat coefs_II_tsfe			
svmat CI_II_tsfe				

	quietly generate str varnames = ""
local varnames = "RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013"
forvalues i = 1/`=_N' {
    local varname : word `i' of `varnames'
    quietly replace varnames = "`varname'" in `i'
}

	quietly generate str varnames2 = ""
local varnames = "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij "
forvalues i = 1/`=_N' {
    local varname : word `i' of `varnames'
    quietly replace varnames2 = "`varname'" in `i'
}
order varnames*

*save Estimates/EU_integration/EU_integration_coefficients.dta, replace
save Estimates/EU_integration/EU_integration_coefficients_paper.dta, replace


di as red "done step: save results coefficients"
}