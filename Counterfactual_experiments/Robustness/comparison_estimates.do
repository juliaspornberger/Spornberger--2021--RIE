args seed Breps
local Breps  1000
local seed 1234

qui{

clear all
*matrix drop _all
*timer clear
clear mata
capture set matsize 10000 
capture log close 
capture set more off
program drop _all
sca drop _all

* Chosen arguments
global Breps = `Breps'
global seed = `seed'
set seed $seed
global sigma = 5
sca number=0

*Load data
use Data/WIOD_long_manu.dta, clear



*******************/
* 3-year interval
*******************
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

qui{	
	
* Define variables for baseline and counterfactual scenario
	* First step
		
xi i.schengen_noEU
		rename _Ischengen* schengen_noEU*
		
	global indepI = "RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013"
	global K: word count $indepI
	di $K
			

	* Second step (CONTROLS explaining endogeneity of RTA (EU): lDISTii lDISTij noCONTIG noCOLONY noCOMLANG)
	global indepII = "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij "
	global K_II: word count $indepII
	di $K_II


	
		
*******************
* Dummies for conditional estimator
*******************
	*Choose reference country
		qui tab im
		global N_j=r(r)

			* Origin-time & destination-time FE
			desmat  im_id.year=ind($N_j) ex_id.year=ind($N_j)
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
		
		* (N-1) pair FE for 2) step
		desmat ex_id im_id=ind($N_j), full  
				global exporter_fe = "$term1"
				global importer_fe = "$term2"
			rename _x_# _x_#, renumber 
		unab vars : _x*
		global FE_II  `: word count `vars''
		di $FE_II
				

* Exclude EU from RTA
replace RTA = RTA_noEU 
forvalues i= 1(1)12 {
	replace RTA_lag`i' =RTA_lag`i'_noEU
}	

	
	save tempfile_0.dta, replace

	
/*******************************************************************************
********* TSFE: PPML and OLS ***************************************************
*******************************************************************************/
display as red "start step: TSFE -  PPML and OLS"	

*Load data
use tempfile_0.dta, clear

	/**********************		
	 - First step: estimate time-variant parameters
	**********************/		
	 qui ppmlhdfe x _I_* $indepI, absorb(i.im_id#i.ex_id) cluster(pair) nocons d(ln_pairFE)
		estimates store first_step
		matrix parameters_initial = e(b)
		predict xp_BLN, mu
		
		** Bias Correction **
		qui ppmlhdfe x  $indepI, absorb(im_id#year ex_id#year i.im_id#i.ex_id) cluster(im_id#ex_id) nocons d
		sca obs=e(N)
		predict lambda
        matrix beta = e(b)
		
        qui ppml_fe_bias x  $indepI, i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) 
		estadd sca N = obs, replace
		estimates store first_step_unbiased
		
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

program define bootprogram4, eclass

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

		**************************************************			
		* Second step: estimate time variant parameters
		**************************************************
		save tempfile_2.dta, replace
		sort   pair_id year
			
		collapse (first) ex im ex_id im_id  RTA schengen* EU* EEA* (mean) emu_ij (max) _x_*   ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type, by(pair_id)	
		gen mu_ij= ln(emu_ij)
	
		reg mu_ij _x* $indepII  ,  nocons vce(robust)
			sca obs=e(N)
			estadd sca N = obs, replace
			estimates store second_step_ols
		*	mata t_CFL = t_CFL + d_mu
		*	mata st_store(., "t_CFL", (t_CFL))
			* Store second step estimates
			mat coef_second_step = e(b)
			mat coef_second_step = coef_second_step[1,$FE_II+1..$FE_II+$K_II]
			* Rename columns 
			local colnames_II
			local colnames_II `colnames_II' "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij "
			matrix colnames coef_second_step = `colnames_II'
			matrix list coef_second_step

		display as result "done step: second step"	
	

		save tempfile_3.dta, replace

		**************************************************			
		* Store coefficient estimates
		**************************************************

		matrix results = (draw_parameters, coef_second_step)
		
		* Store in ereturn b (drops all other values in e)
		ereturn post results

*********************
*********************
*********************
* end Program
*********************
*********************
	restore
end

}
********************************************************************************

simulate _b,  reps($Breps): bootprogram4

	**************************************************
	* Generate matrices with coefs and CI for first & second step
	**************************************************
	missings dropobs
	rename _b_* *

	* Gen empty matrizes
		matrix coefs_II = J(1,$K_II,.)	
		matrix CI_l_II = J(1,$K_II,.)
		matrix CI_u_II = J(1,$K_II,.)
		
	* Store values for mean, percentile CI(2.5) & CI(9.75)	
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
	matrix colnames coefs_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij		
	matrix colnames CI_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij
	
	
* Check data
clear
mat li coefs_II

mat coefs_II_ols = coefs_II'
mat CI_II_ols = CI_II'

display as red "done step: TSFE - PPML and OLS"	



	
*/


/*******************************************************************************
********************** Weights following Mayer et al. (2019) ********************
********************************************************************************
display as red "start step: Shares following Mayer et al. (2019)"	

use tempfile_0.dta, clear

	****************************************************************************
	* Step 1) estimate time-variant and time-invariant parameters using a 
	*         two step procedure
	****************************************************************************
	/**********************		
	 - First step: estimate time-variant parameters
	**********************/		
	 qui ppmlhdfe share _I_* $indepI, absorb(i.im_id#i.ex_id) cluster(pair) nocons d(ln_pairFE)
		estimates store first_step
		matrix parameters_initial = e(b)
		predict xp_BLN, mu
		
		** Bias Correction **
		qui ppmlhdfe share  $indepI, absorb(im_id#year ex_id#year i.im_id#i.ex_id) cluster(im_id#ex_id) nocons d
		sca obs=e(N)
		predict lambda
        matrix beta = e(b)
		
        qui ppml_fe_bias share  $indepI, i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) 
		estadd sca N = obs, replace
		estimates store first_step_unbiased_shares
		
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

program define bootprogram1, eclass

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

		/**************
		* Load mata procedures
		**************/				
		qui do Counterfactual_experiments/Mata/mata_procedures.do

		/**************************************************
		* Save trade costs in mata
		**************************************************/
		*qui do Counterfactual_experiments/Mata/save_trade_costs.do
		display as result "done step: save trade costs"	

		**************************************************			
		* Second step: estimate time variant parameters
		**************************************************
		do Counterfactual_experiments/Mata/second_step.do
			sca obs=e(N)
			estadd sca N = obs, replace
			estimates store second_step_shares
		*	mata t_CFL = t_CFL + d_mu
		*	mata st_store(., "t_CFL", (t_CFL))
			* Store second step estimates
			mat coef_second_step = e(b)
			mat coef_second_step = coef_second_step[1,$FE_II+1..$FE_II+$K_II]
			* Rename columns 
			local colnames_II
			local colnames_II `colnames_II' "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij "
			matrix colnames coef_second_step = `colnames_II'
			matrix list coef_second_step

		display as result "done step: second step"	
	

		save tempfile_3.dta, replace

		**************************************************			
		* Store coefficient estimates
		**************************************************

		matrix results = (draw_parameters, coef_second_step)
		
		* Store in ereturn b (drops all other values in e)
		ereturn post results

*********************
*********************
*********************
* end Program
*********************
*********************
	restore
end

}
********************************************************************************

simulate _b,  reps($Breps): bootprogram1

qui{
	**************************************************
	* Generate matrices with coefs and CI for first & second step
	**************************************************
	missings dropobs
	rename _b_* *

	* Gen empty matrizes
	matrix coefs_I = J(1,$K,.)	
	matrix CI_l_I = J(1,$K,.)
	matrix CI_u_I = J(1,$K,.)
		matrix coefs_II = J(1,$K_CS,.)	
		matrix CI_l_II = J(1,$K_CS,.)
		matrix CI_u_II = J(1,$K_CS,.)
		
	* Store values for mean, percentile CI(2.5) & CI(9.75)	
	local i 0	
	foreach var in $indepI {
	local ++ i
		qui sum `var'
		matrix coefs_I[1, `i'] = r(mean)
			qui centile `var', centile(2.5 97.5)
			matrix CI_l_I[1, `i'] = r(c_1)
			matrix CI_u_I[1, `i'] = r(c_2)
   }			
	matrix CI_I = (CI_l_I\CI_u_I)	
	
		local i 0	
		foreach var in $indepII {
		local ++ i
			qui sum `var'
			matrix coefs_II[1, `i'] = r(mean)
				qui centile `var', centile(2.5 97.5)
				matrix CI_l_II[1, `i'] = r(c_1)
				matrix CI_u_II[1, `i'] = r(c_2)
		}			
	matrix CI_II = (CI_l_II\CI_u_II)	
	
	
	* Rename cols as variables
	matrix colnames coefs_I = RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013
	matrix colnames CI_I =RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013
	
	matrix colnames coefs_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij		
	matrix colnames CI_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij
* Check data
clear
mat li coefs_I
mat li coefs_II

mat coefs_I_share = coefs_I'
mat CI_I_share = CI_I'
mat coefs_II_share = coefs_II'
mat CI_II_share = CI_II'

display as red "done step: Weights following Mayer et al. (2019)"	

* Three-way as suggested in Baier and Bergstrand (2007), Yotov et al. (2016), used in Felbermayr et al. (2018, 2018), Baier et al. 2019
*/
*******************************************************************************
						* Two-way (similar as in Borchert & Yotov 2016)
*******************************************************************************
}

qui{

**** Two-way  *****
		 ppmlhdfe x  $indepI $indepII, absorb(im_id#year ex_id#year) cluster(im_id#ex_id) nocons d
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
			ppml_fe_bias x  $indepI $indepII, i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) twoway
			estadd sca N = obs, replace
			estimates store tw_unbiased
			drop lambda
			
			* Store parameters b (without FE and constant)	
			matrix b_initial = e(b)	
			matrix V_initial = e(V)
			

* Start Program
program define bootprogram2, eclass

	clear mata 
	preserve	
		/**************************************************
		* Random draw from estimated parameter distribution 
		**************************************************/
		* Set seed
		drawnorm $indepI $indepII, means(b_initial) cov(V_initial) n(1)  clear
	
		mkmat _all, matrix(draw_parameters)
				matrix results = (draw_parameters)
						* Store in ereturn b (drops all other values in e)
		ereturn post results


end	
}		
simulate _b,  reps($Breps): bootprogram2
	
qui{
	**************************************************
	* Generate matrices with coefs and CI for first & second step
	**************************************************
	missings dropobs
	rename _b_* *

	* Gen empty matrizes
	matrix coefs_I = J(1,$K,.)	
	matrix CI_l_I = J(1,$K,.)
	matrix CI_u_I = J(1,$K,.)
		matrix coefs_II = J(1,$K_CS,.)	
		matrix CI_l_II = J(1,$K_CS,.)
		matrix CI_u_II = J(1,$K_CS,.)
		
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

mat coefs_I_tw = coefs_I'
mat CI_I_tw = CI_I'
mat coefs_II_tw = coefs_II'
mat CI_II_tw = CI_II'

display as red "done step: Two-way"	
*/
*******************************************************************************
						* Cross-Section similar as in Anderson and Yotov (2010), Anderson et al. (2018) GE PPML, Anderson et al. (2018)
*******************************************************************************
}
qui{
display as red "start step: Cross-Section"	

*Load data
use tempfile_0.dta, clear

**** CS  *****
* assumes unobservable trade costs are related to observables
* Anderson & yotov (2010) (take the midyear)
		ppmlhdfe x  $indepII if year ==1995, absorb(im_id ex_id) cluster(im_id#ex_id) nocons d
			estimates store cs
			sca obs=e(N)
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
			*ppml_fe_bias x $indepII, i(ex_id) j(im_id)  lambda(lambda) beta(beta) twoway
			estadd sca N = obs, replace
			estimates store cs_unbiased
			drop lambda
			
			* Store parameters b (without FE and constant)	
			matrix b_initial = e(b)	
				matrix b_initial = b_initial[1,1..7]		
			matrix V_initial = e(V)

* Start Program
program define bootprogram3, eclass

	clear mata 
	preserve	
		/**************************************************
		* Random draw from estimated parameter distribution 
		**************************************************/
		* Set seed	
		drawnorm  $indepII, means(b_initial) cov(V_initial) n(1)  clear
	
		mkmat _all, matrix(draw_parameters)
				matrix results = (draw_parameters)
						* Store in ereturn b (drops all other values in e)
		ereturn post results


end	
}		
simulate _b,  reps($Breps): bootprogram3

qui{
	**************************************************
	* Generate matrices with coefs and CI for first & second step
	**************************************************
	missings dropobs
	rename _b_* *

	* Gen empty matrizes
		matrix coefs_II = J(1,$K_CS,.)	
		matrix CI_l_II = J(1,$K_CS,.)
		matrix CI_u_II = J(1,$K_CS,.)
		
	* Store values for mean, percentile CI(2.5) & CI(9.75)	
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
	matrix colnames coefs_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij		
	matrix colnames CI_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij
* Check data
clear
mat li coefs_II

mat coefs_II_cs = coefs_II'
mat CI_II_cs = CI_II'

display as red "done step: Cross-Section"	


/*???????******************************************************************************
						* Cross-Section OLS similar as in Anderson and Yotov (2010), Egger and Nigai (2015)
*******************************************************************************
}
qui{
display as red "start step: Cross-Section"	

*Load data
use tempfile_0.dta, clear

**** CS  *****
* assumes unobservable trade costs are related to observables
* Anderson & yotov (2010) (take the midyear)
		reg x _x* $indepII  if year==1995, nocons  vce(robust)
			estimates store ols
			sca obs=e(N)
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
			*ppml_fe_bias x $indepII, i(ex_id) j(im_id)  lambda(lambda) beta(beta) twoway
			estadd sca N = obs, replace
			estimates store ols_unbiased
			drop lambda
			
			* Store parameters b (without FE and constant)	
			matrix b_initial = e(b)	
				matrix b_initial = b_initial[1,1..7]		
			matrix V_initial = e(V)

* Start Program
program define bootprogram5, eclass

	clear mata 
	preserve	
		/**************************************************
		* Random draw from estimated parameter distribution 
		**************************************************/
		* Set seed	
		drawnorm  $indepII, means(b_initial) cov(V_initial) n(1)  clear
	
		mkmat _all, matrix(draw_parameters)
				matrix results = (draw_parameters)
						* Store in ereturn b (drops all other values in e)
		ereturn post results


end	
}		
simulate _b,  reps($Breps): bootprogram5

qui{
	**************************************************
	* Generate matrices with coefs and CI for first & second step
	**************************************************
	missings dropobs
	rename _b_* *

	* Gen empty matrizes
		matrix coefs_II = J(1,$K_CS,.)	
		matrix CI_l_II = J(1,$K_CS,.)
		matrix CI_u_II = J(1,$K_CS,.)
		
	* Store values for mean, percentile CI(2.5) & CI(9.75)	
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
	matrix colnames coefs_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij		
	matrix colnames CI_II = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij
* Check data
clear
mat li coefs_II

mat coefs_II_ols = coefs_II'
mat CI_II_ols = CI_II'

display as red "done step: Cross-Section"	
*/



***** Store all estimates
clear
	svmat coefs_I_tw		
	svmat CI_I_tw		
	svmat coefs_II_tw		
	svmat CI_II_tw		
svmat coefs_II_cs		
svmat CI_II_cs		
	svmat coefs_II_ols		
	svmat CI_II_ols		

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
save "Estimates\Robustness\comparison_estimates_paper.dta", replace
	
display as red "done step: store"	

}

