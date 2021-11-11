args seed Breps
*local Breps  1
*local seed 123450


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
sca number = 0

*Load data
use Data/WIOD_long_manu.dta, clear


/*******************
* 3-year interval
*******************/
keep if year==1995 | year==1998 | year==2001 | year==2004  | year==2007 | year==2010 | year==2013

/*******************
* Generate id's
*******************/
egen ex_id=group(ex)
egen im_id=group(im)
egen ex_time_id=group(ex year)
egen im_time_id=group(im year)
egen pair_id=group(ex im) 
egen time=group(year)
egen id=group(ex im year)
xtset pair_id time
}

* Define variables for baseline and counterfactual scenario
	* First step	
xi i.schengen_noEU
		rename _Ischengen* schengen_noEU*
		
	global indepI = "RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013"
	global K: word count $indepI
			

	* Second step (CONTROLS explaining endogeneity of RTA (EU): lDISTii lDISTij noCONTIG noCOLONY noCOMLANG)
	global indepII = "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij "
	global K_II: word count $indepII
	di $K_II

qui{		
*******************
* Fixed effects for conditional estimator
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
			
						

	****************************************************************************
	* Step 1) estimate time-variant and time-invariant parameters using a 
	*         two step procedure
	****************************************************************************
	* Exclude EU & Schengen from RTA
	replace RTA = RTA_noEU 
	forvalues i= 1(1)12 {
		replace RTA_lag`i' =RTA_lag`i'_noEU
	}	


	/**********************		
	 - First step: estimate time-variant parameters
	**********************/		
	qui ppmlhdfe x _I_* $indepI, absorb(i.im_id#i.ex_id) cluster(pair) nocons d(ln_pairFE)
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
		*sca number=number+1
		*local run_number=number
		*local seed = $seed`run_number'
	
		drawnorm $indepI, means(b_initial) cov(V_initial) n(1) clear
	
		mkmat _all, matrix(draw_parameters)

		use tempfile_1.dta,clear 
		
		/**************************************************
		* Reparametrization
		**************************************************/
		qui do Counterfactual_experiments/Mata/reparametrization.do	
			mata b = st_matrix("e(b)")'						// drawn parameters 
			mata b = b[$FE+1..$FE+$K,.]	
		display as result "done step: reparametrization"	

		/**************
		* Load mata procedures
		**************/				
		qui do Counterfactual_experiments/Mata/mata_procedures.do

		/**************************************************
		* Save trade costs mata
		**************************************************/
		qui do Counterfactual_experiments/Mata/save_trade_costs_II.do
		display as result "done step: save trade costs"	

		/**************************************************			
		* Second_step
		**************************************************/
		qui do Counterfactual_experiments/Mata/second_step_II.do
			mata t_CFL = t_CFL + d_mu
			mata st_store(., "t_CFL", (t_CFL))
		
		display as result "done step: second step"	
	
		save tempfile_3.dta, replace

		/***********************
		*Results
		************************/
		display as result "start step: solve "	
			qui do Counterfactual_experiments/Mata/mata_solver_II.do
			mata: results(p, pge)
			
		/***********************
		* Store results by country and type
		***********************/
		matrix GE_type_2013 = (d_xp_type_2013)'
		matrix GE_country_2013 = (d_xp_ex_2013)'
		matrix W_type_2013 = (d_W_type_2013)'
		matrix W_country_2013 = (d_W_ex_2013)'
		matrix C = C_2013
		
		* Rename columns trade
		local colnames
			forvalues type = 1/8 {
				local colnames `colnames' " xp_t`=`type'2013 '" 
			}
		matrix colnames GE_type_2013 = `colnames'							
			* Rename columns trade
			local colnames
				forvalues ex = 1/43 {
					local colnames `colnames' " xp_c`=`ex'2013 '" 
				}
			matrix colnames GE_country_2013 = `colnames'
			matrix list GE_country_2013
		* Rename columns welfare
		local colnames
			forvalues type = 1/3 {
				local colnames `colnames' " w_t`=`type'2013 '" 
			}
		matrix colnames W_type_2013 = `colnames'							
			* Rename columns welfare
			local colnames
				forvalues ex = 1/43 {
					local colnames `colnames' "w_c`=`ex'2013 '" 
				}
			matrix colnames W_country_2013 = `colnames'
			* colnames C
			local C_colnames
			forvalues year = 2013(3)2013 {
					local C_colnames `C_colnames' " c`=`year ''" 
			}
			matrix colnames C = `C_colnames'
			
	
		* Store coefficient estimates
		matrix results = (GE_type_2013, GE_country_2013, W_type_2013, W_country_2013, C)
		
		* Store in ereturn b (drops all other values in e)
		ereturn post results

*********************
*********************
*********************
* end Program
*********************
*********************
di as red "done run" " number= " `seed'  ", rep= " `run_number' " of " $Breps
	restore
end

}

simulate _b, reps($Breps): bootprogram

*save Estimates/EU_trade_potential/EU_trade_potential.dta, replace
save Estimates/EU_trade_potential/EU_trade_potential_paper.dta, replace
erase tempfile_1.dta
erase tempfile_2.dta
erase tempfile_3.dta


********************************************************************************
* Organize estimates
********************************************************************************
qui{
****************************************************
* CFL experiment II: estimates by type
****************************************************
*use Estimates/EU_trade_potential/EU_trade_potential.dta, clear
use Estimates/EU_trade_potential/EU_trade_potential_paper.dta, clear

missings dropobs
keep _b_*t* 


*** Rename variables ***
forvalues year=2013(3)2013{
	forvalues type = 1(1)8{
		rename _b_xp_t`type'`year' xp`type'
	}
	forvalues type = 1(1)3{
		rename _b_w_t`type'`year' w`type'
	}
}

* Reshape simulated results
gen sim_id= _n
reshape long xp w , i(sim_id) j(type)
gen year=2015
global types=8

* Generate and store means and confidence intervals 
	* Gen empty matrizes
	matrix xp_mean = J(1,$types,.)	
	matrix xp_CI_l = J(1,$types,.)
	matrix xp_CI_u = J(1,$types,.)	
		matrix w_mean = J(1,$types,.)	
		matrix w_CI_l = J(1,$types,.)
		matrix w_CI_u = J(1,$types,.)
		
	* Store values for mean, percentile CI(2.5) & CI(9.75)	by type
	replace xp =xp+1
	replace w = w+1	
	*trade
	forvalues t = 1(1)$types {
	qui ameans xp if type==`t'
		matrix xp_mean[1, `t'] = r(mean_g)
			centile  xp if type==`t', centile(2.5 97.5)
			matrix xp_CI_l[1, `t'] = r(c_1)
			matrix xp_CI_u[1, `t'] = r(c_2)
   }	
   * welfare
	forvalues t = 1(1)3 {
	qui ameans w if type==`t'
		matrix w_mean[1, `t'] = r(mean_g)
			centile  w if type==`t', centile(2.5 97.5)
			matrix w_CI_l[1, `t'] = r(c_1)
			matrix w_CI_u[1, `t'] = r(c_2)
   }	
  
   * Combine
	mat xp_CI_potential = (xp_CI_l\xp_CI_u)'	
	mat xp_potential = xp_mean'
	mat w_CI_potential = (w_CI_l\w_CI_u)'	
	mat w_potential = w_mean'
	* Generate and store yearly mean
	*mean b_ if b_!=0
	*	mat ave_potential = r(table)
	*	mat ave_potential = ave_potential[1,1]
	
* Load data from matrix into stata
clear
svmat xp_potential
svmat xp_CI_potential
svmat w_potential
svmat w_CI_potential
	* Percentage changes
	foreach v of varlist  _all {
		replace `v' = (`v'-1)*100
	}


gen type = _n
gen type_w = type
	replace type_w = 4 if type_w==1
	replace type_w = 5 if type_w==2
	replace type_w = 8 if type_w==3

	*Define labels	
		label define type_label 1 "EUnew (domestic)" 2 "EU15 (domestic)" ///
								3 "ROW (domestic)" 4 "EUnew" ///
								5 "EU15" 6 "EUnew-EU15, EU15-EUnew" ///
								7 "EU-ROW, ROW-EU" 8 " ROW", modify 
		label values type type_label
		label values type_w type_label



*save Estimates/EU_trade_potential/EU_trade_potential_type.dta, replace
save Estimates/EU_trade_potential/EU_trade_potential_type_paper.dta, replace

di as red "done step: save results by groups of trade flows"

****************************************************
* CFL experiment II: estimates by country
****************************************************
*use Estimates/EU_trade_potential/EU_trade_potential.dta, clear
use Estimates/EU_trade_potential/EU_trade_potential_paper.dta, clear
missings dropobs

*** Rename variables ***
forvalues year = 2013(3)2013{
forvalues country = 1(1)43{
	rename _b_xp_c`country'`year' xp`country'
	rename _b_w_c`country'`year' w`country'
}
}

keep xp* w*
global N = 43

* Reshape simulated results
gen sim_id= _n
reshape long xp w , i(sim_id) j(country)
gen year=2015

* Generate and store means and confidence intervals 
	* Gen empty matrizes
	matrix xp_mean = J(1,$N,.)	
	matrix xp_CI_l = J(1,$N,.)
	matrix xp_CI_u = J(1,$N,.)
		matrix w_mean = J(1,$N,.)	
		matrix w_CI_l = J(1,$N,.)
		matrix w_CI_u = J(1,$N,.)		
		
	* Store values for mean, percentile CI(2.5) & CI(9.75) by country
	replace xp =xp+1
	replace w = w+1	
		
	forvalues c = 1(1)$N {
		qui ameans xp if country==`c'
		matrix xp_mean[1, `c'] = r(mean)
			qui centile  xp if country==`c', centile(2.5 97.5)
			matrix xp_CI_l[1, `c'] = r(c_1)
			matrix xp_CI_u[1, `c'] = r(c_2)
		qui ameans w if country==`c'
		matrix w_mean[1, `c'] = r(mean)
			qui centile  w if country==`c', centile(2.5 97.5)
			matrix w_CI_l[1, `c'] = r(c_1)
			matrix w_CI_u[1, `c'] = r(c_2)
   }
  	qui ameans xp 
 	matrix xp_yearly = r(mean_g)
  	qui ameans w
 	matrix w_yearly = r(mean_g)

	*mat xp_CI_potential = (xp_CI_l\xp_CI_u)'
	mat xp_potential = (xp_mean)'
		*mat w_CI_potential = (w_CI_l\w_CI_u)'
		mat w_potential = (w_mean)'

clear	
svmat xp_potential
svmat w_potential
	* Percentage changes
	foreach v of varlist  _all {
		replace `v' = (`v'-1)*100
	}

	gen country_id = _n

label define country_label 1 "Australia" 2 "Austria" 3 "Belgium" 4 "Bulgaria" 5 "Brazil" 6 "Canada" 7 "Switzerland" 8 "China" 9 "Cyprus" 10 "CzechRepublic" 11 "Germany" 12 "Denmark" 13 "Spain" 14 "Estonia" 15 "Finland" 16 "France" 17 "GreatBritain" 18 "Greece" 19 "Croatia" 20 "Hungary" 21 "India" 22 "Indonesia" 23 "Ireland" 24 "Italy" 25 "Japan" 26 "SouthKorea" 27 "Lithuania" 28 "Luxembourg" 29 "Latvia" 30 "Mexico" 31 "Malta" 32 "Netherlands" 33 "Norway" 34 "Poland" 35 "Portugal" 36 "Romania" 37 "Russia" 38 "Slovakia" 39 "Slovenia" 40 "Sweden" 41 "Turkey" 42 "Taiwan" 43 "USA"   , modify
label values country country_label


*save Estimates/EU_trade_potential/EU_trade_potential_country.dta, replace
save Estimates/EU_trade_potential/EU_trade_potential_country_paper.dta, replace

di as red "done step: save results countrywise"
}




