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
global sigma = 5
sca number=0

*Load data
use Data/WIOD_long_manu.dta, clear

*******************/
* 3-year interval
*******************
	keep if year==1995 | year==1998 | year==2001 | year==2004  | year==2007 | year==2010 | year==2013
	drop *lead1 *lead2

* Define id's
egen ex_id=group(ex)
egen im_id=group(im)
egen pair_id=group(ex im) 

egen ex_time_id=group(ex year)
egen im_time_id=group(im year)
egen time=group(year)
egen id=group(ex im year)
xtset pair_id time


}

	
* Define variables for baseline and counterfactual scenario
	* First step
	global indep_twoway = "RTA lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER"
	global indep_threeway = "RTA "
	global indep_phase_in = "RTA RTA_lag3 RTA_lag6 RTA_lag9 RTA_lag12"
	global indep_globalization = "RTA RTA_lag3 RTA_lag6 RTA_lag9 RTA_lag12 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013"
	global indep_exogeneity = "RTA RTA_lag3 RTA_lag6 RTA_lag9 RTA_lag12 RTA_lead3 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013"
	global indep_consecutive = "RTA RTA_lag1 RTA_lag2 RTA_lag3 RTA_lag4 RTA_lag5 RTA_lag6 RTA_lag7 RTA_lag8 RTA_lag9 RTA_lag10 RTA_lag11 RTA_lag12 RTA_lead1 RTA_lead2 RTA_lead3 BORDERx1996 BORDERx1997 BORDERx1998 BORDERx1999 BORDERx2000 BORDERx2001 BORDERx2002 BORDERx2003 BORDERx2004 BORDERx2005 BORDERx2006 BORDERx2007 BORDERx2008 BORDERx2009 BORDERx2010 BORDERx2011 BORDERx2012 BORDERx2013 BORDERx2014"
	*global K: word count $indepI
	*di $K

	* Second step
	global indepII = "BORDER EUxBORDER_ij lDISTii lDISTij noCONTIG noCOLONY noCOMLANG "
	global K_II: word count $indepII
	di $K_II
	
				
		
	
/****************************
		RTA with observable variables to proxy bilateral (time-invariant) trade costs
 as in Borchert and Yotov (2016) (eq 2) & estimation as in Yotov et al. (2016) (table 3 col 3) (UN comtrade and cepii tradeprod and indstat)
-> Egger and Nigai & Agnosteva show endogeneity problems 
*******************************/
 ppmlhdfe x $indep_twoway , absorb(i.im_id#i.year i.ex_id#i.year )  cluster(pair) nocons 
		estimates store two_way

/****************************
	 RTA (endogeneity)
as recommended by Baier and Bergstrand (2007)  
&  estimated in Yotov et al. (2016) (table 3 col 4) RTA= 0.557**
-> absorb time invariant trade cost variables
*******************************/
 ppmlhdfe x $indep_threeway , absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons 
 	estimates store three_way
	
/****************************
	Phase-in effects 
as suggested by Baier and Bergstrand (2007),  & Anderson and Yotov (2016) 
&  estimated in Yotov et al. (2014) (table 3 col 5) RTA= 0.291**
&  estimated in Bergstrand et al. (2015): RTA= 0.241**, 
-> comfirm phase-in effects (stat. significant)
*******************************/
 ppmlhdfe x $indep_phase_in , absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons 
	estimates store phase_in

	
/****************************
		Border (inter) 
as in Bergstrand, Larch and Yotov (2015) table 2 : RTA= 0.111*, inter1994= 0.099**, inter2002= 0.292**
& Yotov et al. (2014): table 3 col 7 : RTA= 0.116 (ns)
-> confirm that RTA's effect diminishes
-> confirm decreasing border effect on trade 
*******************************/
 ppmlhdfe x  $indep_globalization , absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair)  
	estimates store globalization

 
/****************************
		Bias Correction 
as in Weidner & Zylkin (2019)
--> small bias
*******************************/

ppmlhdfe x  $indep_globalization , absorb(im_id#year ex_id#year i.im_id#i.ex_id) cluster(im_id#ex_id) nocons d
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
 ppml_fe_bias x $indep_globalization , i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact
			drop lambda
			estadd sca N = obs, replace
			estimates store bias

/****************************
	Shares following Mayer et al. (2019)
*******************************
ppmlhdfe share  $indep_globalization , absorb(im_id#year ex_id#year i.im_id#i.ex_id) cluster(im_id#ex_id) nocons d
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
 ppml_fe_bias share $indep_globalization , i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact
			drop lambda
			estadd sca N = obs, replace
			estimates store share
*/	
/****************************
		RTA lead 
to test for reversed causality as recommended in Wooldridge (2010)
&  estimated in Yotov et al. (2014) (table 3 col 5) RTA_lead= 0.077
&  estimated in Bergstrand et al. (2015): RTA_lead= 0.005, 
-> ok
*******************************/
ppmlhdfe x $indep_exogeneity , absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons d
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
 ppml_fe_bias x  $indep_exogeneity , i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact
			drop lambda
			estadd sca N = obs, replace
	estimates store exogeneity
		
			
/****************************
		Full sample 
as suggested by Egger, Larch and Yotov (2020): Time-interval data may lead to significant biases in the duration and timing of policy effects. The pattern of anticipation effects cannot be captured with interval samples. 
---> Parameter estimates not directly comparable, since the DGP of the underlying variable (RTA) is unknown. Following Greene (fifth edition Fig 19.1) we can compare the total effect for "lead" coefs of the consecutive model with the interval sample. Perform a Wald test to test whether the sum of "lead" coefficients in the consecutive model is zero. Cannot reject the Null. The total effect is similar, only the timing is different. 
*******************************/
preserve
	use Data/WIOD_long_manu.dta, clear
	drop BORDERx*
	* Generate BORDERxYEAR Dummies with controlling the reference dummy
	forvalues i=1996(1)2014{
		gen BORDERx`i'=1 if BORDER==1 & year==`i'
			replace BORDERx`i'=0 if BORDERx`i'==.
	}
	
	egen ex_id=group(ex)
	egen im_id=group(im)
	egen pair_id=group(ex im) 
	
	* Regress
	ppmlhdfe x $indep_consecutive , absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons d
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
	ppml_fe_bias x  $indep_consecutive , i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact
		drop lambda
		estadd sca N = obs, replace	
			* Wald test (test the sum of coefs = 0)  Prob > chi2 less than 0.05?
			test RTA_lead1 + RTA_lead2 + RTA_lead3=0
			estadd scalar  test_lead = r(p), replace
				test RTA_lag1 + RTA_lag2 + RTA_lag3=0
				estadd scalar  test_lag3 = r(p), replace
			test RTA_lag4 + RTA_lag5 + RTA_lag6=0
			estadd scalar  test_lag6 = r(p), replace
				test RTA_lag7 + RTA_lag8 + RTA_lag9=0
				estadd scalar  test_lag9 = r(p), replace
			test RTA_lag10 + RTA_lag11 + RTA_lag12=0
			estadd scalar  test_lag12 = r(p), replace

	estimates store consecutive
restore	

* joint size of RTA plus lags and leads: 
estimates restore exogeneity
	di _b[RTA] + _b[RTA_lead3] + _b[RTA_lag3] + _b[RTA_lag6] + _b[RTA_lag9] + _b[RTA_lag12]
	estadd scalar joint_size = _b[RTA] + _b[RTA_lead3] + _b[RTA_lag3] + _b[RTA_lag6] + _b[RTA_lag9] + _b[RTA_lag12], replace

estimates restore consecutive
	di _b[RTA] + _b[RTA_lead1] + _b[RTA_lead2] + _b[RTA_lead3] + _b[RTA_lag1]  + _b[RTA_lag2]  + _b[RTA_lag3]  + _b[RTA_lag4]  + _b[RTA_lag5] + _b[RTA_lag6]  + _b[RTA_lag7]  + _b[RTA_lag8]  + _b[RTA_lag9]  + _b[RTA_lag10]  + _b[RTA_lag11]  + _b[RTA_lag12]
	estadd scalar joint_size = _b[RTA] + _b[RTA_lead1] + _b[RTA_lead2] + _b[RTA_lead3] + _b[RTA_lag1]  + _b[RTA_lag2]  + _b[RTA_lag3]  + _b[RTA_lag4]  + _b[RTA_lag5] + _b[RTA_lag6]  + _b[RTA_lag7]  + _b[RTA_lag8]  + _b[RTA_lag9]  + _b[RTA_lag10]  + _b[RTA_lag11]  + _b[RTA_lag12], replace

* Table: Robustness RTA
*esttab two_way three_way phase_in  globalization bias exogeneity consecutive, b(%7.3f) se(%7.3f) nogaps 	order(RTA RTA_lag1 RTA_lag2 RTA_lag3 RTA_lag4 RTA_lag5 RTA_lag6 RTA_lag7 RTA_lag8 RTA_lag9 RTA_lag10 RTA_lag11 RTA_lag12 RTA_lead1 RTA_lead2 RTA_lead3 BORDERx* ) drop(_cons BORDER*1996 BORDER*1997 BORDER*1999 BORDER*2000 BORDER*2002 BORDER*2003 BORDER*2005 BORDER*2006 BORDER*2008 BORDER*2009 BORDER*2011 BORDER*2012 BORDER*2014 ) mtitle("Two-way" "Three-way" "Phase-in" "Globalization" "Bias" "Exogeneity" "Consecutive") scalars("test_lag3 Wald test lag 1-3" "test_lag6 Wald test lag 4-6" "test_lag9 Wald test lag 7-9" "test_lag12 Wald test lag 10-12" "test_lead Wald test lead" "joint_size joint size RTA"  "N N")
/*
consecutive: Use consecutive years, show only border effects of every third. 
*/

	
	
	

********************************************************************************
******************** EU ********************************************************
********************************************************************************



	global indep_EU = "RTA RTA_lag3 RTA_lag6 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 BORDERx1998xEU BORDERx2001xEU BORDERx2004xEU BORDERx2007xEU BORDERx2010xEU BORDERx2013xEU"
	
	global indep_EURO = "RTA RTA_lag3 RTA_lag6 EURO BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 BORDERx1998xEU BORDERx2001xEU BORDERx2004xEU BORDERx2007xEU BORDERx2010xEU BORDERx2013xEU"

	
	global indep_EU_hetero = "RTA RTA_lag3 RTA_lag6 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013"
	
xi i.schengen_noEU
		rename _Ischengen* schengen_noEU*
		
	global indep_controls = "RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013"


	global indepII = "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij "
	global K_II: word count $indepII
	di $K_II

	
/*******************************************************************************
Specification EU
*******************************************************************************/
	forvalues i=1998(3)2014{
		gen BORDERx`i'xEU=1 if BORDER==1 & year==`i' & EU==1
			replace BORDERx`i'xEU=0 if BORDERx`i'xEU==. 
			*replace BORDERx`i'=0 if BORDERx`i'xEU==1 // include / exclude!! globalization
	}
	
tsset pair_id year
forvalues year= 1998(3)2013{
	gen BORDERx`year'xEU_lead3 = F3.BORDERx`year'xEU
	replace BORDERx`year'xEU_lead3 = 0 if BORDERx`year'xEU_lead3 ==.  
}
xtset pair_id year

	
* Exclude EU from RTA
replace RTA = RTA_noEU
forvalues i= 3(3)12 {
	replace RTA_lag`i' =RTA_lag`i'_noEU
}	
	*table year , c(sum RTA sum RTA_noEU ) row col



/****************************
EU homogeneous
*****************************/

ppmlhdfe x $indep_EU, absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons d
		predict lambda
        matrix beta = e(b)
 		sca obs= e(N)
ppml_fe_bias x $indep_EU, i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact	
			estadd sca N = obs, replace	
	estimates store EU
	estimates store first_step
		drop lambda
		
	preserve	
		/****************************
		recover pair fixed effects
		****************************/
		predict z_delta_BLN, xb								// linear prediction of regressors without pair fe
		gen double ez_delta_BLN = exp(z_delta_BLN)
			egen double ez_delta_BLN_bar = mean(ez_delta_BLN), by(pair_id)
			egen double x_bar = mean(x), by(pair_id)
			gen double emu_ij = (x_bar)/(ez_delta_BLN_bar)


		**************************************************			
		* Second step: estimate time variant parameters
		**************************************************
		collapse (first) ex im ex_id im_id  RTA schengen* EU* EEA* (mean) emu_ij (max)  ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type, by(pair_id)	
	
		glm emu_ij i.ex_id ib(43).im_id $indepII  , family(poisson) irls nocons vce(robust)
		estimates store EU_homo_II
	restore
		
/****************************
EU heterogeneous
-> globalization captured by ROW (cannot be identified simulaneously with BORDERxYEAR)
*****************************/
	
ppmlhdfe x $indep_EU_hetero , absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons d		
		predict lambda
        matrix beta = e(b)
		sca obs= e(N)
 ppml_fe_bias x $indep_EU_hetero, i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact	
			estadd sca N = obs, replace	
	estimates store EU_hetero
		drop lambda

	preserve	
		/****************************
		recover pair fixed effects
		****************************/
		predict z_delta_BLN, xb								// linear prediction of regressors without pair fe
		gen double ez_delta_BLN = exp(z_delta_BLN)
			egen double ez_delta_BLN_bar = mean(ez_delta_BLN), by(pair_id)
			egen double x_bar = mean(x), by(pair_id)
			gen double emu_ij = (x_bar)/(ez_delta_BLN_bar)


		**************************************************			
		* Second step: estimate time variant parameters
		**************************************************
		collapse (first) ex im ex_id im_id  RTA schengen* EU* EEA* (mean) emu_ij (max)  ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type, by(pair_id)	
	
		glm emu_ij i.ex_id ib(43).im_id $indepII  , family(poisson) irls nocons vce(robust)
		estimates store EU_hetero_II
	restore
	
/****************************
EURO zone
-> no suggestion for deeper integration of EURO countries
*****************************/

 ppmlhdfe x $indep_EURO, absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons d
		predict lambda
        matrix beta = e(b)
 		sca obs= e(N)
ppml_fe_bias x $indep_EURO, i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact	
			estadd sca N = obs, replace	
	estimates store euro
		drop lambda
		
	preserve
		/****************************
		recover pair fixed effects
		****************************/
		predict z_delta_BLN, xb								// linear prediction of regressors without pair fe
		gen double ez_delta_BLN = exp(z_delta_BLN)
			egen double ez_delta_BLN_bar = mean(ez_delta_BLN), by(pair_id)
			egen double x_bar = mean(x), by(pair_id)
			gen double emu_ij = (x_bar)/(ez_delta_BLN_bar)


		**************************************************			
		* Second step: estimate time variant parameters
		**************************************************
		collapse (first) ex im ex_id im_id  RTA schengen* EU* EEA* (mean) emu_ij (max)  ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type, by(pair_id)	
		
		glm emu_ij i.ex_id ib(43).im_id $indepII  , family(poisson) irls nocons vce(robust)
		estimates store euro_II
	restore	
/****************************
Controls - (Schengen, EEA, EU-CHE, EU-TUR) 
-> control for schengen effects outside EU, because not part of EU integration
*****************************/
/* Exclude EU & Schengen from RTA
replace RTA = RTA_nocontrols 
forvalues i= 1(1)12 {
	replace RTA_lag`i' =RTA_lag`i'_nocontrols 
}	
*/

ppmlhdfe x $indep_controls , absorb(i.im_id#i.year i.ex_id#i.year i.im_id#i.ex_id) cluster(pair) nocons d
		predict lambda
        matrix beta = e(b)
			sca obs= e(N)
			global K = e(rank)
ppml_fe_bias x $indep_controls , i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact	
			estadd sca N = obs, replace	
estimates store controls
		drop lambda
		
	preserve	
		/****************************
		recover pair fixed effects
		****************************/
		predict z_delta_BLN, xb								// linear prediction of regressors without pair fe
		gen double ez_delta_BLN = exp(z_delta_BLN)
			egen double ez_delta_BLN_bar = mean(ez_delta_BLN), by(pair_id)
			egen double x_bar = mean(x), by(pair_id)
			gen double emu_ij = (x_bar)/(ez_delta_BLN_bar)


		**************************************************			
		* Second step: estimate time variant parameters
		**************************************************
		collapse (first) ex im ex_id im_id  RTA schengen* EU* EEA* (mean) emu_ij (max)  ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type, by(pair_id)	
		
		glm emu_ij i.ex_id ib(43).im_id $indepII  , family(poisson) irls nocons vce(robust)
		estimates store controls_II
	restore
	
/****************************
	Shares following Mayer et al. (2019)
	-> similar coefs
*******************************/
ppmlhdfe share  $indep_controls , absorb(im_id#year ex_id#year i.im_id#i.ex_id) cluster(im_id#ex_id) nocons d
		predict lambda
        matrix beta = e(b)
		sca obs=e(N)
 ppml_fe_bias share $indep_controls , i(ex_id) j(im_id) t(year) lambda(lambda) beta(beta) exact
			drop lambda
			estadd sca N = obs, replace
			estimates store share
	preserve	
		/****************************
		recover pair fixed effects
		****************************/
		predict z_delta_BLN, xb								// linear prediction of regressors without pair fe
		gen double ez_delta_BLN = exp(z_delta_BLN)
			egen double ez_delta_BLN_bar = mean(ez_delta_BLN), by(pair_id)
			egen double x_bar = mean(x), by(pair_id)
			gen double emu_ij = (x_bar)/(ez_delta_BLN_bar)


		**************************************************			
		* Second step: estimate time variant parameters
		**************************************************
		collapse (first) ex im ex_id im_id  RTA schengen* EU* EEA* (mean) emu_ij (max)  ///
		lDIST* noCONTIG* noCOLONY* noCOMLANG* LANDLO* BORDER*  type, by(pair_id)	
		
		glm emu_ij i.ex_id ib(43).im_id $indepII  , family(poisson) irls nocons vce(robust)
		estimates store share_II
	restore
	
	
	

