clear

* Input your directory
global directory "C:\Users\Julia\Desktop\EU_integration_effects\Replicationfiles"
cd "$directory"


* Install programs
*ssc install ppmlhdfe
*ssc install ppml_fe_bias
*ssc install outreg
*ssc install gtools
*ssc install rowmat_utils
*ssc install missings
*ssc install spmap
*ssc install shp2dta
*ssc install palettes,replace
*ssc install grstyle, replace

/*******************************************************************************
Data
*******************************************************************************/
qui do Data/make_dummies.do

qui do Data/make_WIOD_long_manu.do


/*******************************************************************************
Counterfactual experiments. Choose "seed" and "number of bootstrap replications"
*******************************************************************************/
* Counterfactual experiment I: EU integration effects (with level)
*do Counterfactual_experiments/EU_integration.do 1234 500
*do Counterfactual_experiments/EU_integration.do 1235 200
*do Counterfactual_experiments/EU_integration.do 1236 300
*do Counterfactual_experiments/EU_integration.do 1237 15
*do Counterfactual_experiments/Robustness/ge_gravity.do 1234 1000

* Counterfactual experiment II: EU integration potential
*do Counterfactual_experiments/EU_trade_potential.do 1234 1000


/*******************************************************************************
 Figures & Tables
*******************************************************************************/
**************************
* Table 1: Data description
**************************
use Data\WIOD_long_manu.dta, clear
	keep if year==1995 | year==1998 | year==2001 | year==2004  | year==2007 | year==2010 | year==2013
	egen ex_id = group(ex)
	egen im_id = group(im)
	global vars = "x ex_id im_id year RTA RTA_noEU schengen_noEU EEA_noEU EU_CHE_noEU EU_TUR_noEU EU15 EU25 EU27 EU28 BORDER lDISTii lDISTij noCOMLANG noCONTIG noCOLONY"
	format $vars %9.3f
	sum $vars, f
	
**************************
* Table 2: Description of groups of trade flows
**************************
use Data\WIOD_long_manu.dta, clear
	keep if year==1995 | year==1998 | year==2001 | year==2004  | year==2007 | year==2010 | year==2013
	table type year, row col
	
**************************
* Figure 1 & Table 3: Coefficient estimates
**************************
*do Counterfactual_experiments\Robustness\comparison_estimates.do 1234 1000
do Figures_and_tables\Figures\coefplot.do
	esttab tsfe twoway_p twoway_cs secondstep_ols,    cell(coefs( fmt(3)) CIL(par(" [") fmt(3)) & CIU(par("" "] "))) incelldelimiter(", ")	mtitle("TSFE (Three-way)" "Two-way (Panel)" "Two-way (Cross-Section)" "OLS (Second step)") scalars("N N(panel)" "obs_II N(cross-section)" "FE_ex_time Exporter-time FE" "FE_im_time Importer-time FE" "FE_pair Pair FE" "FE_ex Exporter FE" "FE_im Importer FE") collabels(none)  
	
**************************		
* Figures 2, Table 4 and 5: GE and Welfare effects by type
**************************
use Estimates/EU_integration/EU_integration_nolevel_type_paper.dta, clear
	table type year, c(m xp) format(%9.1f) 
	table type year, c(m w) format(%9.1f) 
use Estimates/EU_integration/EU_integration_level_type_paper.dta, clear
	table type year, c(m xp) format(%9.1f) 		
	table type year, c(m w) format(%9.1f) 	
keep if year==2013
	qui merge 1:1 type  using Estimates/EU_trade_potential/EU_trade_potential_type_paper.dta
	qui gen xp_potential_wrt_CFL = (( ((xp/100)+1) * ((xp_potential1/100)+1) )-1)*100 if xp ~=. | xp ~=0
	qui gen w_potential_wrt_CFL =  (( ((w/100)+1) * ((w_potential1/100)+1) )-1)*100 if w ~=0
	format xp* w* %9.1f
	list type xp_potential1 xp_potential_wrt_CFL, compress abbreviate(9) table separator(43)  clean 	
	list type w_potential1 w_potential_wrt_CFL, compress abbreviate(9) table separator(43)  clean 	
do Figures_and_tables\Figures\GE_effects_scatterplot.do


**************************	
* Table 6 & 7: GE & welfare effects by country of EU integration (without level)
**************************
*use Estimates/EU_integration/EU_integration_nolevel_country.dta, clear
use Estimates/EU_integration/EU_integration_nolevel_country_paper.dta, clear
	format xp* w* %9.1f
	* sort xp2013
	list country_id xp1* xp2*, compress abbreviate(9) table separator(43) clean 
	list country_id w*, compress abbreviate(9) table separator(43) clean

**************************	
* Table 8 & 9: GE & welfare effects by country of EU integration (with level) and of EU trade potential
**************************
use Estimates/EU_integration/EU_integration_level_country_paper.dta, clear
merge  1:1 country_id using Estimates/EU_trade_potential/EU_trade_potential_country_paper.dta
	gen xp_potential_wrt_BLN = xp_potential1
	gen xp_potential_wrt_CFL = (( ((xp2013/100)+1) * ((xp_potential1/100)+1) )-1)*100 if xp2013 ~=0
		gen w_potential_wrt_BLN = w_potential1
		gen w_potential_wrt_CFL =   (( ((w2013/100)+1) * ((w_potential1/100)+1) )-1)*100 if w2013 ~=0
		format xp* w* %9.1f
		rename *potential1 _*potential
		list country_id xp* , compress abbreviate(9) table separator(43)  clean 
		list country_id w* , compress abbreviate(9) table separator(43)  clean


	
/*******************************************************************************
 Robustness. 
*******************************************************************************/
do Counterfactual_experiments/Robustness/gravity_estimates.do 

esttab two_way three_way phase_in globalization bias exogeneity consecutive, b(%7.3f) se(%7.3f) nogaps 	order(RTA RTA_lag1 RTA_lag2 RTA_lag3 RTA_lag4 RTA_lag5 RTA_lag6 RTA_lag7 RTA_lag8 RTA_lag9 RTA_lag10 RTA_lag11 RTA_lag12 RTA_lead1 RTA_lead2 RTA_lead3 BORDERx* ) drop(_cons BORDER*1996 BORDER*1997 BORDER*1999 BORDER*2000 BORDER*2002 BORDER*2003 BORDER*2005 BORDER*2006 BORDER*2008 BORDER*2009 BORDER*2011 BORDER*2012 BORDER*2014 ) mtitle("Two-way" "Three-way" "Phase-in" "Globalization" "Bias" "Exogeneity" "Consecutive") scalars("test_lag3 Wald test lag 1-3" "test_lag6 Wald test lag 4-6" "test_lag9 Wald test lag 7-9" "test_lag12 Wald test lag 10-12" "test_lead Wald test lead" "joint_size joint size RTA"  "N N") compress staraux

* Robustness: EU specification	
esttab EU EU_homo_II euro euro_II  EU_hetero EU_hetero_II  controls controls_II share share_II , b(%7.3f) se(%7.3f) nogaps 	mtitle( "EU" "EU II" "EURO" "euro II"  "Heterogeneous" "Hetero II" "Controls" "Controls II" "Share" "Share II") order(RTA* EURO EEA* EU_*  schengen* BORDERx* EU15* EUnew*  ) drop(*ex_id* *im_id*)  scalars( "N N") compress staraux 


	
