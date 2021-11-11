*use "Estimates/EU_integration/EU_integration_coefficients.dta", clear
*merge 1:1 _n  using "Estimates/Robustness/comparison_estimates.dta", force

use "Estimates/EU_integration/EU_integration_coefficients_paper.dta", clear
merge 1:1 _n  using "Estimates/Robustness/comparison_estimates_paper.dta"
mat drop _all

rename CI*1 CI*_l
rename CI*2 CI*_u
rename *1 *

*estimates save "Estimates/Robustness/main_estimates.dta", replace

* Define independent variables 		
	global indepI = "RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013"
	global K: word count $indepI
	di $K

	global indepII = "lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij"
	global K_CS: word count $indepII
	di $K_CS
	
	local varlist_I = "coefs_I_tsfe coefs_I_tw CI_I_tsfe_l  CI_I_tw_l CI_I_tsfe_u  CI_I_tw_u"
	local varlist_II ="coefs_II_tsfe  coefs_II_tw coefs_II_cs coefs_II_ols CI_II_tsfe_l  CI_II_tw_l  CI_II_cs_l CI_II_ols_l CI_II_tsfe_u  CI_II_tw_u  CI_II_cs_u CI_II_ols_u"
	
foreach var of varlist `varlist_I'	{
	mkmat `var', matrix(`var')
	matrix `var' = `var''
	matrix colnames `var' = RTA RTA_lag3 RTA_lag6 EEA_noEU EEA_noEU_lag3 EEA_noEU_lag6 EU_CHE_noEU EU_CHE_noEU_lag3 EU_CHE_noEU_lag6 EU_TUR_noEU EU_TUR_noEU_lag3 EU_TUR_noEU_lag6 schengen_noEU__1 schengen_noEU__2 schengen_noEU__3 schengen_noEU__4 schengen_noEU__5 BORDERx1998 BORDERx2001 BORDERx2004 BORDERx2007 BORDERx2010 BORDERx2013 EU15x1998 EU15x2001 EU15x2004 EU15x2007 EU15x2010 EU15x2013 EUnewx2004 EUnewx2007 EUnewx2010 EUnewx2013 EUnewxEU15x2004 EUnewxEU15x2007 EUnewxEU15x2010 EUnewxEU15x2013
}

foreach var of varlist `varlist_II'	{
	mkmat `var', matrix(`var')
	matrix `var' = `var'[1..$K_CS,.]
	matrix `var' = `var''
	matrix colnames `var' = lDISTii lDISTij noCONTIG noCOLONY noCOMLANG BORDER EUxBORDER_ij	
}


matrix CI_I_tsfe = (CI_I_tsfe_l\CI_I_tsfe_u)
matrix CI_I_tw = (CI_I_tw_l\CI_I_tw_u)
	matrix CI_II_tsfe = (CI_II_tsfe_l\CI_II_tsfe_u)
	matrix CI_II_tw = (CI_II_tw_l\CI_II_tw_u)
	matrix CI_II_cs = (CI_II_cs_l\CI_II_cs_u)
	matrix CI_II_ols = (CI_II_ols_l\CI_II_ols_u)


matrix coefs_tsfe = (coefs_I_tsfe, coefs_II_tsfe)	
matrix CI_tsfe = (CI_I_tsfe, CI_II_tsfe)	
matrix coefs_tw = (coefs_I_tw, coefs_II_tw)	
matrix CI_tw = (CI_I_tw, CI_II_tw)	
			
	**************************************************
	* Coefplot combined
	**************************************************
	grstyle clear
	grstyle init
	grstyle set plain, grid
	grstyle set graphsize 15cm 12cm

*colorpalette: FloralWhite dkorange black, ipolate(15)
*colorpalette FloralWhite dkorange black, ipolate(15)  stylefiles(,replace)
	grstyle color background white	
	*grstyle anglestyle vertical_tick 30
grstyle linewidth major_grid vthin
grstyle color major_grid dimgray
grstyle set color black

	grstyle set linewidth .3pt: axisline
	grstyle set linewidth .3pt: xyline
	grstyle color xyline black	
	grstyle set size 3.5pt: tick_label 
	grstyle set linewidth .4pt: pmark  p
	
grstyle set symbolsize 3,pt
grstyle set symbol circle


grstyle color p1markfill  "ebblue%80"
grstyle color p1markline  "black"
grstyle color p2markfill  "238 180 121"
grstyle color p2markline  "gray"
grstyle color p3markfill  "orange_red%80"
grstyle color p3markline  "gray"
grstyle color p4markfill  "magenta%80"
grstyle color p4markline  "gray"

grstyle set legend 11, inside 


	coefplot (matrix(coefs_tsfe), label("TSFE")   ci(CI_tsfe) ciopts(recast(pcspike)) ) ///    
			 (matrix(coefs_tw),  label("Two-way FE (Panel)" ) ci(CI_tw) ciopts(recast(pcspike) lcolor(gray)) msymbol(D) msize(2pt)) ///    
	         (matrix(coefs_II_cs), label("Two-way FE (Cross-Section)" )   ci(CI_II_cs) ciopts(recast(pcspike) lcolor(gray)) msymbol(S) msize(2pt) ) ///
			 			 (matrix(coefs_II_ols),  label("OLS (Second step)" ) ci(CI_II_ols) ciopts(recast(pcspike) lcolor(gray)) msymbol(T) msize(2pt)) ,    ///
		 byopts(rows(1)  graphregion(fcolorwhite) lcolor(white))      	 ///
		 keep(BORDER* EUnew* EU15* EUxBORDER_ij ) ///
		  coeflabels(  BORDERx1998="BORDER{subscript:ij,1998}" ///
		 BORDERx2001="BORDER{subscript:ij,2001}" ///
		 BORDERx2004="BORDER{subscript:ij,2004}" ///
		 BORDERx2007="BORDER{subscript:ij,2007}" ///		
		 BORDERx2010="BORDER{subscript:ij,2010}" ///		
		 BORDERx2013="BORDER{subscript:ij,2013}" ///	
			 EU15x1998="BORDER{subscript:ij,1998} x EU15_EU15{subscript:ij,1998}" ///
			 EU15x2001="BORDER{subscript:ij,2001} x EU15_EU15{subscript:ij,2001}" ///
			 EU15x2004="BORDER{subscript:ij,2004} x EU15_EU15{subscript:ij,2004}" ///
			 EU15x2007="BORDER{subscript:ij,2007} x EU15_EU15{subscript:ij,2007}" ///		
			 EU15x2010="BORDER{subscript:ij,2010} x EU15_EU15{subscript:ij,2010}" ///		
			 EU15x2013="BORDER{subscript:ij,2013} x EU15_EU15{subscript:ij,2013}" ///	
			    EUnewx2004="BORDER{subscript:ij,2004} x EUnew_EUnew{subscript:ij,2004}" ///
			    EUnewx2007="BORDER{subscript:ij,2007} x EUnew_EUnew{subscript:ij,2007}" ///
			    EUnewx2010="BORDER{subscript:ij,2010} x EUnew_EUnew{subscript:ij,2010}" ///
			    EUnewx2013="BORDER{subscript:ij,2013} x EUnew_EUnew{subscript:ij,2013}"  ///
			 EUnewxEU15x2004="BORDER{subscript:ij,2004} x EUnew_EU15{subscript:ij,2004}" ///
			 EUnewxEU15x2007="BORDER{subscript:ij,2007} x EUnew_EU15{subscript:ij,2007}" ///
			 EUnewxEU15x2010="BORDER{subscript:ij,2010} x EUnew_EU15{subscript:ij,2010}" ///
			 EUnewxEU15x2013="BORDER{subscript:ij,2013} x EUnew_EU15{subscript:ij,2013}"  ///
		BORDER="BORDER{subscript:ij} " ///
		EUxBORDER_ij="BORDER{subscript:ij} x EU15_EU15{subscript:ij} ", labsize(tiny))	///
		headings( EU15x1998 = "{bf:deeper EU integration effects}" BORDERx1998 = "{bf:Average integration effects}"  BORDER = "{bf:Initial integration effects}" , labsize(tiny)) ///
		  xlabel(-3(1)1.5, labsize(vsmall)) legend(size(vsmall) cols(1)) ///  ///
		xline(0, lpattern(longdash))   ///
		text( 10 -5.5 "time-variant", orientation(vertical) size(vsmall)) ///
		text( 12 -5.2 "{subscript:|}{superscript:_____________________________________________________________________________________}{subscript:|}", orientation(vertical) ) ///
		text( 26 -5.5 "time-invariant", orientation(vertical) size(vsmall)) ///
		text( 26 -5.2 "{subscript:|}{superscript:___________}{subscript:|}", orientation(vertical)) graphregion(margin(l+7))

	*graph export Figures_and_tables\Figures\coefficients.pdf, as(pdf) replace
	*graph export "C:\Users\Julia\Desktop\EU_integration_effects\Review I\Revision I\Figures\coefficients.pdf", as(pdf) replace
	graph export Figures_and_tables\Figures\coefficients_paper.pdf, as(pdf) replace
	*graph export "C:\Users\Julia\Desktop\EU_integration_effects\Review I\Revision I\Figures\coefficients_paper.pdf", as(pdf) replace
	
	
qui reg coefs_I_tsfe coefs_I_tsfe, vce(bootstrap)
estimates use Estimates/Robustness/main_estimates.dta

estimates store tsfe
estadd matrix coefs = coefs_tsfe, replace
estadd matrix CIL = CI_tsfe[1,1...], replace	
estadd matrix CIU = CI_tsfe[2,1...], replace	
estadd scalar  N = 12455, replace
estadd scalar  obs_II = 1849, replace
estadd scalar  FE_ex_time = 1, replace
estadd scalar  FE_im_time = 1, replace
estadd scalar  FE_pair = 1, replace
estadd scalar  FE_ex = ., replace
estadd scalar  FE_im = ., replace

	estimates store twoway_p
	estadd matrix coefs = coefs_tw, replace
	estadd matrix CIL = CI_tw[1,1...], replace	
	estadd matrix CIU = CI_tw[2,1...], replace	
	estadd scalar  N = 12455, replace
	estadd scalar  obs_II = ., replace
	estadd scalar  FE_ex_time = 1, replace
	estadd scalar  FE_im_time = 1, replace
	estadd scalar  FE_pair = ., replace
	estadd scalar  FE_ex = ., replace
	estadd scalar  FE_im = ., replace

estimates store twoway_cs
estadd matrix coefs = coefs_II_cs, replace
estadd matrix CIL = CI_II_cs[1,1...], replace	
estadd matrix CIU = CI_II_cs[2,1...], replace	
estadd scalar  N = ., replace
estadd scalar  obs_II = 1600, replace
estadd scalar  FE_ex_time = ., replace
estadd scalar  FE_im_time = ., replace
estadd scalar  FE_pair = ., replace
estadd scalar  FE_ex = 1, replace
estadd scalar  FE_im = 1, replace

	estimates store secondstep_ols
	estadd matrix coefs = coefs_II_ols, replace
	estadd matrix CIL = CI_II_ols[1,1...], replace	
	estadd matrix CIU = CI_II_ols[2,1...], replace	
	estadd scalar  N = ., replace
	estadd scalar  obs_II = 1849, replace
	estadd scalar  FE_ex_time = ., replace
	estadd scalar  FE_im_time = ., replace
	estadd scalar  FE_pair = ., replace
	estadd scalar  FE_ex = 1, replace
	estadd scalar  FE_im = 1, replace


