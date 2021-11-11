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
* Descriptives
********************************************************************************

*Load data
use Data/WIOD_long_manu.dta, clear


*******************/
* 3-year interval
*******************
	keep if year==1995 | year==1998 | ///
			year==2001 | year==2004  | year==2007 | year==2010 | year==2013


egen ex_id=group(ex)
egen im_id=group(im)

* Table 1: Data description
tabstat ex_id im_id year, s(n min max) format(%5.0f) columns(s)
tabstat RTA EU15 EU25 EU27 EU28 BORDER_ij, s(n mean min max) format(%9.3f) columns(s)
	
* Table 2: Description of groups of trade flows
table type year, row col
	
	
	
