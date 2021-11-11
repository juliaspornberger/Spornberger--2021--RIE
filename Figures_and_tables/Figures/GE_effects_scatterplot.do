***********
*Plot Style
***********
	grstyle clear
	grstyle init
	grstyle set plain, grid
	*grstyle set graphsize 7cm 15cm

*colorpalette: FloralWhite dkorange black, ipolate(15)
*colorpalette FloralWhite dkorange black, ipolate(15)  stylefiles(,replace)
	grstyle color background white	
	grstyle anglestyle vertical_tick 30
	grstyle linewidth major_grid vthin
	grstyle color major_grid dimgray

	grstyle set linewidth .3pt: axisline xyline
	grstyle color xyline black	
	grstyle set size 6pt: tick_label title  heading subheading axis_title
	grstyle set linewidth .4pt: pmark  p
	
grstyle set symbolsize 1,pt
grstyle set symbol circle_hollow

grstyle color p1 "edkblue%50"
grstyle color p2 "ebblue%80"
grstyle color p3 "orange_red%80"
grstyle color p4 "edkblue%80"
grstyle color p5 "ebblue"
grstyle color p6 "midblue%80"
grstyle color p7 "orange%80"
grstyle color p8 "orange_red%80"
grstyle color p9 "edkblue%50"
grstyle color p10 "ebblue%80"
grstyle color p11 "orange_red%80"
grstyle set linewidth 0pt: parea
	grstyle color p12area "edkblue%10"
	grstyle color p13area "ebblue%10"
	grstyle color p14area "orange_red%10"
	grstyle color p15area "edkblue%10"
	grstyle color p16area "ebblue%10"
	grstyle color p17area "midblue%10"
	grstyle color p18area "orange%10"
	grstyle color p19area "orange_red%10"
	grstyle color p20area "edkblue%10"
	grstyle color p21area "ebblue%10"
	grstyle color p22area "orange_red%10"

grstyle set size 5pt : plabel 

/*******************************************************************************
Figure 3: GE trade changes of EU integration without level
*******************************************************************************/
*use Estimates/EU_integration/EU_integration_nolevel_type.dta, clear
use Estimates/EU_integration/EU_integration_nolevel_type_paper.dta, clear
gen level = 0
append using "Estimates/EU_integration/EU_integration_level_type_paper.dta"
replace level = 1 if level ==.
rename xp res1
rename xp_CI_l res_CI_l1
rename xp_CI_u res_CI_u1
rename w res2
rename w_CI_l res_CI_l2
rename w_CI_u res_CI_u2

egen id = group(level type year )

reshape long res res_CI_l res_CI_u , i(id) j(result)
sort result level type year 
drop if result ==2 & type >3 
drop id 

forvalues i = 4(1)8{
	replace type=10+`i' if result==2 & type_w ==`i'
}

foreach var of varlist res* {
	replace `var'=. if `var' ==0
}

	*Define labels	
	replace type=10 if type==7 & level==0 & year==2013
		label define type_label 1 "EUnew (domestic)" 2 "EU15 (domestic)" ///
								3 "ROW (domestic)" 4 "EUnew {&harr} EUnew" ///
 								5 "EU15 {&harr} EU15" 6 "EU15 {&harr} EUnew" ///
								7 "EU {&harr} ROW" 10 "EU {&harr} ROW, ROW {&harr} ROW" 8 " ROW {&harr} ROW" 14 "EUnew" 15 "EU15" 18 "ROW", modify 
		label values type type_label
		
		replace level = 3 if level ==0 & result==2
		replace level = 4 if level ==1 & result==2
		label define level_label 0 "a) without level" 1 "b) with level" 3"c) without level" 4 "d) with level", modify
		label values level level_label

		label define result_label 1 "trade" 2 "welfare" , modify
		label values result result_label

			
	* Connect lines with reference year 2001
	forvalues level=0(1)1{
	   insobs 3, before(1)
	   replace level = `level' in 1/3
	   replace year = 2001 in 1/3 
	   replace result = 1 in 1/3 
			replace type = 1 in 1 if level==`level'
			replace type = 4 in 2 if level==`level'
			replace type = 6 in 3 if level==`level'
		qui sum res if type==3 & year==2001 & result==1 & level==`level'
		replace res = r(mean)  if  year==2001 & res==. & type==1 & level==`level'
		qui sum res if type==8 & year==2001 & result==1 & level==`level'
		replace res = r(mean)  if  year==2001 & res==. & type==4 & level==`level'
		qui sum res if type==7 & year==2001 & result==1 & level==`level'
		replace res = r(mean)  if  year==2001 & res==. & type==6 & level==`level'
			qui sum res_CI_l if type==3 & year==2001 & result==1 & level==`level'
			replace res_CI_l = r(mean)  if  year==2001 & res_CI_l==. & type==1 & level==`level'
			qui sum res_CI_l if type==8 & year==2001 & result==1 & level==`level'
			replace res_CI_l = r(mean)  if  year==2001 & res_CI_l==. & type==4 & level==`level'
			qui sum res_CI_l if type==7 & year==2001 & result==1 & level==`level'
			replace res_CI_l = r(mean)  if  year==2001 & res_CI_l==. & type==6 & level==`level'
		qui sum res_CI_u if type==3 & year==2001 & result==1 & level==`level'
		replace res_CI_u = r(mean)  if  year==2001 & res_CI_u==. & type==1 & level==`level'
		qui sum res_CI_u if type==8 & year==2001 & result==1 & level==`level'
		replace res_CI_u = r(mean)  if  year==2001 & res_CI_u==. & type==4 & level==`level'
		qui sum res_CI_u if type==7 & year==2001 & result==1 & level==`level'
		replace res_CI_u = r(mean)  if  year==2001 & res_CI_u==. & type==6 & level==`level'
	}
		
	forvalues level=3(1)4{
	   insobs 1, before(1)
	   replace level = `level' in 1
	   replace year = 2001 in 1 
	   replace result = 2 in 1 
			replace type = 14 in 1 if level==`level'
		qui sum res if type==18 & year==2001 & result==2 & level==`level'
		replace res = r(mean)  if  year==2001 & res==. & type==14 & level==`level'
			qui sum res_CI_l if type==18 & year==2001 & result==2 & level==`level'
			replace res_CI_l = r(mean)  if  year==2001 & res_CI_l==. & type==14 & level==`level'
			qui sum res_CI_u if type==18 & year==2001 & result==2 & level==`level'
			replace res_CI_u = r(mean)  if  year==2001 & res_CI_u==. & type==14 & level==`level'
	}
	
	* Visually enlarge lower values
	replace res = res *3 if type==1 & result==1 | type==2 & result==1 | type==3 & result==1
	replace res_CI_l = res_CI_l *3 if type==1& result==1 | type==2 & result==1| type==3 & result==1
	replace res_CI_u = res_CI_u *3 if type==1& result==1 | type==2 & result==1| type==3 & result==1

	
	* Duplicate last and first values for labels		
			gen last = .
			gen first = 0 if year ==1995 
			replace first = res if year ==1995 & level==1 | year ==1995 & level==4
	forvalues level=0(1)4{
			replace last = res if year ==2013 & level==`level'
	}	
	
*** Plot **
	tw (scatter res year if type==1, sort lpattern(dash) msymbol(none) c(line)) ///
	   (scatter res year if type==2, sort lpattern(dash) msymbol(none) c(line)) ///
	   (scatter res year if type==3, sort lpattern(dash)  msymbol(none)c(line)) ///
	   (scatter res year if type==4, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter res year if type==5, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter res year if type==6, sort lpattern(solid)  msymbol(none) c(line)) ///
	   (scatter res year if type==7, sort lpattern(solid) msymbol(none) c(line)) ///	  
  	   (scatter res year if type==8, sort lpattern(solid)  msymbol(none) c(line))  ///
  	   (scatter res year if type==14, sort lpattern(solid)  msymbol(none) c(line))  ///
  	   (scatter res year if type==15, sort lpattern(solid)  msymbol(none) c(line))  ///
  	   (scatter res year if type==18, sort lpattern(solid)  msymbol(none) c(line))  ///
		(rarea res_CI_u res_CI_l year if type ==1 ,color(edkblue%10) )		///
		(rarea res_CI_u res_CI_l year if type ==2 , color(edkblue%10))		///
		(rarea res_CI_u res_CI_l year if type ==3 , color(orange_red%10))		///
		(rarea res_CI_u res_CI_l year if type ==4 , color(edkblue%10))		///
		(rarea res_CI_u res_CI_l year if type ==5 , color(ebblue%10))		///
		(rarea res_CI_u res_CI_l year if type ==6 , color(midblue%10))		///
		(rarea res_CI_u res_CI_l year if type ==7 , color(orange%10))		///
		(rarea res_CI_u res_CI_l year if type ==8 , color(orange_red%10))		///
		(rarea res_CI_u res_CI_l year if type ==14 , color(edkblue%10))		///
		(rarea res_CI_u res_CI_l year if type ==15 , color(ebblue%10))		///
		(rarea res_CI_u res_CI_l year if type ==18 , color(orange_red%10))		///
			(scatter last year if type==14, sort lpattern(dash) c(line) msymbol(none) mlabel(type) mlabcolor(edkblue) mlabposition()   ) ///
			(scatter last year if type==15, sort lpattern(dash) c(line) msymbol(none) mlabel(type) mlabcolor(ebblue) mlabposition()  ) ///
			(scatter last year if type==18, sort lpattern(solid) c(line) msymbol(none) mlabel(type) mcolor(orange_red)  mlabcolor(orange_red) msymbol(none)) ///
			(scatter last year if type==1, sort lpattern(dash) c(line) msymbol(none) mlabel(type) mlabcolor(edkblue) mlabposition()   ) ///
			(scatter last year if type==2, sort lpattern(dash) c(line) msymbol(none) mlabel(type) mlabcolor(ebblue) mlabposition()  ) ///
			(scatter last year if type==4, sort lpattern(solid) c(line) msymbol(none) mlabel(type) mcolor(edkblue)  mlabcolor(edkblue) msymbol(none)) ///
			(scatter last year if type==5, sort lpattern(solid) c(line) msymbol(none) mlabel(type) mcolor(ebblue)  mlabcolor(ebblue) msymbol(none) ) ///
			(scatter last year if type==6, sort lpattern(solid) c(line) msymbol(none) mlabel(type) mcolor(midblue)  mlabcolor(midblue) msymbol(none) ) ///
			(scatter last year if type==7  & level ==1, sort lpattern(solid) c(line) msymbol(none) mlabel(type) mcolor(orange)  mlabcolor(orange) msymbol(none)  ) ///
			(scatter last year if type==8 & level ==1, sort lpattern(solid) c(line) msymbol(none) mlabel(type) mcolor(gray)  mlabcolor(orange_red) msymbol(none) mlabposition(1)) ///
			(scatter last year if type==10  & level ==0, sort lpattern(solid) c(line) msymbol(none) mlabel(type) mcolor(orange)  mlabcolor(orange) msymbol(none) mlabposition(3)  ) ///		
		(scatter first year if type==15 & level ==4, sort lpattern(solid) mcolor(ebblue%80))  ///
		(scatter first year if type==18 & level ==4, sort lpattern(solid) mcolor(orange_red%80))  ///
		(scatter first year if type==2 & level ==1, sort lpattern(solid) mcolor(ebblue%80))  ///
		(scatter first year if type==3 & level ==1, sort lpattern(solid) mcolor(orange_red%80))  ///
		(scatter first year if type==5 & level ==1, sort lpattern(solid) mcolor(ebblue%80))  ///
		(scatter first year if type==7 & level ==1, sort lpattern(solid) mcolor(orange%80)) ///
	  (scatter first year if type==1 & level ==0 | type==18 & level ==3, by(level, yrescale   legend(off) note("")   title("{bf:Trade} ") subtitle("{bf:Welfare} ")	)   lpattern(solid) c(line) mcolor(black)    ///
	  		 xscale(range(2001(3)2013))    ///
		 xline(2004 2007 2013, lpattern(longdash) ) ///
		xlabel(1995(6)2013)  ytitle(percent) xtitle("") ///
		note("") ///
		  legend(off) plotregion(margin(5 20 5 5))  ysize(5)  saving("GE_I", replace) ) 

	graph save Figures_and_tables\Figures\GE_I.gph, replace
	
	graph export Figures_and_tables\Figures\GE_I.pdf, as(pdf) replace


	
*/
/*******************************************************************************
Figure 5: GE trade changes of EU integration and EU trade potential
*******************************************************************************/
*use	Estimates/EU_trade_potential/EU_trade_potential_type.dta, clear
use	Estimates/EU_trade_potential/EU_trade_potential_type_paper.dta, clear
	rename *_CI_potential1 *_CI_l	
	rename *_CI_potential2 *_CI_u	
	rename *_potential1 *

global year = 2020
gen year = $year

*append using Estimates/EU_integration/EU_integration_level_type.dta
append using Estimates/EU_integration/EU_integration_level_type_paper.dta
	drop if xp ==0 & year!=2001
	
	
	* Connect lines with reference year 2001
				qui sum xp if type==3 & year==2001
				replace xp = r(mean) if type==1 & year==2001
				qui sum xp if type==8 & year==2001
				replace xp = r(mean) if type==4 & year==2001
				qui sum xp if type==7 & year==2001
				replace xp = r(mean) if type ==6 & year==2001	
					*welfare
					qui sum w if type_w==8 & year==2001
					replace w = r(mean) if type_w==4 & year==2001 & type==1
				
	* Connect CI with references in 2001
			qui sum xp_CI_l if type==3 & year==2001
			replace xp_CI_l = r(mean) if type==1 & year==2001
				qui sum xp_CI_u if type==3 & year==2001
				replace xp_CI_u = r(mean) if type==1 & year==2001
			qui sum xp_CI_l if type==8 & year==2001
			replace xp_CI_l = r(mean) if type==4 & year==2001
				qui sum xp_CI_u if type==8 & year==2001
				replace xp_CI_u = r(mean) if type==4 & year==2001
			qui sum xp_CI_l if type==7 & year==2001
			replace xp_CI_l = r(mean) if type ==6 & year==2001
				qui sum xp_CI_u if type==7 & year==2001
				replace xp_CI_u = r(mean) if type ==6 & year==2001
				*welfare
					qui sum w_CI_l if type_w==8 & year==2001
					replace w_CI_l = r(mean) if type_w==4 & year==2001 & type==1
				qui sum w_CI_u if type_w==8 & year==2001
				replace w_CI_u = r(mean) if type_w==4 & year==2001 & type==1
				
	* Visualize CFL II with reference to CFL I [ (( ((xp%/100)+1) * ((xp_potential1%/100)+1) )-1)*100]
	 gen xp_potential = xp if year == $year
	 gen xp_CI_l_potential = xp_CI_l if year == $year
	 gen xp_CI_u_potential = xp_CI_u if year == $year
	 gen w_potential = w if year == $year
	 gen w_CI_l_potential = w_CI_l if year == $year
	 gen w_CI_u_potential = w_CI_u if year == $year
 
	 
		forvalues type=1(1)8 {
			sum xp if year==2013 & type==`type' 
			replace xp_potential = (( ((r(mean)/100)+1) * ((xp_potential/100)+1) )-1)*100  if type==`type'
				sum xp_CI_l if year==2013 & type==`type' 
				replace xp_CI_l_potential = (( ((r(mean)/100)+1) * ((xp_CI_l_potential/100)+1) )-1)*100  if type==`type'
					sum xp_CI_u if year==2013 & type==`type' 
				replace xp_CI_u_potential = (( ((r(mean)/100)+1) * ((xp_CI_u_potential/100)+1) )-1)*100 if type==`type'
			sum w if year==2013 & type==`type' 
			replace w_potential = (( ((r(mean)/100)+1) * ((w_potential/100)+1) )-1)*100 if type==`type'
				sum w_CI_l if year==2013 & type==`type' 
				replace w_CI_l_potential = (( ((r(mean)/100)+1) * ((w_CI_l_potential/100)+1) )-1)*100  if type==`type'
					sum w_CI_u if year==2013 & type==`type' 
				replace w_CI_u_potential = (( ((r(mean)/100)+1) * ((w_CI_u_potential/100)+1) )-1)*100  if type==`type'
		}
		replace xp_potential = xp if year ==2013 
		replace xp = . if year >2013
		replace xp_CI_l = . if year >2013
		replace xp_CI_u = . if year >2013
			replace w_potential = w if year ==2013
			replace w = . if year >2013
			replace w_CI_l = . if year >2013
			replace w_CI_u = . if year >2013
						
	* Duplicate last value for labels		
			gen last = .
			replace last = xp_potential if year == $year
			gen last_w = . 
			replace last_w = w_potential if year == $year
			gen first = xp if year==1995
			gen first_w = w if year==1995
			
	*Define labels	
		label define type_label 1 "EUnew (domestic)" 2 "EU15 (domestic)" ///
								3 "ROW (domestic)" 4 "EUnew {&harr} EUnew" ///
 								5 "EU15 {&harr} EU15" 6 "EU15 {&harr} EUnew" ///
								7 "EU {&harr} ROW" 8 "ROW {&harr} ROW", modify 
		label values type type_label
		label define type_w_label 4 "EUnew" 5 "EU15" 8 "ROW" , modify 
		label values type_w type_w_label




		*!! label!
		* EUnew di (log(2)/(exp(( 0.197 + 1.289) /18)-1)) ~8
		replace year = 2013+8 if type ==1  & year==$year
		replace year = 2013+8 if type ==4  & year==$year
		* EU15 di (log(2)/(exp(( 0.197 ) /18)-1)) ~63
		replace year = 2013+63-30 if type ==2  & year==$year
		replace year = 2013+63-30 if type ==5  & year==$year
		* EU15-EUnew di (log(2)/(exp(( 0.197 + .9886279) /18)-1)) ~10
		replace year = 2013+10+7 if type ==6  & year==$year // label !!2023
		* ROW .di (log(2)/(exp(( 0.197 ) /18)-1)) ~63
		replace year = 2013+63-30 if type ==3  & year==$year
		replace year = 2013+63-30 if type ==7  & year==$year
		replace year = 2013+63-30 if type ==8  & year==$year
		

*** Plot **
grstyle set size 7pt :  tick_label title  
grstyle set size 9pt : heading subheading axis_title plabel

grstyle color p1 "edkblue%80"
grstyle color p2 "ebblue%80"
grstyle color p3 "orange_red%80"
grstyle color p4 "edkblue%80"
grstyle color p5 "ebblue"
grstyle color p6 "midblue%80"
grstyle color p7 "orange%80"
grstyle color p8 "orange_red%80"
grstyle set linewidth 0pt: parea
	grstyle color p9area "edkblue%10"
	grstyle color p10area "ebblue%10"
	grstyle color p11area "orange_red%10"
	grstyle color p12area "edkblue%10"
	grstyle color p13area "ebblue%10"
	grstyle color p14area "midblue%10"
	grstyle color p15area "orange%10"
	grstyle color p16area "orange_red%10"
grstyle color p17 "edkblue%80"
grstyle color p18 "ebblue%80"
grstyle color p19 "orange_red%80"
grstyle color p20 "edkblue%80"
grstyle color p21 "ebblue"
grstyle color p22 "midblue%80"
grstyle color p23 "orange%80"
grstyle color p24 "orange_red%80"
	grstyle color p25area "edkblue%10"
	grstyle color p26area "ebblue%10"
	grstyle color p27area "orange_red%10"
	grstyle color p28area "edkblue%10"
	grstyle color p29area "ebblue%10"
	grstyle color p30area "midblue%10"
	grstyle color p31area "orange%10"
	grstyle color p32area "orange_red%10"


	tw (scatter xp year if type==1, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter xp year if type==2, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter xp year if type==3, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter xp year if type==4, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter xp year if type==5, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter xp year if type==6, sort lpattern(solid) msymbol(none) c(line)) ///
	   (scatter xp year if type==7, sort lpattern(solid) msymbol(none) c(line)) ///	  
  	   (scatter xp year if type==8, sort lpattern(solid) msymbol(none) c(line)) ///
		(rarea xp_CI_u xp_CI_l  year if type ==1 & year<=2013)		///
		(rarea xp_CI_u xp_CI_l  year if type ==2 & year<=2013)		///
		(rarea xp_CI_u xp_CI_l  year if type ==3 & year<=2013)		///
		(rarea xp_CI_u xp_CI_l  year if type ==4 & year<=2013)		///
		(rarea xp_CI_u xp_CI_l  year if type ==5 & year<=2013)		///
		(rarea xp_CI_u xp_CI_l  year if type ==6 & year<=2013)		///
		(rarea xp_CI_u xp_CI_l  year if type ==7 & year<=2013)		///
		(rarea xp_CI_u xp_CI_l  year if type ==8 & year<=2013)		///
			(scatter xp_potential year if type==1, sort lpattern(dash) c(line) msize(vsmall) lcolor(edkblue%80) mcolor(edkblue%80) ) ///
			(scatter xp_potential year if type==2, sort lpattern(dash) c(line) msize(vsmall) lcolor(ebblue%80) mcolor(ebblue%80) ) ///
			(scatter xp_potential year if type==4, sort lpattern(dash) c(line) msize(vsmall) lcolor(edkblue%80) mcolor(edkblue%80) ) ///
			(scatter xp_potential year if type==5, sort lpattern(dash) c(line) msize(vsmall) lcolor(ebblue%80) mcolor(ebblue%80)) ///
			(scatter xp_potential year if type==6, sort lpattern(dash) c(line) msize(vsmall) lcolor(midblue%80) mcolor(midblue%80)) ///
			(scatter xp_potential year if type==7, sort lpattern(dash)  c(line) msize(vsmall) lcolor(orange%80) mcolor(orange%80) ) ///	  
			(scatter xp_potential year if type==8, sort lpattern(dash)  c(line) msize(vsmall) lcolor(orange_red%80) mcolor(orange_red%80) ) ///	  
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==1 & year>2013, lcolor(edkblue%50) lwidth(vsmall) msymbol(none)  ) ///
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==2 & year>2013, lcolor(ebblue%50) lwidth(vsmall) msymbol(none)  ) ///
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==3 & year>2013,  lcolor(gray%50) lwidth(vsmall) msymbol(none) ) ///
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==4 & year>2013,  lcolor(edkblue%50) lwidth(vsmall) msymbol(none) )		///
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==5 & year>2013,  lcolor(ebblue%50) lwidth(vsmall) msymbol(none) )		///
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==6 & year>2013,  lcolor(midblue%50) lwidth(vsmall) msymbol(none) )		///
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==7 & year>2013,  lcolor(orange%50) lwidth(vsmall) msymbol(none) )		///
		(rcapsym xp_CI_u_potential xp_CI_l_potential year if type==8 & year>2013,  lcolor(orange_red%50) lwidth(vsmall) msymbol(none) )		///
			(scatter last year if type==1, sort lpattern(dash)  msize(vsmall) mlabel(type) mcolor(edkblue%80)  mlabcolor(edkblue%80) msymbol(none)mlabposition(6)) ///
			(scatter last year if type==2, sort lpattern(dash)  msize(vsmall) mlabel(type) mcolor(ebblue%80)  mlabcolor(ebblue%80) msymbol(none) mlabposition(3)) ///
			(scatter last year if type==4, sort lpattern(solid)  msize(vsmall) mlabel(type) mcolor(edkblue%80)  mlabcolor(edkblue%80) msymbol(none) mlabposition(3) ) ///
			(scatter last year if type==5, sort lpattern(solid)  msize(vsmall) mlabel(type) mcolor(ebblue%80)  mlabcolor(ebblue%80) msymbol(none)  mlabposition(1) ) ///
			(scatter last year if type==6, sort lpattern(solid)  msize(vsmall) mlabel(type) mcolor(midblue%80)  mlabcolor(midblue%80) msymbol(none) mlabposition(2)) ///
			(scatter last year if type==7, sort lpattern(solid) msize(vsmall) mlabel(type) mcolor(orange%80)  mlabcolor(orange%80) msymbol(none) mlabposition()) ///
			(scatter last year if type==8, sort lpattern(solid)  msize(vsmall) mlabel(type)  mlabcolor(orange_red%80) msymbol(none) mlabposition(3)) ///
			(scatter first year if type==2, sort lpattern(solid) mcolor(ebblue%80))  ///
		(scatter first year if type==3, sort lpattern(solid) mcolor(orange_red%80))  ///
		(scatter first year if type==5, sort lpattern(solid) mcolor(ebblue%80))  ///
		(scatter first year if type==7, sort lpattern(solid) mcolor(orange%80) ///
	  title("{bf:a) Trade} ")	subtitle("") ///
		 xscale(range(2001(1)2013)) yscale(range(-50 350)) ylabel(0(100)350 -50, labsize(small)) ///
		 ///
		xlabel(1995 2013 2021 2030 2046,  valuelabel labsize(small))  ytitle(percent) graphregion(fcolor(white) lcolor(white)) ///
		note("") ///
		  legend(off) plotregion(margin(5 25 1 5))  saving("xp_potential", replace) ) 
	  
	graph save Figures_and_tables\Figures\GE_II.gph, replace
	graph export Figures_and_tables\Figures\GE_II.pdf, as(pdf) replace
	*graph save Figures_and_tables\Figures\GE_II_paper.gph, replace
	*graph export Figures_and_tables\Figures\GE_II_paper.pdf, as(pdf) replace

*** Plot GE welfare changes **
grstyle color p1 "edkblue%50"
grstyle color p2 "ebblue%80"
grstyle color p3 "orange_red%80"
grstyle set linewidth 0pt: parea
	grstyle color p4area "edkblue%10"
	grstyle color p5area "ebblue%10"
	grstyle color p6area "orange_red%10"

	tw (scatter w year if type_w==4, sort lpattern(solid)  msymbol(none) c(line)) ///
	   (scatter w year if type_w==5, sort lpattern(solid)  msymbol(none) c(line)) ///
	   (scatter w year if type_w==8, sort lpattern(solid)  msymbol(none) c(line)) ///
		(rarea w_CI_u w_CI_l year if type_w ==4 & year<=2013)		///
		(rarea w_CI_u w_CI_l year if type_w ==5 & year<=2013)		///
		(rarea w_CI_u w_CI_l year if type_w ==8 & year<=2013)		///
			(scatter w_potential year if type_w==4, sort lpattern(dash) c(line) msize(vsmall) lcolor(edkblue%50) mcolor(edkblue%50)) ///
			(scatter w_potential year if type_w==5, sort lpattern(dash)  c(line) msize(vsmall) lcolor(ebblue%80) mcolor(ebblue%80) ) ///	  
			(scatter w_potential year if type_w==8, sort lpattern(dash)  c(line) msize(vsmall) lcolor(orange_red%80) mcolor(orange_red%80) ) ///	  
		(rcapsym w_CI_u_potential w_CI_l_potential year if type_w==4 & year>2013, lcolor(edkblue%50) lwidth(vsmall) msymbol(none)  ) ///
		(rcapsym w_CI_u_potential w_CI_l_potential year if type_w==5 & year>2013, lcolor(ebblue%50) lwidth(vsmall) msymbol(none)  ) ///
		(rcapsym w_CI_u_potential w_CI_l_potential year if type_w==8 & year>2013,  lcolor(orange_red%10) lwidth(vsmall) msymbol(none) ) ///
			(scatter last_w year if type_w==4, sort lpattern(solid) c(line) msymbol(none) mlabel(type_w) mcolor(edkblue)  mlabcolor(edkblue) msymbol(none) mlabposition(1)) ///
			(scatter last_w year if type_w==5, sort lpattern(solid) c(line) msymbol(none) mlabel(type_w) mcolor(edkblue)  mlabcolor(ebblue) msymbol(none) mlabposition()) ///
			(scatter last_w year if type_w==8, sort lpattern(solid) c(line) msymbol(none) mlabel(type_w) mcolor(ebblue)  mlabcolor(orange_red) msymbol(none) mlabposition(3) ) ///
		(scatter first_w year if type_w==5, sort lpattern(solid) mcolor(ebblue%80))  ///
		(scatter first_w year if type_w==8, sort lpattern(solid) mcolor(orange_red%80) ///
	  title("{bf: b) Welfare} ")	subtitle("") ///
	  		 xscale(range(2001(3)2013)) yscale(range())  ///
		 ///
		xlabel(1995 2013 2021 2030 2046)  ytitle(percent) graphregion(fcolor(white) lcolor(white)) ///
		note("") ///
		  legend(off) plotregion(margin(5 5 5 5)) saving("W_potential", replace))	
 

graph combine "xp_potential" "W_potential",cols(2) xcommon graphregion(margin(l+10))    ysize(2)
	graph save Figures_and_tables\Figures\GE_potential.gph, replace
	graph export Figures_and_tables\Figures\GE_potential.pdf, as(pdf) replace
	
