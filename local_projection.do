cls 
clear all

cd "" 


// saves the stuff printed on the screen into a separate file 
capture log close 
log using lpirf_inflation_expectations.smcl, replace 

import excel "",  firstrow clear

// browse

label variable CAPB "Cyclically Adjusted Primary Balance"
label variable RBEER "Real Broad Effective Exchange Rate"
label variable RNEER "Real Narrow Effective Exchange Rate"
label variable inflation_expectations "Inflation Expectations"
label variable long_rates "Long term government bond yield"
label variable policy_rates "Central bank policy rates"
label variable short_rates "Short term nominal rates"
label variable unemployment_rate "Unemployment Rate"
la var ZLB "Zero Lower Bound"
la var FR "Fixed Exchange Rate Regime"

// describe
// summarize 
// correlate

// creating a new variable to represent countries with numbers
egen country_id = group(Country)


// declaring that the data is panel data 
xtset country_id Year, yearly


// commented out for including Pre-COVID-19 data
// drop if Year > 2019


// xtdescribe 
// xtsum


// table of descriptive statistics
// dtable CAPB unemployment_rate inflation_expectations RBEER debt_GDP


// creating a new variable to represent the lag of real GDP for local projection regression  
// the log was applied on real_GDP during the data cleaning process 
gen log_real_GDP = 100 * real_GDP

/*
The codes below are based on the replication codes provided by Òscar Jordà. 

Sources: 
https://sites.google.com/site/oscarjorda/home/local-projections
https://github.com/ojorda/JEL-Code
*/


// Long differences: real_GDP(t+h) - real_GDP(t−1)
// d henceforth represents difference in variables
forvalues h=0/5 {
    gen D`h'_real_GDP = F`h'.log_real_GDP - L1.log_real_GDP
}


// sum of long differences
forv h=0/5 {
	egen S`h'_real_GDP = rowtotal(D0_real_GDP-D`h'_real_GDP)
}


// fiscal shock though CAPB
gen capb_pp = 100 * CAPB 
gen dtCAPB = capb_pp - L.capb_pp

// flip the sign so increase = expansion rather than tightening
gen dCAPB = -dtCAPB


// HP filters that first blocks the cyclical component of GDP and finds trends
gen y_hpcyc=.
gen y_hptrend=. 
levelsof country_id, local(ctys) 
foreach c in `ctys' {
	disp "`c'"
	gen temp0 = log_real_GDP if country_id==`c'
	gen ifsyear = country_id*10000+Year
	tsset ifsyear
	tsfilter hp temp1 = temp0  , smooth(400) trend(temp2)
	replace y_hpcyc = temp1 if country_id==`c'
	replace y_hptrend = temp2 if country_id==`c'
	
// 	list country_id Year y_hpcyc if temp0~=.

	drop temp0 temp1 temp2 ifsyear

}
xtset country_id Year, yearly 
// summarize dCAPB unemployment_rate inflation_expectations RBEER debt_GDP, detail


// changing the numbers into percentage points 
gen inflation_pct = inflation_expectations * 100
gen unemployment_pct = unemployment_rate * 100 
gen long_rates_pct = long_rates * 100
gen debt_GDP_pct = debt_GDP * 100 
gen policy_rates_pct = policy_rates * 100


gen _x1 = L.D.log_real_GDP			// lagged GDP growth
gen _x2 = L2.D.log_real_GDP
gen _x3 = L.dCAPB					// lagged fiscal shock
gen _x4 = L2.dCAPB
gen _x5 = L.y_hpcyc
gen _x6 = L.unemployment_pct		// current year's unemployment
gen _x7 = L.debt_GDP_pct			
gen _x8 = L.D.RBEER					// real broad exchange rate 
gen _x9 = L.D.long_rates_pct		// long‐term gov bond yield
gen _x10 = inflation_pct
gen _x11 = L.policy_rates_pct		// central bank policy rates


// Horizon-by-horizon OLS
eststo clear
forvalues h = 0/5 {
	xtreg S`h'_real_GDP dCAPB _x*, ///
	fe vce(cluster country_id)
    estimates store OLS`h'
}
	
// Tabulate
esttab OLS* , b(2) se keep(dCAPB*)

correlate dCAPB* _x*


quietly {
    // Slicing by horizon and saving files 
    forvalues h = 0/4 {
        preserve

		gen Y = S`h'_real_GDP
        gen T_h`h' = dCAPB
        gen I = country_id

        forvalues k = 1/11 {
            rename _x`k' X`k'_h`h'
        }

        gen H = `h'
        keep Y T_h`h' X*_h`h' I H Year
        save "stack`h'.dta", replace
        restore
    }

    // Appending stacked files
    preserve
    clear
    forvalues h = 0/4 {
        append using "stack`h'.dta"
    }

    // Replacing missing values with 0 
	// https://www.statalist.org/forums/forum/general-stata-discussion/general/1690218-converting-missing-values-into-0-for-some-columns-with-conditions
    ds T_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
    ds X*_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
    save "stack.dta", replace

    // Stacked OLS for joint inference of estimates 
    gen FE = I*1000 + H
    tsset FE Year
    regress Y T_h0 T_h1 T_h2 T_h3 T_h4 X* i.FE, vce(cluster I)

    estimates store stackedOLS

    lincom (T_h0 + T_h1 + T_h2 + T_h3 + T_h4)/5
    test T_h0 T_h1 T_h2 T_h3 T_h4

   
    forvalues h = 0/4 {
        erase "stack`h'.dta"
    }
}
	
esttab stackedOLS , b(2) se keep(*T_*)
	lincom (T_h0+T_h1+T_h2+T_h3+T_h4)/5


// output summary for Table 1 
esttab OLS* , b(2) se keep(*CAPB*)
esttab stackedOLS , b(2) se keep(*T_*)
lincom (T_h0+T_h1+T_h2+T_h3+T_h4)/5


// p-value  0.0549
test T_h0 T_h1 T_h2 T_h3 T_h4

// p-value 0.0202
test T_h0 T_h1 T_h2 > 0



/* 
Code for exchange regime effects
*/
eststo clear 
forvalues h = 0/5 {
	xtreg  S`h'_real_GDP				///
           c.dCAPB##i.FR  _x*,	///
           fe vce(cluster country_id)
    estimates store FR_OLS`h'
}
	

quietly {
	forvalues h = 0/4 {
		preserve

		gen Y = S`h'_real_GDP
		gen I = country_id
		
		gen T_h`h' = dCAPB          // floating regime
		gen TFR_h`h' = dCAPB * FR      // fixed regime


     
		forvalues k = 1/11 {
			rename _x`k' X`k'_h`h'
        }

		gen H = `h'
		// keep only the variables we need
		keep Y T_h`h' TFR_h`h' X*_h`h' I H FR Year
		save "FR_stack`h'.dta", replace
		restore
	}
  
	// Appending stacked files
    preserve
    clear
    forvalues h = 0/4 {
        append using "FR_stack`h'.dta"
    }

    // Replacing missing values with 0 
	// https://www.statalist.org/forums/forum/general-stata-discussion/general/1690218-converting-missing-values-into-0-for-some-columns-with-conditions
    ds T_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
	
	ds TFR_h*
	foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
	
    ds X*_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
    save "FR_stack.dta", replace
	
	
	use FR_stack.dta, clear
	gen FE = I*1000 + H
	xtset FE Year
	
	xtreg Y T_h* TFR_h* ///
	X*_h*,	///
	fe vce(cluster I)
	

	estimates store FR_stackedOLS

  
    forvalues h = 0/4 {
        erase "FR_stack`h'.dta"
    }
}
	
	
esttab FR_stackedOLS , b(2) se keep(*T_*)

// average multiplier - floating regime 
lincom (T_h0+T_h1+T_h2+T_h3+T_h4)/5    


esttab FR_stackedOLS , b(2) se keep(T_*) 
esttab FR_stackedOLS , b(2) se keep(TFR_*) 


// fixed regime multiplier by horizon 
forvalues h = 0/4 {
    lincom T_h`h' + TFR_h`h', 
}

// average multiplier - fixed regime
lincom ((T_h0+TFR_h0)+(T_h1+TFR_h1)+(T_h2+TFR_h2)+(T_h3+TFR_h3)+(T_h4+TFR_h4))/5 
	


/* 
Now will code for ZLB period with interaction effects 
               /* SCAPB`h'    */                      ///
               /* SCAPB`h'#ZLB*/                      ///
*/	
	

// https://www.statalist.org/forums/forum/general-stata-discussion/general/1582318-interaction-terms-in-regression

// central bank policy rates no longer needed due to ZLB 
drop _x11 



eststo clear
forvalues h = 0/5 {
	xtreg  S`h'_real_GDP                        ///
           c.dCAPB##i.ZLB _x*,              ///
           fe vce(cluster country_id)
    estimates store ZLB_OLS`h'
}

// esttab ZLB_OLS* , b(2) se keep(dCAPB)

quietly {
	forvalues h = 0/4 {
		preserve

		gen Y = S`h'_real_GDP
		gen I = country_id
		
		gen T_h`h' = dCAPB          
		
		gen ZLB_h`h' = dCAPB * ZLB      // ZLB interaction effect
		
		gen ZLB_dummy = ZLB  
		
		
		forvalues k = 1/9 {
			rename _x`k' X`k'_h`h'
        }
		
		gen H = `h'
		// keep only the variables we need
		keep Y T_h`h' ZLB_h`h' X*_h`h' I H Year ZLB_dummy
		save "ZLB_stack`h'.dta", replace
		restore
	}
  
	// Appending stacked files
    preserve
    clear
    forvalues h = 0/4 {
        append using "ZLB_stack`h'.dta"
    }

    // Replacing missing values with 0 
	// https://www.statalist.org/forums/forum/general-stata-discussion/general/1690218-converting-missing-values-into-0-for-some-columns-with-conditions
    ds T_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
	
	ds ZLB_h*
	foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
	
    ds X*_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
    save "ZLB_stack.dta", replace
	
	use ZLB_stack.dta, clear
	gen FE = I*1000 + H
	xtset FE Year
	
	xtreg Y T_h* ZLB_h* ///
	X*_h*,	///
	fe vce(cluster I)
	
	estimates store ZLB_stackedOLS

  
    forvalues h = 0/4 {
        erase "ZLB_stack`h'.dta"
    }

}

esttab ZLB_stackedOLS , b(2) se keep(*T_*)

// average multiplier non ZLB
lincom (T_h0+T_h1+T_h2+T_h3+T_h4)/5    


// ZLB  multiplier by horizon 
forvalues h = 0/4 {
    lincom T_h`h' + ZLB_h`h', 
}

// average multiplier
lincom ((T_h0+ZLB_h0)+(T_h1+ZLB_h1)+(T_h2+ZLB_h2)+(T_h3+ZLB_h3)+(T_h4+ZLB_h4))/5 
	
	
	
	
// Inflation expectations and ZLB 		
drop _x10
	
eststo clear
forvalues h = 0/5 {
	xtreg  S`h'_real_GDP                             ///
           c.dCAPB##i.ZLB##c.inflation_pct _x*,              ///
           fe vce(cluster country_id)
    estimates store IE_ZLB_OLS`h'
}	
	
	
	
quietly {
	forvalues h = 0/4 {
		preserve

		gen Y = S`h'_real_GDP
		gen I = country_id
		
		gen T_h`h' = dCAPB          
		
		gen ZLB_h`h' = dCAPB * ZLB      // ZLB interaction effect
		gen ZLB_dummy = ZLB  
		
		gen inflation = inflation_pct
		
		gen IE_h`h'      = dCAPB * inflation_pct	// inflation effect on CAPB
		gen ZLB_IE_h`h'  = dCAPB * ZLB * inflation_pct 	// inflation effect on multiplier at ZLB 
		

		
		forvalues k = 1/9 {
			rename _x`k' X`k'_h`h'
        }
		
		gen H = `h'
		// keep only the variables we need
		keep Y T_h`h' ZLB_h`h' X*_h`h' I H Year ZLB_dummy IE_h`h' ZLB_IE_h`h' inflation
		save "IE_ZLB_stack`h'.dta", replace
		restore
	}
  
	// Appending stacked files
    preserve
    clear
    forvalues h = 0/4 {
        append using "IE_ZLB_stack`h'.dta"
    }

    // Replacing missing values with 0 
    ds T_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
	
	ds ZLB_h*
	foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
	
	ds IE_h*
	foreach var of varlist `r(varlist)' {
		replace `var' = 0 if missing(`var')
	}
	
	ds ZLB_IE_h*
	foreach var of varlist `r(varlist)' {
		replace `var' = 0 if missing(`var')
	}
	
    ds X*_h*
    foreach var of varlist `r(varlist)' {
        replace `var' = 0 if missing(`var')
    }
    save "IE_ZLB_stack.dta", replace
	

	use IE_ZLB_stack.dta, clear
	gen FE = I*1000 + H
	xtset FE Year
	
	xtreg Y T_h* ZLB_h* IE_h*  ZLB_IE_h* ///
	X*_h*,	///
	fe vce(cluster I)
	
	estimates store IE_ZLB_stackedOLS

  
    forvalues h = 0/4 {
        erase "IE_ZLB_stack`h'.dta"
    }
}	
	

// CAPB, fiscal multiplier in normal times
esttab IE_ZLB_stackedOLS , b(2) se keep(*T_*) 
lincom (T_h0+T_h1+T_h2+T_h3+T_h4)/5   


// ZLB interaction effect, fiscal multiplier when IE = 0
esttab IE_ZLB_stackedOLS , b(2) se keep(ZLB_h*)	
lincom (ZLB_h0+ZLB_h1+ZLB_h2+ZLB_h3+ZLB_h4)/5   


// inflation expectations effect in non-ZLB
esttab IE_ZLB_stackedOLS , b(2) se keep(IE_h*)
lincom (IE_h0+IE_h1+IE_h2+IE_h3+IE_h4)/5   


// inflation expectations effect on multiplier at ZLB 
esttab IE_ZLB_stackedOLS , b(2) se keep(ZLB_IE_h*)
lincom (ZLB_IE_h0+ZLB_IE_h1+ZLB_IE_h2+ZLB_IE_h3+ZLB_IE_h4)/5  	
	
	
	
	
	
// plotting the IRF for Table 1
clear
estimates restore stackedOLS
set obs 5
gen h = _n
gen irf = .
gen se = .

forvalues h = 0/4 {
    local coefname = "T_h`h'"
    replace irf = _b[`coefname'] in `=`h'+1'
    replace se  = _se[`coefname'] in `=`h'+1'
}

gen ub = irf + 1.96 * se    // upper 95% CI
gen lb = irf - 1.96 * se    // lower 95% CI


twoway (rarea ub lb h, color(gs12)) ///
       (line irf h, lwidth(medthick) lcolor(blue)), ///
       yline(0, lpattern(dash)) ///
       xtitle("Horizon (Years)") ///
       ytitle("Percent change in GDP") ///
       title("IRF to 1pp CAPB Shock Pre-COVID-19") ///
       legend(off)
