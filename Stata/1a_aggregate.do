* AGGREGATE TRANSITION RATES

* data
use "$working_dir\Dta\CPS_data_2", clear

** DUMMIES

** from E
gen E0 = L.E if n == 1
gen E_U = ( L.E == 1 & U == 1 ) if n == 1
gen E_I = ( L.E == 1 & I == 1 ) if n == 1

** from U
gen U0 = L.U if n == 1
gen U_E = ( L.U == 1 & E == 1 ) if n == 1
gen U_I = ( L.U == 1 & I == 1 ) if n == 1

** from I
gen I0 = L.I if n == 1
gen I_E = ( L.I == 1 & E == 1 ) if n == 1
gen I_U = ( L.I == 1 & U == 1 ) if n == 1


** EE (employer to employer) dummies
tab empsame if E0 == 1 & E == 1
by year, sort: tab empsame if E0 == 1 & E == 1
sort cpsidp date

* gen dummies for employment stayers with non missing samemp
gen E0b = L.E if n == 1 & empsame >= 1 & empsame <= 2 
gen E_E  = (L.E == 1 & E == 1 & empsame == 1) if n == 1 & empsame >= 1 & empsame <= 2


* list of variables
#delimit ; 

global varlist =  " (sum) E (sum) U (sum) I
				    (sum) E0 (sum) U0 (sum) I0
				    (sum) E_U (sum) E_I 
				    (sum) U_E (sum) U_I
				    (sum) I_E (sum) I_U 
					(sum) E0b (sum) E_E ";

#delimit cr

* Collapse to get aggregates
collapse $varlist [w=panelweights], by( date year month )
sort date

* compute ratio to pop
gen e = E / (E + U + I)
gen u = U / (E + U + I)
gen i = 1 - e - u

* compute transition probabilities between E, U, I
gen eu = E_U / E0
gen ei = E_I / E0

gen ue = U_E / U0
gen ui = U_I / U0

gen ie = I_E / I0
gen iu = I_U / I0

* compute EE transitions
gen ee = E_E / E0b

* cleanup
keep date year month e u i eu ei ue ui ie iu ee


*******
*******


** RMA seasonal adjustment
tsset date

* loop over transition probabilities
foreach var in e u i eu ei ue ui ie iu ee {

* 13-month centered moving average
tssmooth ma moving_average = `var', window(6 1 6) 

* actual-to-smoothed ratio
gen ma_ratio = `var' / moving_average

* compute calendar-month average
egen month_average = mean( ma_ratio ), by( month )

* center the ratio index around one
egen index_average = mean( month_average )
gen month_average_centered = 1 + ( month_average - index_average )

* compute the seasonnaly adjusted series
gen `var'_sa = `var' / month_average_centered

* cleanup
drop moving_average ma_ratio month_average index_average month_average_centered

}

* rename sa variables and drop raw series
drop e u i eu ei ue ui ie iu ee
rename *_sa *


*******
*******


** quarterly averages

* quarter
gen quarter = ceil(month/3)

* time variables
gen year_str = string(year)
gen quarter_str = string(quarter)
gen date_str = year_str + "q" + quarter_str
gen date_ = quarterly( date_str, "YQ" )
drop date *_str
rename date_ date

* collapse to get unweighted quarterly averages
collapse e u i eu ei ue ui ie iu ee, by( date year quarter )
tsset date, quarter
order date year quarter
sort date year quarter


*******
*******

* export in R folder
export delimited "$R_directory\Data\CPS_flows", replace
