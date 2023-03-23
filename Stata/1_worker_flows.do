* program name: "2_worker_flows.do"
* Project: Heterogeneity and unemployment flows across countries
* Jonathan Créchet (2022)
* Created: 2020
* Last updated: March 2023

* Compute aggregate worker-flow statistics and transition rates by age, u duration, job tenure

*****
*****

capture log close
log using Log\1_worker_flows, replace text

* data
use "$working_dir\Dta\CPS_data", clear

* variables
drop earnwt earnweek hourwage paidhour uhrswork1 ahrswork1 classwkr 

***
***

* 1. PANEL PREPARATION

* date
gen year_str = string(year)
gen month_str = string(month)
gen date_str = year_str + "m" + month_str
gen date = monthly( date_str, "YM" )
drop *_str

* declare panel
xtset cpsidp date, monthly
order cpsidp date
sort cpsidp date


* sample of individuals in itw month two, three, four, six, seven, eight with available 
* info in previous month (and that can be longitudinally matched based on cpsidp)
gen n = .

foreach i of numlist 2 3 4 6 7 8 {
replace n = 0 if mish==`i'
replace n = 1 if mish==`i' & L.empstat!=.
}

* show number of longitudinally matched id by month in sample
foreach i of numlist 2 3 4 6 7 8 {

disp " "
disp "month in sample:" `i'
tab n if mish == `i', miss

}


** check id mismatch

* dummies for race
decode race, gen(race_string)
gen white = ( strpos(race_,"white")>0 )
gen black = ( strpos(race_,"black")>0 )
gen asian = ( strpos(race_,"asian")>0 )
gen pacific = ( strpos(race_,"pacific")>0 )
gen native  = ( strpos(race_,"hawaiian") + strpos(race_,"american") + strpos(race_,"aleut") + strpos(race_,"eskimo") > 0 )
gen tmp = white + black + asian + native + pacific
gen no_race = tmp == 0

* check
tab race if white
tab race if black
tab race if asian
tab race if native
tab race if pacific
tab race if no_race


* identify mismatched based on sex
gen mismatch = .
foreach i of numlist 2 3 4 6 7 8 {
replace mismatch = 0 if mish == `i' & n == 1
replace mismatch = 1 if mish == `i' & n == 1 & ( L.sex != sex )
}

* age
foreach i of numlist 2 3 4 6 7 8 {
replace mismatch = 1 if mish == `i' & n == 1 & ( L.age > age | age > L.age+1 )
}

* race
foreach i of numlist 2 3 4 6 7 8 {
replace mismatch = 1 if mish == `i' & n == 1 & ( L.race != race ) & ( L.no_race == 0 | no_race == 0 )
}

* count number of mismatched per id
by cpsidp, sort: egen sum_mismatch = sum(mismatch)

* sample: remove mismatched and id with three or more mismatched in total
foreach i of numlist 2 3 4 6 7 8 {
replace n = 0 if mish == `i' & ( mismatch == 1 | sum_mismatch > 2 )
}

* final sample
tab n, miss


* tabulate by interview month
foreach i of numlist 2 3 4 6 7 8 {

disp " "
disp "month in sample:" `i'
tab mismatch if mish == `i', miss
tab sum_mism if mish == `i'

}

* cleanup
drop sex *race* mismatch sum_mism white black asian pacific native tmp 

**
**

** panel weight
gen panelweights = (1/2)*( L.wtfinl + wtfinl ) if n==1
sum panelweights


******
******

* 2. LABOR MARKET STATUS VARIABLES

* Employment 
gen E = (empstat >= 10 & empstat <= 12)

* Unemployment
gen U = (empstat >= 20 & empstat <= 22)

* Non-participant
gen I = (empstat >= 30 & empstat <= 36)

* Non-employed
gen N = U+I 			//non-employment state


* check
tab empstat if E, miss
tab empstat if U, miss
tab empstat if I, miss
tab empstat if N, miss

* save intermediate datasets (with panel variables and LM dummies)
save "$working_dir\Dta\CPS_data_2", replace

******
******

* AGGREGATE WORKER FLOWS 
do 1a_aggregate

******
******

* TRANSITION RATES BY AGE
do 1b_age

******
******

* TRANSITION RATES  TENURE
do 1c_job_tenure.do


******
******

* TRANSITION RATES  U DURATION
do 1d_unemployment_duration.do

******
******

* erase dta file used in worker flows analysis
erase "$working_dir\Dta\CPS_data_2.dta"

******
******

log close
