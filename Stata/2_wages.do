* program name: "2_worker_flows.do"
* Project: Heterogeneity and unemployment flows across countries
* Jonathan Créchet (2022)
* Created: 2020
* Last updated: March 2023

* Compute wage statistics using ORG sample files

***
***

capture log close
log using Log\2_wages, replace text

* data
use "$working_dir\Dta\CPS_data", clear

* variables
keep year month cpsidp mish earnwt empstat age hourwage paidhour earnweek uhrswork1 ahrswork1 classwkr 

* ORG samples
keep if mish == 4 | mish == 8

* Merge with CPI99
merge m:1 year using "$working_dir\Dta\cpi99"
rename cpi cpi99

* age restrictions
keep if age >= 20 & age <= 60

** LF status restriction
 * employed
keep if empstat >= 10 & empstat <= 12

* exclude self-employed 
keep if classwkr >= 22 & classwkr <= 28


**
**

* weekly earnings descriptives
sum earnweek, d

**
**

** identify top codes 
* (see https://cps.ipums.org/cps-action/variables/EARNWEEK#codes_section)
gen topcode = 0

* 1989-1997
replace topcode = 1 if earnweek >= 1923 & year >= 1989 & year <= 1997

* 1998-onward
replace topcode = 1 if earnweek >= 2884.61 & year >= 1998 & year <= 2023


**
**

* nomimal hourly earnings: worker paid by hourly wage
gen nominal_wage = hourwage if paidhour == 2

* salary workers
replace nominal_wage = earnweek / uhrswork1 if paidhour == 1 & uhrswork1 >= 1 & uhrswork1 <= 99

* workers declaring variable hours: use actual hours
replace nominal_wage = earnweek / ahrswork1 if paidhour == 1 & uhrswork1 == 997 & ahrswork1 >= 1 & ahrswork1 <= 99 

* count missing
count if nominal_wage == .

* distribution
sum nominal_wage, d

* real wage
gen real_wage = nominal_wage * cpi99
sum real_wage, d

**
**

* top coded: times 1.4. (Lemieux, 2006)
replace real_wage = real_wage * 1.4 if topcode == 1

* "trim" bottom and top 1% (by year) to eliminate outliers and top coded
by year, sort: egen bottom = pctile(real_wage), p(1)
by year: tab bottom

by year, sort: egen top = pctile(real_wage), p(99)
by year: tab top

* check
gen n = (real_wage >= bottom & real_wage <= top)
tab topcode if n
keep if n

drop n bottom top

**
**

* normalize weights in remaining sample
by year month, sort: egen total_earnwt = sum( earnwt )
gen normalized_weights = earnwt / total_earnwt
sum normalized_weights, d

* log wage distribution
gen log_wage = log( real_wage )
sum log_wage [w = normalized_weights], d

* close log file
log close

**
**

* Reference

* Lemieux, Thomas. "Increasing residual wage inequality: Composition effects, noisy data, or rising demand for skill?." 
* American economic review 96, no. 3 (2006): 461-498.
