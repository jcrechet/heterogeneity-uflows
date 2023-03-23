* program name: "cross_country.do"
* Project: Heterogeneity and unemployment flows across countries
* Jonathan Créchet (2022)
* Created: 2020
* Last updated: March 2023

* Compute cross-country statistics for producing figure 3 and statistics used in calibration

**
**

clear

capture log close
log using Log/3_cross-country, text replace

**
**

** Construct dataset using OECD raw data

**
**

* LFS data (unemployment)
import delimited Raw_data/LFS.csv, clear

* rename variable names (for consistency with other data files)
rename ï ctry
rename time year

* keep population and variables of interest
keep if age == 900000 // age: 15 to 64
keep if sex == "MW" // men and women
keep ctry country year series value //variables of interest
keep if series == "U" | series == "E"

* reshape data (wide)
reshape wide value, i(ctry country year) j(series) string

* rename variables of interest
rename valueE E 
rename valueU U

* save temporary dta file
save tmp, replace

**
**

** Unemployment by duration
import delimited Raw_data/unemployment_by_duration.csv, clear

* renaming key data
rename ï ctry
rename time year

* population
keep if sex == "MW" // men and women
keep if age == 900000 // all ages

* variables
keep ctry country year duration value //variables of interest

* duration less than one month
keep if duration == "UN1"

* cleanup
drop duration
rename value u_short_term

* merging data
merge 1:1 ctry year using tmp
drop _merge
save tmp, replace

**
**

** Elsby, Hobijn, Sahin (2013) flows

* country codes
use tmp, clear
keep country ctry
duplicates drop 
save Dta/ctry_codes.dta, replace

* import excel file
import excel Raw_data\UnemploymentDynamicsInTheOECD.xlsm, sheet("Calculations") cellrange(A3:IV717) firstrow clear

* variables
keep country year u f s

* country codes
merge m:m country using Dta/ctry_codes.dta
drop _merge

* merge data with master file
drop if year == .
merge 1:1 ctry year using tmp
drop _merge
order country ctry year

* Labels
lab var f "U outflows -- Elsby et al. (2013)"
lab var s "U inflows -- Elsby et al. (2013)"

* Monthly transition probabilities
gen ue_ehs = 1-exp(-f)
gen eu_ehs = 1-exp(-s)

save tmp, replace

**
**

* Taxes
import delimited Raw_data\tax_wages.csv, clear
rename cou ctry

* employee and employer's social security
keep if ïindicator == "2_2" | ïindicator == "2_3"

* single, 100% of average earnings, no child
keep if fam_type == "SINGLE2" 

* cleanup
keep ctry country year ïindicator value  
rename ïindicator series
reshape wide value, i(ctry country year) j(series) string

* cleanup
rename value2_2 ssc_employee
rename value2_3 ssc_employer

* merge data
drop if year == .
merge 1:1 ctry year using tmp
drop _merge

* order variables and sort by country, year
order country ctry year
save tmp, replace

**
**

* Unemployment benefits

* UB (net) replacement ratio
import delimited Raw_data\Net_replacement_rate_ub.csv, clear
rename ï ctry
keep ctry country year family duration earnings value includesocialassistancebenefits satopups hbtopups includehousingbenefits

* single
keep if family == "SINGLE"

* Average earnings
keep if earnings == "AW"

* duration less than six months
keep if duration >=0 & duration <= 6

* no social assistance/housing benef
keep if includesocialassistancebenefits == "No"
keep if includehousingbenefits == "No"

* satopups, hbtopups??
keep if satopups == 0
keep if hbtopups == 0

* cleanup
drop family earnings include* sato hbto
order ctry country year duration value
sort ctry country year duration value
rename value uemployment_benefits

* one year average over fduration
collapse uemployment_benefits, by(ctry country year)

* merge data with master file
drop if year == .
merge 1:1 ctry year using tmp
drop _merge

* order variables and sort by country, year
order country ctry year
save tmp, replace

**
**

* Employment Protection Legislation (EPL) index

* regular, individual
import delimited Raw_data\EPL_individual_regular.csv, clear

* EPL v1
keep if series == "EPR_V1"

* cleanup
rename ï ctry
rename time year
keep ctry country year value  
rename value EPR_V1

* merge with master file
merge 1:1 ctry year using tmp
drop _merge
save tmp, replace


***
***

* Figure 3: EPL index vs Elsby et al. (2013) probabilities
preserve

keep if ue_ehs!=. & eu_ehs!=. & EPR_V1 != . 
keep if year >= 1985 & year <= 2009

* collapse
collapse ue_ehs eu_ehs EPR_V1, by(ctry)

* transition rates in log deviation from the US
egen ue_US = max(ue_ehs)
egen eu_US = max(eu_ehs)

gen ue_rel = log(ue_ehs) - log(ue_US)
gen eu_rel = log(eu_ehs) - log(eu_US)

* Scatter plot: unemployment outflows vs EPL
scatter ue_rel EPR_V1, mlabel(ctry) || lfit ue_rel EPR_V1 || ///
lfit ue_rel EPR_V1 if ctry != "PRT" & ctry != "IRL" & ctry != "NOR" & ctry != "SWE", ///
ytitle("UE rate (rel. to the U.S.)") xtitle("EPL index") legend(off) graphregion(color(white))
graph export "$results_path\Figures\Figure_3a.eps", as(eps) replace
graph export "$results_path\Figures\Figure_3a.png", as(png) replace

* Scatter plot: unemployment inflows vs EPL
scatter eu_rel EPR_V1, mlabel(ctry) || lfit eu_rel EPR_V1 || ///
lfit eu_rel EPR_V1 if ctry != "PRT" & ctry != "IRL" & ctry != "NOR" & ctry != "SWE", ///
ytitle("EU rate (rel. to the U.S.)") xtitle("EPL index") legend(off) graphregion(color(white))
graph export "$results_path\Figures\Figure_3b.eps", as(eps) replace
graph export "$results_path\Figures\Figure_3b.png", as(png) replace

restore

***
***

** Data for calibration

* UE and EU rates based on Rogerson, Shimer (2011)
gen ue = u_short_term/U
gen eu = 1-(1-ue)^(U/E)

* sample
keep if year >= 1995 & year <= 2018
keep if ctry == "DEU" | ctry == "ESP" | ctry == "FRA" | ctry == "ITA" | ctry == "PRT" | ctry == "USA"

* mean 
collapse ue eu uemployment_benefits ssc_employee ssc_employer, by( country ctry )

* export to Matlab
export delimited using "$Matlab_directory\Data\cross_country_data.csv", replace

* erase temporary datafile
erase tmp.dta

**

** References

* Elsby, Michael WL, Bart Hobijn, and Ayşegül Şahin. 
* "Unemployment Dynamics in the OECD." Review of Economics and Statistics 95, no. 2 (2013): 530-548.

* Rogerson, Richard, and Robert Shimer. 
* "Search in macroeconomic models of the labor market." In Handbook of labor economics, vol. 4, pp. 619-700. Elsevier, 2011.






