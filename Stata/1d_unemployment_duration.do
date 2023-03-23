* TRANSITION PROBABILITIES BY U DURATION

* data
use "$working_dir\Dta\CPS_data_2", clear



* tabulate distribution
tab durunemp if U == 1, miss

* gen lag unemployment duration in months
gen duration0 = floor( L.durunemp/(52/12) ) if L.U == 1 & L.durunem < 999
tab duration0 if L.U == 1, miss

* dummy for lagged state
gen U0 = L.U if n == 1 & duration0 != .

* dummy for transition using longitudinal cpsidp linkage and the number of weeks in unemployment
gen U_E = ( L.U == 1 & ( E == 1 | ( U == 1 & durunemp < 4 ) ) ) if n == 1 & duration0 != .

* collapse by unemployment duration
collapse (sum) U0 (sum) U_E [w=panelweight], by( date year month duration )

* probabilities
gen ue = U_E / U0

* cleanup
rename duration0 duration
keep date year month duration ue
keep if duration != .
keep if duration < 25
order date month duration



** RMA seasonal adjustment
xtset duration date

* 13-month centered moving average
tssmooth ma moving_average = ue, window(6 1 6) 

* actual-to-smoothed ratio
gen ma_ratio = ue / moving_average

* compute calendar-month average (by u-duration)
egen month_average = mean( ma_ratio ), by( month duration )

* center the ratio index around one (normalize)
egen index_average = mean( month_average ), by( duration )
gen month_average_centered = month_average + ( 1 - index_average )

* compute the seasonnaly adjusted series
gen ue_sa = ue / month_average_centered

* cleanup
drop moving_average ma_ratio month_average index_average month_average_centered

drop ue
rename ue_sa ue


********
********


* Quarterly averages

* quarter
gen quarter = ceil(month/3)

* time variables
gen year_str = string(year)
gen quarter_str = string(quarter)
gen date_str = year_str + "q" + quarter_str
gen date_ = quarterly( date_str, "YQ" )
drop date *_str
rename date_ date

* collapse to get unweighted quarterly averages by duration groups
collapse ue, by( date year quarter duration )
order date year quarter duration
sort date year quarter duration

xtset duration date, quarterly

* export to R folder
export delimited "$R_directory\Data\CPS_flows_uduration", replace
