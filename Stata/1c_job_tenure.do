* TRANSITION PROBABILITIES BY JOB TENURE

* data
use "$working_dir\Dta\CPS_data_2", clear

* weights
sort cpsidp date
gen ten_weights = L.jtsuppwt

* lagged tenure in years
gen tenure0 = floor(L.jtyears) if L.jtyears < 99 & n == 1

* Employment dummy, lagged
gen E0 = L.E if L.jtresp == 2 & L.jtyears < 99 & n == 1

* Employment dummy, lagged, with interview dependant info 
gen E0b = L.E if L.jtresp == 2 & L.jtyears < 99 & empsame >= 1 & empsame <= 2 & n == 1

* Employment separation flow dummies
gen E_U = ( L.E == 1 & U == 1 ) if L.jtresp == 2 & L.jtyears < 99 & n == 1
gen E_N = ( L.E == 1 & ( U == 1 | I == 1 ) ) if L.jtresp == 2 & L.jtyears < 99 & n == 1

* EE dummy
gen E_E = ( L.E == 1 & E == 1 & empsame == 1 ) if L.jtresp == 2 & L.jtyears < 99 & empsame >= 1 & empsame <= 2 & n == 1

* collapse by tenure
collapse (sum) E0 (sum) E0b (sum) E_U (sum) E_N (sum) E_E [w=ten_weights], by( date tenure )

* probabilities
gen eu = E_U / E0
gen en = E_N / E0
gen ee = E_E / E0b

* cleanup
keep date tenure eu en ee
rename tenure0 tenure

**
**

* export to R folder
export delimited "$R_directory\Data\CPS_flows_tenure", replace
