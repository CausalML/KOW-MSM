



*cd "~/Desktop/case study/data/stata"
cd "~/research/Research/projects/Cornell/Longitudinal/Rcode/Case Study/MACS/data/stata/"
cd "/Volumes/biostatistik$/STAFF/Michele Santacatterina/Research/projects/Cornell/Longitudinal/Rcode/Case Study/MACS/data/stata/"

global PATH2 = "`c(pwd)'"



*lab results
use "$PATH2/lab_rslt.dta", clear

bysort caseid ldaty (visit): gen order = _n
order order, after(ldaty)

bysort caseid ldaty: gen totalv = _N
order totalv, after(order)

*ANDREA

gen month = .
order month, after(totalv)
replace month = 6 if totalv == 1
replace month = 1 if totalv == 2 & order == 1
replace month = 6 if totalv == 2 & order == 2
replace month = 1 if totalv == 3 & order == 1
replace month = 4 if totalv == 3 & order == 2
replace month = 8 if totalv == 3 & order == 3

gen datenew = mdy(month ,1 , ldaty)
format datenew %td
order datenew, after(ldaty)

bysort caseid (datenew): gen timestart = 0 if _n ==1
order timestart, after(datenew)
bysort caseid (datenew): replace timestart = datenew[_n]-datenew[1] if _n > 1

bysort caseid (datenew): gen timeend = timestart[_n+1] if _n != _N
order timeend, after(timestart)
bysort caseid (datenew): replace timeend = mdy(12, 31, year(datenew)) - datenew[1] if _n == _N


keep caseid visit leu3n wbc rbc plate ldaty timestart timeend
reshape wide leu3n wbc rbc plate ldaty timestart timeend, i(caseid) j(visit)

saveold "$PATH2/analysis/lab_rslt.dta", replace


*treatments
use "$PATH2/drugf1.dta", clear

sort caseid visit

quietly by caseid visit:  gen dup = cond(_N==1,0,_n)
drop if dup>1 & dup!=.

keep caseid visit avqy drgav dup
drop dup

reshape wide avqy drgav, i(caseid) j(visit)

saveold "$PATH2/analysis/drugf1.dta", replace


*symptoms
use "$PATH2/section4.dta", clear

keep caseid visit fev2w trs2w dia2w wtlos oleuk herpz 

gen symp = 0 
replace symp = 1 if fev2w==2 | trs2w==2 | dia2w==2 | wtlos==2 | oleuk==2 | herpz==2 

keep caseid visit symp

reshape wide symp, i(caseid) j(visit)

saveold "$PATH2/analysis/section4.dta", replace


*info
use "$PATH2/section2.dta", clear

sort caseid visit

keep caseid visit borny

reshape wide borny, i(caseid) j(visit)

saveold "$PATH2/analysis/section2.dta", replace


*outcome
/*
1= AIDS-prior dx
2= AIDS-no prior dx
3= Not AIDS
4= Unknown
Blank= Missing*/

use "$PATH2/outcome.dta", clear

keep caseid death dthdateyy

sort caseid

saveold "$PATH2/analysis/outcome.dta", replace



use "$PATH2/analysis/outcome.dta", clear

merge 1:1 caseid using "$PATH2/analysis/lab_rslt.dta", gen(_mergelab)
merge 1:1 caseid using "$PATH2/analysis/drugf1.dta", gen(_mergedrug)
merge 1:1 caseid using "$PATH2/analysis/section4.dta", gen(_mergesec4)
merge 1:1 caseid using "$PATH2/analysis/section2.dta", gen(_mergesec2)

drop _merge*


reshape long ldaty wbc rbc plate leu3n avqy drgav symp borny timestart timeend, i(caseid) j(visit)

saveold "$PATH2/analysis/complete_temp.dta", replace




use "$PATH2/analysis/complete_temp.dta", clear

gen treat=1
replace treat=0 if drgav==.

keep if ldaty!=. & wbc!=. & rbc!=. &  plate!=. &  leu3n!=. &  treat!=. &  symp!=. &  borny!=. 

rename death death_temp
replace death_temp = 0 if death_temp==.

gen death = 0
bysort caseid: replace death = 1 if _n == _N & death_temp>=1

gen censor = 0
bysort caseid: replace censor = 1 if _n == _N & death==0


by caseid: egen firstvisityy = min(ldaty)
by caseid: egen lastvisityy = max(ldaty)

gen lastyy = lastvisityy if lastvisityy<=dthdateyy

gen time = lastyy - firstvisityy
gen age = firstvisityy - borny

drop if age <0

by caseid: gen tempv = _n


by caseid, sort: gen treat_ipw = treat[_n-1] if _n>1
by caseid, sort: gen censor_ipw = censor[_n-1] if _n>1

sort caseid visit
replace treat_ipw = treat if treat_ipw == .
replace censor_ipw = censor if censor_ipw == .

sort caseid tempv

by caseid, sort: gen tempt0 = _n if treat==0
by caseid, sort: gen tempt1 = _n if treat==1

by caseid, sort: gen tempCD4 = leu3n if treat==1
by caseid, sort: egen CD4start = min(tempCD4) 


by caseid, sort: gen sumt = sum(treat)

*drop _merge*
sort tempv caseid 

*keep if sumt >0

keep if firstvisityy >=2001

*keep if CD4start >500 
*keep if CD4start >500 & CD4start<=750
*keep if CD4start >200 & CD4start<=350
*keep if CD4start >350 & CD4start<=500



saveold "$PATH2/analysis/complete_2001_180418.dta", replace

*saveold "$PATH2/analysis/complete2001_350500.dta", replace














