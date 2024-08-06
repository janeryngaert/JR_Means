//concatenate.do
//This file concatenates the Stata datasets created by importing the individual SCE spreadsheet files into a single data set, saved as raw_all.dta 
//Merge the data.
use raw_13_16
append using raw_17_19
append using raw_latest

//Create a monthly Stata date variable.
gen date1 = ym(year, month)
format date1 %tm

//Annotate the data.
label data "Survey of Consumer Expectations: All Raw Data"
notes: raw_all.dta: Created from concatenate.do running Stata `c(edition_real)' Version `c(stata_version)' by `c(username)' on `c(hostname)', a `c(machine_type)' running `c(os)' version `c(osdtl)'.
notes: raw_all.dta: Created on $S_DATE at $S_TIME.

//The data set is saved as all.dta. Its directory address should indicate to users that it is raw SCE data.
save raw_all, replace

