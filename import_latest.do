//This file imports the SCE public microdata from 2020 to the end of the sample period and saves them to raw_latest.dta
clear all

//Note that the first row of this spreadsheet contains bibliographic information about the SCE. Therefore, we specify the range of data to be imported explicitly so that the ``firstrow'' option points to the cells containing variable names.
import excel ".\sce_raw_data\frbny-sce-public-microdata-latest.xlsx", sheet("Data") cellrange(A2:HL30155) firstrow case(lower)

//The date variable is imported as a long integer. turn it into a string so that we can more easily extract the year and month.
tostring date, replace force

//The year is the first four characters. The month is the last two characters.
gen year = substr(date,1,4)
gen month = substr(date,5,2)

//Convert the year and month into integers.
destring year month, replace force

label data "Survey of Consumer Expectations: 2020-latest Raw Data"
notes: Created from import_latest.do running Stata `c(edition_real)' Version `c(stata_version)' by `c(username)' on `c(hostname)', a `c(machine_type)' running `c(os)' version `c(osdtl)'.
notes: Created on $S_DATE at $S_TIME.
save raw_latest, replace
