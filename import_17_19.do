//This file imports the SCE public microdata from 2017-2019 and saves them to raw_17_19.dta
clear all

//Note that the first row of this spreadsheet contains bibliographic information about the SCE. Therefore, we specify the range of data to be imported explicitly so that the ``firstrow'' option points to the cells containing variable names.
import excel ".\sce_raw_data\FRBNY-SCE-Public-Microdata-Complete-17-19.xlsx", sheet("Data") cellrange(A2:HL47683) firstrow case(lower)

//The date variable is imported as a long integer. turn it into a string so that we can more easily extract the year and month.
tostring date, replace force

//The year is the first four characters. The month is the last two characters.
gen year = substr(date,1,4)
gen month = substr(date,5,2)

//Convert the year and month into integers.
destring year month, replace force

//Annotate the data.
label data "Survey of Consumer Expectations: 2017-2019 Raw Data"
notes: Created from import_17_19.do running Stata `c(edition_real)' Version `c(stata_version)' by `c(username)' on `c(hostname)', a `c(machine_type)' running `c(os)' version `c(osdtl)'.
notes: Created on $S_DATE at $S_TIME.

save raw_17_19, replace
