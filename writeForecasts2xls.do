//writeForecasts2xls.do
//This file extracts the point forecasts and bin-probability forecasts and writes them to .xls files for use by the Matlab files that implement JMR's distributional fitting procedure.


//Load the data
use raw_all

//Declare the data to be a panel with userid indicating households and date1 indicating time.
tsset userid date1

//Generate a variable containing the number of records in the data set.
gen n = _n

//JANE, IS THERE A CODEBOOK WHICH DEFINES THE VARIABLES? THE QUESTIONAIRE DOES NOT LIST EVERY VARIABLE.
// JR 5/28 - to my knowledge and search of the website there is not

//Drop variables generated by the NY Fed's distribution-fitting routines, so that they do not appear in wildcard-generated question lists.
drop q9_mean q9_cent25 q9_cent50 q9_cent75 q9_iqr q9_probdeflation q9_var q9c_mean q9c_cent25 q9c_cent50 q9c_cent75 q9c_iqr q9c_probdeflation q9c_var

//Calculate the number of bins containing positive probability for each response.
foreach question in q9 q9c {

	//Create a macro containing the names of the bin-probability question.
	local binProbabilityQuestions 
	forvalues x = 1/10{
		local binProbabilityQuestions = "`binProbabilityQuestions' `question'_bin`x'"
	}
	disp "`binProbabilityQuestions'"

	//How many non-missing answers to the bin questions are there? 
	egen _rm=rownonmiss(`question'*)

	//Cycle through the reported probabilities.
	foreach var in `binProbabilityQuestions'
					{	
						//Replace the data with a zero if it is missing and there is at least one non-missing probability
						replace `var'=0 if `var'==. & _rm>0 & _rm~=. 
						gen f_`var'=(`var'==0) //Create a variable to indicate whether this particular variable equals zero.
					}
	drop _rm //Drop the count of non-missing probability variables.

	//Generate the count of non-zero probability variables. 
	//(This now equates missing values with zeros.)
	egen _rm_`question'=rowtotal(f_`question'*) 
	
}

//Write the xls files for one-year forecasts.
//Keep only the data we need for the excel files. 
keep q9_bin1 q9_bin2 q9_bin3 q9_bin4 q9_bin5 q9_bin6 q9_bin7 q9_bin8 q9_bin9 q9_bin10 q8v2part2 _rm_q9

//These files get used to fit an implied density with the point forecasts and the bin probabilities. If two or more records share the same data, then we only need to calculate the implied density once. Therefore, we remove duplicates by ``collapsing'' the data by the relevant variables.
//Somewhat surprisingly, this operation reduces the number of observations by a factor of at least 2. 
collapse (mean) q8v2part2, by(q9_bin1 q9_bin2 q9_bin3 q9_bin4 q9_bin5 q9_bin6 q9_bin7 q9_bin8 q9_bin9 q9_bin10 q8v2part2 _rm_q9)

//Drop any records with no bin frequency data.
drop if _rm_q9 == 10

//Drop any records without point estimates
drop if mi(q8v2part2) 

//Drop any records with bin frequency data that do not sum to 100 percent. (The CAR system should ensure that this never happens, but... .)
gen sum_q9 = q9_bin1 + q9_bin2 + q9_bin3 + q9_bin4 + q9_bin5 + q9_bin6 + q9_bin7 + q9_bin8 + q9_bin9 + q9_bin10
drop if sum_q9 != 100

//Give the variables a fixed order for the xls files. 
order q8v2part2 q9_bin1 q9_bin2 q9_bin3 q9_bin4 q9_bin5 q9_bin6 q9_bin7 q9_bin8 q9_bin9 q9_bin10 _rm_q9

//Next, we generate three excel spreadsheets.  
preserve

//The first spreadsheet contains only those records with all probability placed into one bin.
keep if _rm_q9 == 9
export excel using "1year_1bin.xls", replace

restore
preserve

//The second spreadsheet contains only those records with all probabilty placed into two bins.
keep if _rm_q9 == 8
export excel using "C:\Users\Jane Ryngaert\Documents\MATLAB\1year_2bin.xls", replace

restore
preserve

//The third spreadsheet contains only those records with positive probability placed into three or more bins.
keep if _rm_q9 < 8

*export excel using "1year_3bins.xlsx", replace
restore

