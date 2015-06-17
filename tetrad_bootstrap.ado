*! v1.5, CTA program, S Bauldry, 01jun2015

capture program drop tetrad_bootstrap
program define tetrad_bootstrap, rclass
	version 13
	syntax varlist(min = 4 numeric) [if] [in], ///
		   icm1(name)         /// implied covariance matrix 1
		   [seed(string)      /// random number seed
		   reps(integer 100)] /// number of bootstrap replications

	marksample touse
	
	*** run tetrad on original data to obtain baseline test statistic
	tetrad `varlist' if `touse', icm1(`icm1') boot(1)
	local T = r(T)
	local df = r(df)
	
	*** save data before transforming
	preserve
	tempfile pres
	qui save `pres'
	
	qui keep `varlist' __000000
	qui keep if `touse'
	
	recast double _all 
	
	*** put observed variables from ICs in same order as SCM
	tempname ICM1
	local k : list sizeof varlist
	mat `ICM1' = J(`k', `k', .)
	foreach row of local varlist {
		foreach col of local varlist {
			local i : list posof "`row'" in varlist
			local j : list posof "`col'" in varlist
			local n = rownumb(`icm1', "`row'")
			local m = colnumb(`icm1', "`col'")
			mat `ICM1'[`i',`j'] = `icm1'[`n',`m']
		}
	}
	
	*** transform data to be consistent with null hypothesis
	mata: _transform("`ICM1'", "`varlist'")
	
	*** bootstrap tetrad test on transformed data and store results
	tempfile bs
	postutil clear
	postfile bs rep T Ts using `bs'
	
	forval i = 1/`reps' {
		_CTABoot `varlist', icm1(`icm1')
		post bs (`i') (`T') (r(Ts))
	}
	
	postclose bs
	
	qui use `bs', clear
	qui gen pv = (Ts > T)
	qui sum pv
	local pv = r(mean)
	
	*** display and return results
	dis as text ""
	dis as text "Confirmatory Tetrad Analysis Results"
	dis as text ""
	dis as text "                     bootstrap"
	dis as text "Chi-sq       df      p-value"
	dis as text "{hline 30}"
	dis as res %8.4f `T' "    " as res %3.0f `df' "    " as res %8.4f `pv'
	dis as text "{hline 30}"

	*** restore data
	qui use `pres', clear
end



*** Bootstrap program
capture program drop _CTABoot
program _CTABoot, rclass
	syntax varlist(min = 4 numeric), icm1(name) 
	preserve
	bsample
	tetrad `varlist', icm1(`icm1') boot(1)
	return scalar Ts = r(T)
	restore
end



*** Mata function to transform data
cap mata: mata drop _transform()
mata:

// transform data
void _transform( string sigma, string vars ) {
	
	// get data from Stata
	st_view(data = ., ., tokens(vars))
	N = rows(data)
	
	// get ICM from Stata
	Sigma = st_matrix(sigma)
	
	// center observed variables and form sample covariance matrix
	means = colsum(data)/N
	S     = ( cross(data,data) - N*means'*means )/(N - 1)
	
	// prepare components to transform data
	SigmaC = cholesky(Sigma)
	SC1    = cholesky(S)
	SC2    = solveupper(SC1', I(rows(S)))
	
	// transform data
	data[.,.] = data*SC2*SigmaC'
	
}

end

/* History
1.3  05.26.15  initial version that corresponds with tetrad v1.3.0 
1.4  06.01.15  update to correspond with tetrad v1.4.0
1.5  06.01.15  update handling of errors
1.6  06.02.15  fixed bug in bootstrap program
*/






	



