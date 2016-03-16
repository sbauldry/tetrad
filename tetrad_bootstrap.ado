*! v1.9, CTA program, S Bauldry, 14mar2016

program define tetrad_bootstrap, rclass
	version 13
	syntax varlist(min = 4 numeric) [if] [in], ///
		   icm1(name)          /// implied covariance matrix 1
		   [icm2(name)         /// implied covariance matrix 2
		   seed(string)        /// random number seed
		   reps(integer 1000)] /// number of bootstrap replications

	marksample touse
	if ("`seed'" != "") {
		set seed `seed'
	}
	
	*** run tetrad on original data to obtain baseline test statistics
	if ("`icm2'" == "") {
		local Nest = 0
		local icm2 = "`icm1'"
		tetrad `varlist' if `touse', icm1(`icm1') seed(`seed') boot(1)
		mat tb   = r(_CTAResults)
		local T  = tb[1,1]
		local df = tb[1,2]
	}
	else if ("`icm2'" != "") {
		local Nest = 1
		tetrad `varlist' if `touse', icm1(`icm1') icm2(`icm2') seed(`seed') boot(1)
		mat tb   = r(_CTAResults)
		local T  = tb[1,7]
		local df = tb[1,8]
	}
	
	*** save data before transforming
	preserve
	tempfile pres
	qui save `pres'
	
	qui keep `varlist' __000000
	qui keep if `touse'
	
	recast double _all 
	
	*** put observed variables from ICM in same order as SCM
	tempname ICM1 ICM2
	local k : list sizeof varlist
	mat `ICM1' = J(`k', `k', .)
	mat `ICM2' = J(`k', `k', .)
	foreach row of local varlist {
		foreach col of local varlist {
			local i : list posof "`row'" in varlist
			local j : list posof "`col'" in varlist
			local n = rownumb(`icm1', "`row'")
			local m = colnumb(`icm1', "`col'")
			local o = rownumb(`icm2', "`row'")
			local p = colnumb(`icm2', "`col'")
			mat `ICM1'[`i',`j'] = `icm1'[`n',`m']
			mat `ICM2'[`i',`j'] = `icm2'[`o',`p']
		}
	}
	
	*** transform data to be consistent with null hypothesis for most restrictive model
	mata: _transform("`ICM1'", "`varlist'")
	
	*** bootstrap tetrad test on transformed data and store results
	if (`Nest' == 0) {
		tempfile _tetrad_bootstrap
		postfile bs rep T Ts using `_tetrad_bootstrap'
	
		forval i = 1/`reps' {
			_CTABoot `varlist', icm1(`icm1')
			post bs (`i') (`T') (r(Ts))
		}
	
		postclose bs
	
		qui use `_tetrad_bootstrap', clear
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
	}
	
		*** bootstrap tetrad test on transformed data and store results
	if (`Nest' == 1) {
		tempfile _tetrad_bootstrap
		postfile bs rep T Ts using `_tetrad_bootstrap'
	
		forval i = 1/`reps' {
			_CTANestBoot `varlist', icm1(`icm1') icm2(`icm2')
			post bs (`i') (`T') (r(Ts))
		}
	
		postclose bs
	
		qui use `_tetrad_bootstrap', clear
		qui gen pv = (Ts > T)
		qui sum pv
		local pv = r(mean)
	
		*** display and return results
		dis as text ""
		dis as text "Nested Confirmatory Tetrad Analysis Results"
		dis as text ""
		dis as text "M1 - M2              bootstrap"
		dis as text "Chi-sq       df      p-value"
		dis as text "{hline 30}"
		dis as res %8.4f `T' "    " as res %3.0f `df' "    " as res %8.4f `pv'
		dis as text "{hline 30}"
	}
	
	*** return bootstrap p-value
	return scalar _bsp = `pv'
		
	*** restore data
	qui use `pres', clear
end



*** Bootstrap programs
program _CTABoot, rclass
	syntax varlist(min = 4 numeric), icm1(name) 
	preserve
	bsample
	tetrad `varlist', icm1(`icm1') boot(1)
	mat tb    = r(_CTAResults)
	local Ts  = tb[1,1]
	return scalar Ts = `Ts'
	restore
end

program _CTANestBoot, rclass
	syntax varlist(min = 4 numeric), icm1(name) icm2(name)
	preserve
	bsample
	tetrad `varlist', icm1(`icm1') icm2(`icm2') boot(1)
	mat tb    = r(_CTAResults)
	local Ts  = tb[1,7]
	return scalar Ts = `Ts'
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
1.7  09.03.15  set default bootstrap replications to 1000
1.8  02.18.16  removed postutil clear
1.9  03.14.16  updated to allow for nested tests
*/
