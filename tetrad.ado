*! v1.6, CTA program, S Bauldry, 02jun2015

capture program drop tetrad
program define tetrad, rclass
	version 13
	syntax varlist(min = 4 numeric) [if] [in], ///
		   icm1(name)       /// implied covariance matrix 1
	       [icm2(name)      /// implied covariance matrix 2
		   reps(integer 1)  /// number of replications
		   seed(string)     /// random number seed
		   tlist(integer 0) /// request for list of vanishing tetrads
		   boot(integer 0)] /// request for bootstrap test statistic
	
	marksample touse
	if ("`seed'" != "") {
		set seed `seed'
	}

	*** if one ICM, set 2nd ICM to 1st for ease of programming
	if ("`icm2'" == "") {
		local Nest = 0
		local icm2 = "`icm1'"
	}
	else if ("`icm2'" != "") {
		local Nest = 1
	}
	
	*** verify that ICMs contain same variables as varlist
	local icm1list1 : colfullnames `icm1'
	local icm1list1 : subinstr local icm1list1 "observed:" "", all
	foreach v1 of local icm1list1 {
		if regexm( "`v1'", "^latent:") {
			local icm1list "`icm1list'"
		}
		else {
			local icm1list "`icm1list' `v1'"
		}
	}
	
	local icm2list1 : colfullnames `icm2'
	local icm2list1 : subinstr local icm2list1 "observed:" "", all
	foreach v2 of local icm2list1 {
		if regexm( "`v2'", "^latent:") {
			local icm2list "`icm2list'"
		}
		else {
			local icm2list "`icm2list' `v2'"
		}
	}

	local var1chck : list varlist === icm1list
	local var2chck : list varlist === icm2list
	if ( `var1chck' == 0 | `var2chck' == 0 ) {
		dis as error "varlist different than ICM variables"
		exit 111
	}
	
	*** obtain sample covariance matrix
	tempname SCM
	qui cor `varlist' if `touse', cov
	mat `SCM' = r(C)
	local N = r(N)
	
	*** put observed variables from ICMs in same order as SCM
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
	
	*** running Mata function for CTA
	mata: _CTA(`N', "`SCM'", "`ICM1'", "`ICM2'", `Nest', `reps', `boot')
	
	*** displaying results for non-nested test
	if (`Nest' != 1 & `boot' != 1) {

		if (`tlist' == 1) {
			local r = rowsof(_ICMTetList)
			
			dis as text ""
			dis as text "Model-Implied Tetrads"
			dis as text "tetrad   residual   AVar       t-value    vanish"
			dis as text "{hline 50}"

			forval i = 1/`r' {
				dis as res %5.0f _ICMTetList[`i',1] "   " ///
				    as res %9.6f _ICMTetList[`i',3] "  " ///
				    as res %9.6f _ICMTetList[`i',4] "  " ///
				    as res %8.4f _ICMTetList[`i',5] "  " ///
                    as res %3.0f _ICMTetList[`i',6]
			}

			dis as text "{hline 50}"
			dis as text ""
			dis as text ""

		}
		
		dis as text ""
		dis as text "Confirmatory Tetrad Analysis Results"
		dis as text "               Model 1"
		dis as text "rep      Chi-sq     df  p-val"
		dis as text "{hline 30}"
		
		forval i = 1/`reps' {
			dis as res %3.0f `i' "     " as res %8.4f _CTAResults[`i',1] "   " ///
			                             as res %3.0f _CTAResults[`i',2] "  " ///
									     as res %5.4f _CTAResults[`i',3]
		}
		dis as text "{hline 30}"
		
	}
	
	
	
	if (`Nest' == 1 & `boot' != 1) {

		if (`tlist' == 1) {
			local r = rowsof(_ICMTetList)
			
			dis as text ""
			dis as text "Model-Implied Tetrads"
			dis as text "                      Model 1"
			dis as text "tetrad   residual   AVar       t-value    vanish"
			dis as text "{hline 50}"

			forval i = 1/`r' {
				dis as res %5.0f _ICMTetList[`i',1] "   " ///
				    as res %9.6f _ICMTetList[`i',3] "  " ///
				    as res %9.6f _ICMTetList[`i',4] "  " ///
				    as res %8.4f _ICMTetList[`i',5] "  " ///
                    as res %3.0f _ICMTetList[`i',6]
			}

			dis as text "{hline 50}"
			dis as text ""
			dis as text ""
			
			
						dis as text "Model-Implied Tetrads"
			dis as text "                      Model 2"
			dis as text "tetrad   residual   AVar       t-value    vanish"
			dis as text "{hline 50}"

			forval i = 1/`r' {
				dis as res %5.0f _ICMTetList[`i',7] "   " ///
				    as res %9.6f _ICMTetList[`i',9] "  " ///
				    as res %9.6f _ICMTetList[`i',10] "  " ///
				    as res %8.4f _ICMTetList[`i',11] "  " ///
                    as res %3.0f _ICMTetList[`i',12]
			}

			dis as text "{hline 50}"
			dis as text ""
			dis as text ""

		}
		
		dis as text ""
		dis as text "Confirmatory Tetrad Analysis Results"
		dis as text "              Model 1                  Model 2                  M1 - M2"
		dis as text "rep     Chi-sq    df  p-val      Chi-sq    df  p-val      Chi-sq    df  p-val"
		dis as text "{hline 80}"
		
		forval i = 1/`reps' {
			dis as res %3.0f `i' "    " as res %8.4f _CTAResults[`i',1] "  " ///
			                            as res %3.0f _CTAResults[`i',2] "  " ///
									    as res %5.4f _CTAResults[`i',3] "    " ///
									    as res %8.4f _CTAResults[`i',4] "  " ///
			                            as res %3.0f _CTAResults[`i',5] "  " ///
									    as res %5.4f _CTAResults[`i',6] "    " ///
										as res %8.4f _CTAResults[`i',7] "  " ///
			                            as res %3.0f _CTAResults[`i',8] "  " ///
									    as res %5.4f _CTAResults[`i',9] "  "
		}
		dis as text "{hline 80}"
	}
	
	*** returning results
	matrix colnames _CTAResults = chi1 df1 pval1 chi2 df2 pval2 chi3 df3 pval3
	if (`Nest' != 1) {
		matrix _CTAResults = _CTAResults[.,1..3]
	}
	
	if (`boot' == 1) {
		return scalar T = _CTAResults[1,1]
		return scalar df = _CTAResults[1,2]
	}
	return matrix _CTAResults = _CTAResults

end



*** Clear mata functions for this program
cap mata: mata drop _CTA() _IdVT() _IdNRVT() _TestStatistic() ///
                    _BSTestStatistic() _ListTetrad() _TransformICM() ///
					_Tetrads() _TetVar() _checkVT() _UniqueCov() _GrdMat() ///
					_NTAsympCov() _sweep() _sweep_all() 


*** CTA Main Mata function
mata:

void _CTA( real scalar N, string m1, string m2, string m3, ///
           real scalar Nest, real scalar reps, real scalar boot ) {

	// Read Stata matrices
	_SCM  = st_matrix(m1)
	_ICM1 = st_matrix(m2)
	_ICM2 = st_matrix(m3)

	// Transform ICMs for identifying vanishing and nonredudant tetrads
	_ICM1t = _TransformICM( _ICM1 )
	_ICM2t = _TransformICM( _ICM2 )

	// Step 1: Identify vanishing tetrads
	_ICM1VT = _IdVT(_SCM, _ICM1t, N)
	_ICM2VT = _IdVT(_SCM, _ICM2t, N)

		// Check for 0 vanishing tetrads and model 2 tetrad-nested in model 1
		_checkVT(_ICM1VT, _ICM2VT)

	// Step 2: Initiate replications (skip if bootstrapping test statistic)
	if (boot == 0) {
	
		_CTAResults = J(reps, 9, .)
		for (r = 1; r <= reps; r++) {
		
			// randomize order of vanishing tetrads
			_rn1 = jumble( 1::rows(_ICM1VT) )
			_rn2 = _rn1[1..rows(_ICM2VT), .]
			_ICM1VTr = sort( (_ICM1VT, _rn1), 14 )
			_ICM2VTr = sort( (_ICM2VT, _rn2), 14 )

			// Step 3A: Identify set of nonredundant vanishing tetrads
			_ICM1NRVT = _IdNRVT(_ICM1VTr, _ICM1t)
			_ICM2NRVT = _IdNRVT(_ICM2VTr, _ICM2t)
	
			// Step 4A: Calculate test statistics
			_CTAResults[r,1..3] = _TestStatistic( _ICM1NRVT, _SCM, N )
			_CTAResults[r,4..6] = _TestStatistic( _ICM2NRVT, _SCM, N )
		
				// test statistics for nested test
				_CTAResults[r,7] = _CTAResults[r,1] - _CTAResults[r,4]
				_CTAResults[r,8] = _CTAResults[r,2] - _CTAResults[r,5]
				_CTAResults[r,9] = chi2tail( _CTAResults[r,8], _CTAResults[r,7])
		}
	}

	else if (boot == 1) {
		
		// Step 4B: calculate test statistic
		_CTAResults = J(1, 9, .)
		_CTAResults[1,1..3] = _BSTestStatistic( _ICM1VT, _SCM, N )
		
	}
	
	// Step 5: Return results to Stata
		
		// List of tetrads
		_ICM1List = _ListTetrad(_ICM1, _ICM1t, N)
		_ICM2List = _ListTetrad(_ICM2, _ICM2t, N)
		_ICMTetList = _ICM1List, _ICM2List
	
	st_matrix("_ICMTetList", _ICMTetList)
	st_matrix("_CTAResults", _CTAResults)
}

	
	
	

// Step 1: Identify model-implied vanishing tetrads
real matrix _IdVT( real matrix scm, real matrix icm, real scalar N ) {

	// create labels and calculate tetrads
	st = _Tetrads(scm)
	it = _Tetrads(icm)

	// calculate covariance matrix for ICM tetrads
	ic = _TetVar( it, icm, N )

	// identify implied vanishing tetrads
	tmp1 = st, it[.,10], ic, it[.,10] :/ sqrt(ic)
	vt = J( 1, cols(tmp1), . )
	for (i = 1; i <= rows(tmp1); i++) {
		if ( abs( round( tmp1[i,13], 0.001) ) < 0.001 ) {
			vt = vt \ tmp1[i,.]
		}
	}

	// return matrix
	if (rows(vt) == 1) {
		return(vt)
	}
	else {
		vt = vt[|2,1 \ .,13|]
		return(vt)
	}
}



// Step 3A: Identify set of nonredundant vanishing tetrads
real matrix _IdNRVT( real matrix vt, real matrix icm ) {

	 // if only one tetrad, set to set of nonredundant tetrads
	 if ( rows(vt) == 1 ) {
		nrvt = vt
	 }

	 else if ( rows(vt) > 1 ) {
		
		// find unique covariances among vanishing tetrads
		lvt = _UniqueCov(vt)
		uc  = uniqrows( vec(lvt) )'

		// calculate gradient matrix
		dvt = _GrdMat( uc, lvt, vt, icm )

		// calculate asymptotic covariance matrix for unique covariances
		acm = _NTAsympCov( uc, icm )

		// calculate asymptotic covariance matrix for ICM vanishing tetrads
		avt = dvt'acm*dvt

		// use sweep operator to identify set of nonredundant vanishing tetrads
		tmp1 = _sweep_all(avt)
		nrvt  = J(1, cols(vt), .)
		
		for (i = 1; i <= rows(vt); i++) {
			if ( rowsum(tmp1[i,1]) != 0 ) {
				nrvt = nrvt \ vt[i,.]
			}
		}
		nrvt = nrvt[|2,1 \ .,13|]
	}
	
	return(nrvt)
}


// Step 4A: Calculate test statistic
real vector _TestStatistic( real matrix nrvt, real matrix scm, real scalar N) {
	
	// find unique covariances among nonredundant vanishing tetrads
	lnrvt = _UniqueCov(nrvt)
	uc    = uniqrows( vec(lnrvt) )'
	
	// calculate gradient matrix
	dvt = _GrdMat( uc, lnrvt, nrvt, scm )
	
	// calculate asymptotic covariance matrix for unique covariances
	acm = _NTAsympCov( uc, scm )
	
	// calculate asymptotic covariance matrix for SCM nonredundant VTs
	avt = dvt'*acm*dvt
	
	// store chi-square, df, and p-value in vector
	ts = J(1, 3, .)
	ts[1,1] = N*( nrvt[.,10]'*luinv(avt)*nrvt[.,10] )
	ts[1,2] = rows(avt)
	ts[1,3] = chi2tail( ts[1,2], ts[1,1] )
	
	return(ts)
}



// Step 4B: Calculate test statistic if bootstrapping
real vector _BSTestStatistic( real matrix vt, real matrix scm, real scalar N) {
	
	// find unique covariances among vanishing tetrads
	lvt = _UniqueCov(vt)
	uc    = uniqrows( vec(lvt) )'
	
	// calculate gradient matrix
	dvt = _GrdMat( uc, lvt, vt, scm )
	
	// calculate asymptotic covariance matrix for unique covariances
	acm = _NTAsympCov( uc, scm )
	
	// calculate asymptotic covariance matrix for SCM nonredundant VTs
	avt = dvt'*acm*dvt
	
	// create diagonal matrix
	tmp1 = diagonal(avt)
	davt = diag(tmp1)
	
	// store chi-square, df, and p-value in vector
	ts = J(1, 3, .)
	ts[1,1] = N*( vt[.,10]'*luinv(davt)*vt[.,10] )
	ts[1,2] = rows(davt)
	ts[1,3] = chi2tail( ts[1,2], ts[1,1] )
	
	return(ts)
}



// Prepare list of tetrads
real matrix _ListTetrad( real matrix icm, real matrix icmt, real scalar N ) {

	// create labels and calculate tetrads
	icm_tet  = _Tetrads(icm)
	icmt_tet = _Tetrads(icmt)

	// calculate covariance matrix for ICM tetrads
	ic = _TetVar( icmt_tet, icmt, N )
	
	// prepare list of tetrads
	tet = icm_tet[.,1], icm_tet[.,10], icmt_tet[.,10], ic, icmt_tet[.,10] :/ sqrt(ic), J( rows(icm_tet), 1, .)
	for (i = 1; i <= rows(tet); i++) {
		if ( abs( round( tet[i,5], 0.001) ) < 0.001 ) {
			tet[i,6] = 1
		}
		else {
			tet[i,6] = 0
		}
	}
	
	return(tet)
}
	



// Transform ICM 
real matrix _TransformICM(real matrix m1) { 

	real matrix m2

	TempDiag = J(cols(m1), cols(m1), 0)
	for (i = 1; i <= cols(m1); i++) {
		TempDiag[i,i] = 1/sqrt( m1[i,i] )
	}
	m2 = TempDiag*m1*TempDiag

	return(m2)
}



// Create labels and calculate tetrad residuals
real matrix _Tetrads(real matrix m1) {

	m2 = J(1, 10, .)
	for (a = 1; a <= cols(m1); a++) {
		for (b = a + 1; b <= cols(m1); b++) {
			for (c = b + 1; c <= cols(m1); c++) {
				for (d = c + 1; d <= cols(m1); d++) {
					m2 = m2 \ ( a*1000 + b*100 + c*10 + d, a, b, c, d, a, c, b, d, m1[a,b]*m1[c,d] - m1[a,c]*m1[b,d] )
					m2 = m2 \ ( a*1000 + c*100 + d*10 + b, a, c, d, b, a, d, c, b, m1[a,c]*m1[d,b] - m1[a,d]*m1[c,b] )
					m2 = m2 \ ( a*1000 + d*100 + b*10 + c, a, d, b, c, a, b, d, c, m1[a,d]*m1[b,c] - m1[a,b]*m1[d,c] )
				}
			}
		}
	}
	
	m2 = m2[|2,.\.,.|]
	return(m2)
}



// Calculate variance of tetrads
real vector _TetVar(real matrix m1, real matrix m2, real scalar N) {
	
	m3 = J(1, 1, .)

	for (i = 1; i <= rows(m1); i++) {
		a = m1[i,2]
		b = m1[i,3]
		c = m1[i,4]
		d = m1[i,5]
	
		if ( mod(i,3) == 0 | mod(i,3) == 2 ) {
			temp = c
			c = d
			d = temp
		}
	
		m3 = m3 \ 1/N*( m2[b,d]^2*m2[a,a]*m2[c,c] + ///
						m2[a,c]^2*m2[d,d]*m2[b,b] + ///
						m2[c,d]^2*m2[a,a]*m2[b,b] + ///
						m2[a,b]^2*m2[d,d]*m2[c,c] + ///
						2*( m2[b,d]*m2[a,c]*m2[a,d]*m2[b,c] - ///
							m2[b,d]*m2[c,d]*m2[a,a]*m2[b,c] - ///
							m2[b,d]*m2[a,b]*m2[a,d]*m2[c,c] - ///
							m2[a,c]*m2[c,d]*m2[a,d]*m2[b,b] - ///
							m2[a,c]*m2[a,b]*m2[d,d]*m2[b,c] + ///
							m2[c,d]*m2[a,b]*m2[a,d]*m2[b,c]) + ///
						2*( (m2[b,d]*m2[a,c] - m2[a,b]*m2[c,d])^2 ) )
	}

	m3 = m3[|2,.\.,.|]
	return(m3)
}



// Check whether 0 vanishing tetrads and tetrad-nested
void _checkVT( real matrix m1, real matrix m2 ) {
		if ( m1[1,1] == . ) {
		 	_error( "Model 1 has no implied vanishing tetrads." )
	    }

	    if ( m2[1,1] == . ) {
		 	_error("Model 2 has no implied vanishing tetrads.")
	    }

	    tet1 = uniqrows( m1[1..rows(m1),1] )
	    tet2 = uniqrows( m2[1..rows(m2),1] )
	    
	    for (i = 1; i <= rows(tet2); i++) {
	    	tet = tet2[i]
	    	if ( anyof(tet1, tet) == 0 ) {
	    		_error("Model 2 is not tetrad-nested in Model 1.")
	    	}
	    }
}



// Find unique covariances among vanishing tetrads
string matrix _UniqueCov(real matrix m1) {
	
	string matrix m2
	m2 = J( rows(m1), 4, "aaaa" )
	
	for (i = 1; i <= rows(m1); i++) {
		k = 0
		for (j = 2; j <= 8; j = j + 2) {
			k = k + 1
			if ( m1[i,j] < 10 & m1[i,j+1] < 10 ) {
				m2[i,k] = "0" + strofreal(m1[i,j+1]) + "0" + strofreal(m1[i,j])
			}
			else if ( m1[i,j] >= 10 & m1[i,j+1] >= 10 ) {
				m2[i,k] = strofreal(m1[i,j+1]) + strofreal(m1[i,j])
			}
			else if ( m1[i,j] >= 10 ) {
				m2[i,k] = "0" + strofreal(m1[i,j+1]) + strofreal(m1[i,j])
			}
			else {
				m2[i,k] = strofreal(m1[i,j+1]) + "0" + strofreal(m1[i,j])
			}
		}
	}
	
	return(m2)
}



// calculate gradient matrix
real matrix _GrdMat( string vector v1, string matrix m1, ///
                     real matrix m2, real matrix m3 ) {

	real matrix m4
	m4 = J(cols(v1), rows(m1), 0 )
	
	for (i = 1; i <= cols(v1); i++) {
		for (j = 1; j <= rows(m1); j++) {
			if ( v1[1,i] == m1[j,1] ) {
				m4[i,j] = m3[ m2[j,4], m2[j,5] ]
			}
			if ( v1[1,i] == m1[j,2] ) {
				m4[i,j] = m3[ m2[j,2], m2[j,3] ]
			}
			if ( v1[1,i] == m1[j,3] ) {
				m4[i,j] = -1*m3[ m2[j,8], m2[j,9] ]
			}
			if ( v1[1,i] == m1[j,4] ) {
				m4[i,j] = -1*m3[ m2[j,6], m2[j,7] ]
			}
		}
	}
	
	return(m4)
}



// calculate normal theory asymptotic covariance matrix for unique covariances
real matrix _NTAsympCov(string vector v1, real matrix m1) {

	real matrix m2
	m2 = J(cols(v1), cols(v1), 0)
	
	for (i = 1; i <= cols(v1); i++) {
		for (j = 1; j <= cols(v1); j++) {
			e = strtoreal( substr( v1[1,i], 1, 2 ) )
			f = strtoreal( substr( v1[1,i], 3, 2 ) )
			g = strtoreal( substr( v1[1,j], 1, 2 ) )
			h = strtoreal( substr( v1[1,j], 3, 2 ) )
		
			// covariance matrix
			if ( trace(m1) != cols(m1) ) {
				m2[i,j] = m1[e,g]*m1[f,h] + m1[e,h]*m1[f,g]
			}
			
			// correlation matrix
			if ( trace(m1) == cols(m1) ) {
				m2[i,j] = (1/2)*m1[e,f]*m1[g,h]* ///
				          (m1[e,g]^2 + m1[e,h]^2 + m1[f,g]^2 + m1[f,h]^2) + ///
						  m1[e,g]*m1[f,h] + m1[e,h]*m1[f,g] - ///
                   		  m1[e,f]*(m1[f,g]*m1[f,h] + m1[e,g]*m1[e,h]) - ///
                   		  m1[g,h]*(m1[f,g]*m1[e,g] + m1[f,h]*m1[e,h])
			}
		}
	}

	return(m2)
}



// Sweep operator
// Goodnight (1979) recommends tolerances between 1e-8 and 1e-12 when looking
// for dependencies (pg 156); epsilon(1e4) = 2.2e-12, epsilon(1e8) = 2.2e-8
real matrix _sweep( real matrix A, real scalar k ) {

	real matrix B
	B = J( rows(A), cols(A), 0)

	if ( abs( A[k,k] ) < epsilon(1e6) ) {
		B = A
		for (i = 1; i <= rows(B); i++) {
			B[i,k] = 0
			B[k,i] = 0
		}
	}

	else if ( A[k,k] != 0 ) {

		for (i = 1; i <= rows(A); i++) {
			for (j = 1; j <= cols(A); j++) {
		
				if (i == k & j == k) {
					B[i,j] = 1/A[k,k]
				}
		
				else if (i == k & j != k) {
					B[i,j] = A[k,j]/A[k,k]
				}
		
				else if (i != k & j == k) {
					B[i,j] = -A[i,k]/A[k,k]
				}
		
				else if (i != k & j != k) {
					B[i,j] = A[i,j] - ( A[i,k]*A[k,j] / A[k,k] )
				}
			}
		}
	}

	return(B)

}



// Sweep recursion for all rows
real matrix _sweep_all(real matrix A) {

	real matrix B
	B = J( rows(A), cols(A), 0)

	for (i = 1; i <= rows(A); i++) {
		B = _sweep(A, i)
		A = B
	}

	return(B)

}

end




/* History
1.0  02.02.15  initial version of program
1.1  05.26.15  reorganized primary mata function into smaller functions
1.2  05.26.15  added option to output list of tetrads
1.3  05.25.15  developed bootstrap component
1.4  06.01.15  added test to ensure same variables in ICMs and SCM
1.5  06.01.15  updated handling of errors to accomodate bootstrapping
1.6  06.02.15  fixed bug in bootstrap program
*/






	



