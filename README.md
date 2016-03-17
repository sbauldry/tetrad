# tetrad
This respository contains three Stata commands for conducting Confirmatory Tetrad Analysis (CTA). In addition, the repository is set up so that the commands can be installed within Stata. Simply type "net install tetrad, from(https://github.com/sbauldry/tetrad/raw/master) replace" without the quotes in Stata's command window. These commands are currently in development.

## Files
1. tetrad.ado: Primary Stata ado file for conducting CTA.
2. tetrad.sthlp: Stata help file for tetrad command.
3. tetrad_bootstrap.ado: Extension that conducts CTA with a bootstrapped chi-square test statistic.
4. tetrad_bootstrap.sthlp: Stata help file for tetrad_bootstrap.
5. tetrad_matrix.ado: Extension that conducts CTA with matrices rather than raw data.
6. tetrad_matrix.sthlp: Stata help file for tetrad_matrix. 
7. stata.toc: Stata toc file for installing tetrad within Stata
8. tetrad.pkg: Stata pkg file with list of files to install within Stata
