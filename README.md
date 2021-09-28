# tetrad
This respository contains three Stata commands for conducting Confirmatory Tetrad Analysis (CTA). In addition, the repository is set up so that the commands can be installed within Stata. Simply type "net install tetrad, from(https://github.com/sbauldry/tetrad/raw/master) replace" without the quotes in Stata's command window. These commands are introduced in Bauldry and Bollen (2016).

Note: This command was updated in 2021 to reflect a change in Stata 17 to the column names assigned to Sigma.

## Files
1. tetrad.ado: Primary Stata ado file for conducting CTA.
2. tetrad.sthlp: Stata help file for tetrad command.
3. tetrad_bootstrap.ado: Extension that conducts CTA with a bootstrapped chi-square test statistic.
4. tetrad_bootstrap.sthlp: Stata help file for tetrad_bootstrap.
5. tetrad_matrix.ado: Extension that conducts CTA with matrices rather than raw data.
6. tetrad_matrix.sthlp: Stata help file for tetrad_matrix. 
7. stata.toc: Stata toc file for installing tetrad within Stata
8. tetrad.pkg: Stata pkg file with list of files to install within Stata

## Citation
Bauldry, Shawn and Kenneth A. Bollen. 2016. "tetrad: A Set of Stata Commands for Confirmatory Tetrad Analysis." *Structural Equation Modeling* 23:921-930.
