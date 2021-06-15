# PSHREG
## A SAS macro for Proportional Subdistribution Hazards Regression


We present a new SAS macro %PSHREG (https://raw.githubusercontent.com/georgheinze/PSHREG/main/pshreg.sas) that can be used to fit a proportional subdistribution hazards (PSH) model (Fine and Gray, 1999) for survival data subject to competing risks. The macro is described in our accompanying publication (Kohl et al, 2015).

Our macro first modifies the input data set appropriately and then applies SAS's standard Cox regression procedure, PROC PHREG, using weights and counting-process style of specifying survival times to the modified data set (Geskus, 2011). The modified data set can also be used to estimate cumulative incidence curves for the event of interest. The application of PROC PHREG has several advantages, e. g., it directly enables the user to apply the Firth correction, which has been proposed as a solution to the problem of undefined (infinite) maximum likelihood estimates in Cox regression, frequently encountered in small sample analyses (Heinze and Schemper, 2001).

In case of non-PSH, the PSH model is misspecified, but offers a time-averaged summary estimate of the effect of a covariate on the subdistribution hazard (Grambauer, Schumacher and Beyersmann, 2010). Random censoring usually distorts this summary estimate compared to its expected value had censoring not occured, as later event times are underrepresented due to earlier censorship. The solution would be upweighting late event times in the estimating equations by the inverse probability of being observed, similarly to Xu and O'Quigley's (2000) proposal for reweighting the summands of the estimating equations in the Cox model. A very appealing interpretation of the average subdistribution hazard ratio as "odds of concordance" can be obtained by weighting the summands by the expected number of patients at risk (Schemper, Wakounig and Heinze, 2009). Both types of weights are available in %PSHREG. More information can be found in our paper (Kohl et al, 2015) and in a technical report (folder <https://github.com/georgheinze/PSHREG/tree/main/docs>).

## References:

* Kohl, M., Plischke,M., Leffondré,K., Heinze,G. (2015): "PSHREG: A SAS macro for proportional and nonproportional subdistribution hazards regression", Computer Methods and Programs in Biomedicine (doi: <https://dx.doi.org/10.1016/j.cmpb.2014.11.009>)
* Fine, J.P., Gray, R.J. (1999): "A proportional hazards model for the subdistribution of a competing risk", Journal of the American Statistical Association 94(446), 496 - 509
* Geskus, R.B. (2011): "Cause-Specific Cumulative Incidence Estimation and the Fine and Gray Model Under Both Left Truncation and Right Censoring", Biometrics 67, 39 - 49
* Grambauer, N., Schumacher, M., Beyersmann, J. (2010): "Proportional subdistribution hazards modeling offers a summary analysis, even if misspecified", Statistics in Medicine 29, 875-884
* Heinze, G., Schemper, M. (2001): "A solution to the problem of monotone likelihood in Cox regression", Biometrics 57, 114 - 119
* Schemper, M., Wakounig, S., Heinze, G. (2009): "The estimation of average hazard ratios by weighted Cox regression", Statistics in Medicine 28(19), 2473 - 2489
* Xu, R., O'Quigley, J. (2000): "Estimating average regression effect under non-proportional hazards", Biostatistics 1, 423 - 439
