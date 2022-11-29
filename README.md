# micmix

Statistical models to determine empirical cut-off points (ECOFF) for antimicrobial resistance from a collection of MIC values. The R function fits a mixture of log-normal distribtions to the reported MIC values by maximum likelihood. The function compares the fit of mixture models with 2, 3, ... k distributions to the data based on Akaike Information Criterion (AIC).
Candidate cutoff points are calculated as the intersects of the normal distributions in order to minimise the predicted probability of misclassification.

The Rmd file shows the calculation applied to S. epidermidis MIC data as presented in the manuscript by X. Ba and collaborators.
