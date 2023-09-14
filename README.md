# micmix

Statistical models to determine empirical cut-off points (ECOFF) for antimicrobial resistance from a collection of MIC values. The R function fits a mixture of log-normal distribtions to the reported MIC values by maximum likelihood. The function compares the fit of mixture models with 2, 3, ... k distributions to the data based on Akaike Information Criterion (AIC).
Candidate cutoff points are calculated as the intersects of the normal distributions in order to minimise the predicted probability of misclassification.

The Rmd file shows the calculation applied to S. epidermidis MIC data as presented in the manuscript by X. Ba and collaborators: 

"Cryptic susceptibility to penicillins and β-lactamase inhibitors in emerging multidrug-resistant, hospital-adapted Staphylococcus epidermidis lineages"

Xiaoliang Ba(1), Claire L. Raisen(1), Olivier Restif(1), Lina Maria Cavaco(2), Carina Vingsbo Lundberg(2), Jean Y. H. Lee(3), Benjamin P. Howden(3), Mette D. Bartels(4,5), Birgit Strommenger(6), Ewan M. Harrison(7,8,9), Anders Rhod Larsen(2), Mark A. Holmes(1) & Jesper Larsen(2)

1 Department of Veterinary Medicine, University of Cambridge, Cambridge, UK.
2 Department of Bacteria, Parasites & Fungi, Statens Serum Institut, Copenhagen, Denmark.
3 Department of Microbiology and Immunology, The University of Melbourne at The Doherty Institute for Infection and Immunity, Melbourne, Australia.
4 Department of Clinical Microbiology, Copenhagen University Hospital - Amager and Hvidovre, Hvidovre, Denmark.
5 Department of Clinical Medicine, University of Copenhagen, Copenhagen, Denmark.
6 National Reference Centre for Staphylococci and Enterococci, Division Nosocomial Pathogens and Antibiotic Resistances, Department of Infectious Diseases, Robert Koch Institute, Wernigerode Branch, Wernigerode, Germany.
7 Department of Medicine, University of Cambridge, Cambridge, UK.
8 Department of Public Health and Primary Care, University of Cambridge, Cambridge, UK.
9 Wellcome Sanger Institute, Hinxton, UK.

