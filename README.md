## HiSSE Analyses for Afrobatrachia


This repository contains the R code and data necessary to conduct the hidden state speciation and extinction analyses of trait-dependent diversification from:

Portik, D.M., Bell, R.C., Blackburn, D.C., Bauer, A.M., Barratt, C.D., Branch, W.R., Burger, M., Channing, A., Colston, T.J., Conradie, W., Dehling, J.M., Drewes, R.C., Ernst, R., Greenbaum, E., Gvoždík, V., Harvey, J., Hillers, A., Hirschfeld, M., Jongsma, G.F.M., Kielgast, J., Kouete, M.T., Lawson, L., Leaché, A.D., Loader, S.P., Lötters, S., van der Meijden, A., Menegon, M., Müller, S., Nagy, Z.T., Ofori-Boateng, C., Ohler, A., Papenfuss, T.J., Rößler, D., Sinsch, U., Rödel, M.-O., Veith, M., Vindum, J., Zassi-Boulou, A.-G., and J.A. McGuire. **2019**. Sexual dichromatism drives diversification within a major radiation of African amphibians. Systematic Biology, 68: 859-875. 

The paper is provided as a PDF in this repository, but can also be found on the publisher's website [here](https://doi.org/10.1093/sysbio/syz023).


The R-script `HiSSE_Afrobatrachia_Dichromatism.R` uses data available in the [data](https://github.com/dportik/HiSSE_for_Afrobatrachia/tree/master/data) folder and also relies on additional R-scripts present in the [dependency-R-scripts](https://github.com/dportik/HiSSE_for_Afrobatrachia/tree/master/dependency-R-scripts) folder. To use the R script, you will need to edit the paths to the various files at the top of the script, and you should also edit the number of cores you have available on your local machine (it is set to 16 cores in the script!).

The complete set of analyses run for this publication (including those here) are available on an Open Science Framework project page ([link](https://osf.io/yeu38/)).