# lpj2char

This repository contains a helper program to convert burned area fraction and carbon emissions from wildfire simulated with the LPJ-LMfire dynamic global vegetation model into Z-scores to compare with paleo-charcoal data. Transformation of model output into z-scores follows the procedure outlined in Brücher et al (2014); see also Marlon et al (2016) and further references at https://pjbartlein.github.io/GCDv3Analysis/. The fortran code for calculating the z-score lambda was adapted from code developed by P.J. Bartlein and the R MASS library.

References

Marlon, J. R., Kelly, R., Daniau, A.-L., Vannière, B., Power, M. J., Bartlein, P., Higuera, P., Blarquez, O., Brewer, S., Brücher, T., Feurdean, A., Romera, G. G., Iglesias, V., Maezumi, S. Y., Magi, B., Courtney Mustaphi, C. J., & Zhihai, T. (2016). Reconstructions of biomass burning from sediment-charcoal records to improve data–model comparisons. Biogeosciences, 13(11), 3225-3244. doi:10.5194/bg-13-3225-2016

Brücher, T., Brovkin, V., Kloster, S., Marlon, J. R., & Power, M. J. (2014). Comparing modelled fire dynamics with charcoal records for the Holocene. Climate of The Past, 10(2), 811-824. doi:10.5194/cp-10-811-2014

