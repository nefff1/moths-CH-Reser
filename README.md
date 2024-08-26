# Analyses of moth community dataset L. Rezbanyai-Reser Switzerland

This repository contains codes and some data that were used in the analyses for the following manuscripts, which are based on a vast moth community dataset from Switzerland collected by Dr. Ladislaus Rezbanyai-Reser:

-   Neff F, Chittaro Y, Korner-Nievergelt F, Litsios G, Rey E, Knop E. **Bimodal seasonal activity of moths and elevation, weather and land use as drivers of their diversity** (referred to as *manuscript 1*)

-   Neff F, Chittaro Y, Korner-Nievergelt F, Litsios G, Martínez-Núñez C, Rey E, Herzog F, Knop E. **Changes in moth communities across 50 years depend on elevation** (referred to as *manuscript 2*)

The following R files are included in the folder *R_Code*:

-   **R_prepare_data.R**: File used to prepare all data frames used for the analyses of the two manuscripts.

-   **R_moth_models_1.R**: Code used for the analyses of manuscript 1.

-   **R_moth_models_2.R**: Code used for the analyses of manuscript 2.

The following Stan code files are included in the folder *Stan_Code*:

-   **Stan_hg_spline_s1_r4.stan**: Stan model code of a regression model with 1 smoothing term, 4 random terms and a hurdle gamma distribution.

-   **Stan_hg_spline_s1p1_r4.stan**: Stan model code of a regression model with 1 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a hurdle gamma distribution.

-   **Stan_hg_spline_s2_r4.stan**: Stan model code of a regression model with 12 smoothing term, 4 random terms and a hurdle gamma distribution.

-   **Stan_hg_spline_s2p1_r4.stan**: Stan model code of a regression model with 2 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a hurdle gamma distribution.

-   **Stan_nb_spline_s1_r4.stan**: Stan model code of a regression model with 1 smoothing term, 4 random terms and a zero-inflated negative binomial distribution.

-   **Stan_nb_spline_s1p1_r4.stan**: Stan model code of a regression model with 1 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a zero-inflated negative binomial distribution.

-   **Stan_nb_spline_s2_r4.stan**: Stan model code of a regression model with 12 smoothing term, 4 random terms and a zero-inflated negative binomial distribution.

-   **Stan_nb_spline_s2p1_r4.stan**: Stan model code of a regression model with 2 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a zero-inflated negative binomial distribution.

Besides data available from the manuscripts directly or from other sources indicated in the code (e.g. GBIF), the folder *Data* contains:

-   **d_samplings.txt**: Details on sampling site-year combinations used in the analyses: Spatio-temporal clusters, sampling pairs, different land-use proportions.
-   **d_nullnights.txt**: List of nights in which nothing was caught (missing from GBIF dataset).
-   **d_taxonomy.txt**: Moth species names according to the taxonomy used in the current analyses. Can be joint to GBIF data through the *taxonID* variable.
-   **d_mass.txt**: Estimated species-level biomass used to estimate community-wide biomass. Based on a set of allometric relationships, see *Manuscript 1*.