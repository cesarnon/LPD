# LPD - Land Productivity Dynamics
Codes co-developed with WOCAT and FAO projects

# Introduction
The dynamics in the land productivity indicator is related to changes in the health and productive capacity of the land and reflects the net effects of changes in ecosystem functioning due to changes in plant phenology and biomass growth, where declining trends are often (but not always) a defining characteristic of land degradation. Understanding changes in the productive capacity of the land is critical for assessing the impact of land management interventions, its long-term sustainability, and the climate-derived impacts which could affect ecosystem resilience and human livelihoods.
The Land Productivity Dynamics indicator (LPD) is one of the three main indicators to monitor progress towards Land Degradation Neutrality and estimate SDG indicator 15.3.1 (UNCCD 2021). Its calculation is based on the analysis of time series of vegetation indices derived from remote sensed imagery, which are a proxy of the total aboveground net primary production (NPP). 
Different calculations based on the same Earth Observation (EO) source data can produce different results (Teich et al., 2019, Paredes-Trejo et al., 2023), highlighting the importance of integrating EO data with other sources of information, such as experts’ knowledge through participative processes (García et al., 2018, Teich et al.,2023). Also, interpretation of results in the local context and with experts is needed to identify false positives and negatives and the drivers of degradation.
The FAO-WOCAT Land Productivity Dynamics approach
The original algorithm was developed in 2021 in the context of a WOCAT-FAO letter of agreement for the CACILM2-Project in Central Asia. It is based on the methodology applied in the World Atlas of Desertification (Cherlet et al., 2018), functions applied in Trends.Earth (T.E codes),   which were updated using Googe Earth Engine (FAO, 2022). Time-series of annual NDVI from MODIS MOD13Q1 (250m resolution) are analyzed to obtain trends, that are further classified according to their current state, considering an initial biomass. Currently a methodological paper is being developed but more information can be found in Teich & Garcia (2023).
The maps produced then were featured in the:
•	FAO Publication “Overview of land degradation neutrality (LDN) in Europe and Central Asia” (FAO 2022 and App)
•	As an option to use in Trends.Earth Plug-in, tbe UNCCD recommended tool for SDG 15.3.1 calculation (datasets available)
•	The Tools4LDN SDG 15.3.1 Comparison Tool – Released in the UNCCD CRIC 2022 Côte d’Ivoire (https://maps.tools4ldn.org/).
•	Supporting many Convergence of evidence Apps and Decision Support Systems in WOCAT and FAO: https://www.wocat.net/en/ldn/wocatapps/ , https://projectgeffao.users.earthengine.app/ 
•	Contributing to many national processes for UNCCD PRAIS4 reporting, side events at CRIC (2022-2023) and COP (2023), international workshops and publications.
The original version script that is shared in this repository contains the possibility to choose different analysis and to parametrize many sections of the model. To produce a global map of FAO-WOCAT LPD, users need to run the script and will obtain a map with the parameters specified for the default version. Users are encouraged to choose for themselves and parametrize the model to their own area of interest.
Please see the Original model version at: 
[https://code.earthengine.google.com/7103e3406cdc7f560f405a229a2cfda4] (https://code.earthengine.google.com/7103e3406cdc7f560f405a229a2cfda4) 
We are currently working on a publication and we are producing different model versions to add new functionalities and adaptations based on countries request for different ecoregions and land process. We will be sharing our new Derived Models and Next-Generation Model in this repo when they become available.

Citation:  Garcia, C. L., & Teich, I. (2022). FAO-WOCAT Land Productivity Dynamics indicator. Zenodo. https://doi.org/10.5281/zenodo.10849368
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10849368.svg)](https://doi.org/10.5281/zenodo.10849368)


# References
Cherlet M., C. Hutchinson, J. Reynolds, J. Hill, S. Sommer, & G. Maltitz (Eds.). (2018). World atlas of desertification. Publication Office of the European Union. https://doi.org/10.2760/06292 
FAO. (2022b). Overview of land degradation neutrality (LDN) in Europe and Central Asia. FAO. https://doi.org/10.4060/cb7986en 
García, C.L.; Teich, I.; Gonzalez-Roglich, M.; Kindgard, A.F.; Ravelo, A.C.; Liniger, H. Land degradation assessment in the Argentinean Puna: Comparing expert knowledge with satellite-derived information. Environ. Sci. Policy 2019, 91, 70–80. https://doi.org/10.1016/j.envsci.2018.10.018
Paredes-Trejo, F.; Barbosa, H.A.; Daldegan, G.A.; Teich, I.; García, C.L.; Kumar, T.V.L.; Buriti, C.d.O. Impact of Drought on Land Productivity and Degradation in the Brazilian Semiarid Region. Land 2023, 12, 954. https://doi.org/10.3390/land12050954
Teich, I.; Gonzalez Roglich, M.; Corso, M.L.; García, C.L. Combining Earth Observations, Cloud Computing, and Expert Knowledge to Inform National Level Degradation Assessments in Support of the 2030 Development Agenda. Remote Sens. 2019, 11, 2918. https://doi.org/10.3390/rs11242918
Teich I. and Garcia C. 2021. The FAO-WOCAT Land Productivity Map. WOCAT. https://www.wocat.net/documents/1143/FAO_WOCAT_LPD.pdf
Teich, I., Harari, N., Caza, P., Henao-Henao, J. P., Lopez, J. C., Raviolo, E., Díaz-González, A. M., González, H., Bastidas, S., Morales-Opazo, C., & García, C. L. (2023). An interactive system to map land degradation and inform decision-making to achieve land degradation neutrality via convergence of evidence across scales: A case-study in Ecuador. Land Degradation & Development, 34(15), 4475–4487. https://doi.org/10.1002/ldr.4645 
Trends.Earth Codes. https://github.com/ConservationInternational/trends.earth 
UNCCD. 2021. Good Practice Guidance. SDG Indicator 15.3.1, Proportion of Land That Is Degraded Over Total Land Area. Version 2.0. Sims, N.C., Newnham, G.J., England, J.R., Guerschman, J., Cox, S.J.D., Roxburgh, S.H., Viscarra Rossel, R.A., Fritz, S. and Wheeler, I. United Nations Convention to Combat Desertification, Bonn, Germany. https://www.unccd.int/resources/manuals-and-guides/good-practice-guidance-sdg-indicator-1531-proportion-land-degraded 

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
