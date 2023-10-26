# Shark-management strategies effect-size analyses
<a href="https://doi.org/10.5281/zenodo.10043996"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10043996.svg" alt="DOI"></a>
<img align="right" src="www/sharkbitenet.png" alt="shark bite icon" width="280" style="margin-top: 20px">

Analyses to assess probability of detecting differences in the number of shark bites on humans with/without shark management strategies in place.
<a href="https://github.com/cjabradshaw/AustralianSharkIncidentDatabase"><img align="left" src="www/ASIDlogo3.png" alt="ASID logo" width="150" style="margin-top: 20px"></a>

See also the <a href="https://github.com/cjabradshaw/AustralianSharkIncidentDatabase">Australian Shark-Incident Database</a> (ASID).

<br>
<br>
Prof <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
July 2023/updated October 2023 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>
contributors: <a href="https://www.flinders.edu.au/people/charlie.huveneers">Charlie Huveneers</a>

## <a href="https://github.com/cjabradshaw/SharkManagementStrategiesPower/tree/main/scripts">Scripts</a>
- <code>SMSeffectSizeAnalysis.R</code>: main R code for analysis
- <code>new_lmer_AIC_tables3.R</code>: source code for information-theoretic algorithms
- <code>r.squared.R</code>: source code for calculating goodness-of-fit for linear models (including mixed-effects models)

## <a href="https://github.com/cjabradshaw/SharkManagementStrategiesPower/tree/main/data">Data</a>
- <em>beachmesh.csv</em>: data of shark bites at beaches with and without beach-mess protection (all interactions)
- <em>beachmeshonlybites.csv</em>: data of shark bites at beaches with and without beach-mess protection (bites only)
- <em>sms.csv</em>: data describing shark bites across different beaches before and after shark-management strategies


## Required R packages
- <code>performance</code>
- <code>sjPlot</code>
- <code>lme4</code>
- <code>doSNOW</code>
- <code>snow</code>
- <code>iterators</code>
- <code>foreach</code>
- <code>parallel</code>

<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="200" style="margin-top: 20px"></a> <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="200" style="margin-top: 20px"></a> <a href="https://twitter.com/SouthernSharkEG"><img align="bottom-left" src="www/SSEG.png" alt="SSEG logo" width="150" style="margin-top: 20px"></a> <a href="https://github.com/cjabradshaw/AustralianSharkIncidentDatabase"><img align="bottom-left" src="www/asid shark.gif" alt="ASID shark" width="120" style="margin-top: 20px"></a>
