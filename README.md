![image](https://github.com/cjabradshaw/SharkManagementStrategiesEffectSize/assets/26937238/8078973f-5873-415b-8e83-80ae782c00e8)# Shark-management strategies effect-size analyses
<a href="https://doi.org/10.5281/zenodo.10043996"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10043996.svg" alt="DOI"></a>
<img align="right" src="www/sharkbitenet.png" alt="shark bite icon" width="280" style="margin-top: 20px">

Analyses to assess probability of detecting differences in the number of shark bites on humans with/without shark management strategies in place.
<a href="https://github.com/cjabradshaw/AustralianSharkIncidentDatabase"><img align="left" src="www/ASIDlogo3.png" alt="ASID logo" width="150" style="margin-top: 20px"></a>

See also the <a href="https://github.com/cjabradshaw/AustralianSharkIncidentDatabase">Australian Shark-Incident Database</a> (ASID).

<br>
<br>
<br>
Prof <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
July 2023/updated October 2023 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>
contributors: <a href="https://www.flinders.edu.au/people/charlie.huveneers">Charlie Huveneers</a>

## Paper
Huveneers, C, C Blount, CJA Bradshaw, PA Butcher, MP Lincoln Smith, WG MacBeth, DP McPhee, N Moltschaniwskyj, VM Peddemors, M Green. 2023. Shifts in the incidence of shark bites and efficacy of beach-focussed mitigation in Australia. <em>Marine Pollution Bulletin</em> 198:115855.
doi:<a href="https://doi.org/10.1016/j.marpolbul.2023.115855">10.1016/j.marpolbul.2023.115855</a>

### Abstract
Shark-human interactions are some of the most pervasive human-wildlife conflicts, and their frequencies are increasing globally. New South Wales (Australia) was the first to implement a broad-scale program of shark-bite mitigation in 1937 using shark nets, which expanded in the late 2010s to include non-lethal measures. Using 196 unprovoked shark-human interactions recorded in New South Wales since 1900, we show that bites shifted from being predominantly on swimmers to 79 % on surfers by the 1980s and increased 2–4-fold. We could not detect differences in the interaction rate at netted versus non-netted beaches since the 2000s, partly because of low incidence and high variance. Although shark-human interactions continued to occur at beaches with tagged-shark listening stations, there were no interactions while SMART drumlines and/or drones were deployed. Our effectsize analyses show that a small increase in the difference between mitigated and non-mitigated beaches could indicate reductions in shark-human interactions. Area-based protection alone is insufficient to reduce shark-human interactions, so we propose a new, globally transferable approach to minimise risk of shark bite more effectively.

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
