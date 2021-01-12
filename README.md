# PSA

This spatial analysis is part of an assignment for UCL GEOG0114 Principles of Spatial Analysis.
This project explores the spatial distribution of burglaries across Manchester (UK) by testing for spatial autocorrelation, and attempts to explain the prevalence of burglaries by modelling it in relation to the Index of Multiple Deprivation (IMD) scores. 

Data:
1) Geographic Boundaries 
Source: https://borders.ukdataservice.ac.uk/ 
Boundary Data Selector: England; Statistical Building Block; 2011 and later;  English Lower Layer Super Output Areas, 2011; Manchester.

2) Crime Data
Source:  https://data.police.uk/data/
Custom Download: January 2018 - December 2018; Greater Manchester Police (GMP); include crime data 
(I used 2018 data because it's latest year for which GMP data was complete for a the whole year)
	
3) English indices of deprivation 2019
Source: GOV.UK https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019
I downloaded File 5, which includes the Index of Multiple Deprivation (IMD) scores, deciles, and national ranking of LSOAs

All the cleaning and pre-processing was done in RStudio, and is included in the .R file
