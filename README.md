This dataset and analysis code are used in the paper "Does Urban Shrinkage Impact Excess Deaths During the COVID-19 Pandemic?" The description and sources of each file are as follows:

[Step 1: Measurement of Monthly Excess Deaths and the Number of Peaks]

1. Data
	- "total_data.xlsx": This file contains data on urban types and monthly observed, expected, and excess deaths for a total of 1,142 counties.
	- "output-data.zip": This ZIP file contains folders for historical-deaths, expected-deaths-models, and excess-deaths. The historical-deaths folder includes raw monthly mortality data for the 1,142 counties studied. The expected-deaths-models folder contains model data for each county calculated in "code_excess_deaths.R". The excess-deaths folder includes the final results for excess deaths calculated using these models.
2. Code
	- "code_excess_deaths.R": This is the code for calculating excess deaths, which utilizes code from The Economist (https://github.com/TheEconomist/covid-19-excess-deaths-tracker).
	- "peak_detection.ipynb": This Python code is used to calculate the number of peaks.


[Step 3: Data and Statistical Modeling]

1. Data
	- "total_data.xlsx": This file contains data on urban types and monthly observed, expected, and excess deaths for a total of 1,142 counties.
	- "county_data.xlsx": This file contains the data ultimately used for the analysis. It includes the average excess deaths per 100,000 population from March 2020 to February 2023, peak deaths, population, the compound annual growth rate (CAGR) of the population and GRDP, median household income (Income), the percentage of the population over the age of 65 (The elderly), the percentage of adults with a bachelor’s degree or higher from 2017 to 2021 (Education), the unemployment rate (Unemployment), the percentage of white people (Race_white), and the percentage of Black or African American people (Race_black) for each county.
2. Code
	- "code_statistics.ipynb": This file contains the code used for the Kruskal-Wallis test (post hoc) and the mixed effect model in this study.
