##This dataset and analysis code are used in the paper "Does Urban Shrinkage Impact Excess Deaths During the COVID-19 Pandemic?" The description and sources of each file are as follows:

#1. Data
	- "total_data: This file contains data on urban types, and monthly observed, expected, and excess deaths for a total of 1,142 counties.
	- "country_data": This file contains the data ultimately used for the analysis. It includes the average excess deaths per 100,000 population from March 2020 to February 2023, peak deaths, population, the compound annual growth rate (CAGR) of the population and GRDP, median household income (Income), the percentage of the population over the age of 65 (The elderly), the percentage of adults with a bachelor’s degree or higher from 2017 to 2021 (Education), the unemployment rate (Unemployment), the percentage of white people (Race_white), and the percentage of black or African American people (Race_black) for each county.

2. Analysis Code
	- "code_excess_deaths.R": This is the code for calculating excess deaths, which utilizes code from The Economist (https://github.com/TheEconomist/covid-19-excess-deaths-tracker).
	- "code_statistics.ipynb": This contains the code for the Kruskal-Wallis test (post hoc) and mixed effect model used in this study.
