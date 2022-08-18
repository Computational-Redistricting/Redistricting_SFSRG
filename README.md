# ParallelRedistricting
This repository contains the code for creating solutions to the political redistricting problem using the SFSR-G algorithm described in: Haas, C., Miller, P., Kimbrough, S.O. (2022). An algorithmic approach to legislative apportionment bases and redistricting. Journal of Electoral Studies. Elsevier. 

# README for Redistricting Program

The Redistricting Python program consists of several building blocks:
1) At the core, we have the Seed-Fill-Shift-Repair and Seed-Fill-Shift-Repair-Greedy heuristics by Haas, Miller, and Kimbrough, which generate different heuristic solutions using constraint optimization.
2) Direct integration of the Census API: Currently, the program supports automated loading of the relevant Census tract-level data for different states.
3) A pareto-dominance approach for various evaluation criteria, comparing different solutions across a multitude of objectives and returning/saving the non-dominated solutions.

# Neccessary libraries:
- libpysal
- pysal (install libspatialindex first for Unix systems to prevent potential issues, then try pip install pysal if conda times out)
- simpledbf
- pandas
- geopandas
- cenpy (for the Census API functionality)

# Running the program
To use the program, follow these steps:
# Running the program
To use the program, follow these steps:
1) Specify the extent of the run in the parameter file 'params.json':
	State, State-ID, and state_abr need to refer to the State for which the analysis is run. Examples: "Pennsylvania", 42, "PA" refers to Pennsylvania. 
	These parameters are used to download census tract-level data from the Census API as well as ACS 2018
	Please refer to the StateCodes.txt for the specific IDs
	
	Further parameters: 
	- level: either 'tract' or 'blockgroup': this will calculate the solutions based on either Census tracts or blockgroups. Note: the underlying Census data needs to be available in the Data/State folder
	- use_SFSR: "True" indicates that the SFSR algorithm is run with the specified parameters. If False, the program will simply calculate the evaluation metrics for existing solutions (if any)
	- num_SFSR_runs : how many times should the SFSR algorithm be run, i.e., how many new solutions should be calculated? 
	- num_processes: number of parallel processes to use. 
	- population_tolerance: The population tolerance parameter. Example: 0.01 indicates that each district needs to be within 1% of the ideal population.
	- use_num_districts_from_census: if "True" the number of districts is calculated from the existing shapefiles, i.e., existing number of congressional districts. Otherwise the provided number will be used.
	- num_districts: the number of districts that should be used. Is only used if use_num_districts_from_census equals "False"
	- shapefile_location: if left empty (""), then the corresponding shapefile for the specified analysis level is used (e.g., tracts, blockgroups, etc.). 
	- data_location: if left empty (""), then the corresponding population data for the specified analysis level is used (e.g., tracts, blockgroups, etc.). 
	- manual_save_location: default is "". Use to specify which directory the solutions should be saved in. If Census data is used this is not required and will be set automatically to Solutions/State-Name-Abbreviation/.
	- shapefile_id_column: the column name of the shapefile data that corresponds to the precinct IDs. Necessary to join with demographic data. For Census data, this is "GEOID10". If Census data is used this is not required and will be set automatically.
	- data_id_column: the column name of the demographic data that corresponds to the precinct IDs. Necessary to join with shapefile data. For Census data, this is "GEOID10". If Census data is used this is not required and will be set automatically.
	- population_column: the column name of the demographic data that is used to calculate precinct and district populations. For Census data, this is "total_population". If Census data is used this is not required and will be set automatically.

2)	In the main level directory of the python program, run 'python Main.py'
First, if no manual shapefiles are provided (shapefile_location = "", data_location = ""), it will create a connection to the Census API and download the Census and ACS 2018 data for the specified State.
Note: If the data for the State already exists in the Data/ folder, this will be skipped to avoid unnecessary downloads / load times.
Also, sometimes this can cause an error message on the side of the Census API, e.g., when too many calls are made or when the website does not respond. 
In this case, try running the command again or at a later time, as it should resolve by itself.

After getting / loading the shapefiles, it will run the main program and the SFSR algorithm. 
At the end, the costs of the different plans will be found in "Solutions/State-Abbrevation/costs.csv" or in the manually specified save location (see parameters).

# Folder structure:
The main folders are:
- Main folder: 
	Includes the CensusImport.py file that implements the Census API
	Main.py: the main file to run the SFSR and postprocessing heuristics based on the specified parameters
	
- Data: Subfolders for each State that contain the Census and ACS information. Will be downloaded automatically if the state information is specified correctly in the params.json file and when the CensusImport.py is executed.
	
- Functions: the individual functions for SFSR and postprocessing, implementing parallelized versions of the original files
	
- Solutions: Subfolders for each state: 
	- SFSR contains the SFSR solutions	
	
# Evaluation criteria:
Following functions are currently implemented:
- Contiguity constraint
- Population tolerance constraint
- Number of counties split into more than one district
- Upper and lower population tolerance of solutions
- Number of minority and majority districts for white, non-white, black, and latino populations. Note: in order for this to be calculated, the dataset needs to include columns white_population, black_population, and latino_population in addition to the total_population column (or whatever name is specified in the population-column parameter)


Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg