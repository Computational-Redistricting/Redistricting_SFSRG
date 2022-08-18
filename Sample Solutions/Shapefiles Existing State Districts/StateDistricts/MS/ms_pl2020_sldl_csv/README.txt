2020 PL 94-171 Data Summary File for Mississippi based on the Decennial Census at the State Legislative District Lower (SLDL) level

Please note that we have NOT validated against the official data used by your state’s redistricting body or bodies. Some states reallocate incarcerated persons and/or exclude non-permanent residents from the PL 94-171 data file for redistricting. Other states may make additional modifications. For more information about state modifications, visit our PL 94-171 Modifications article included in https://redistrictingdatahub.org/data/about-our-data/

Please note that the Census Bureau does not account for changes to congressional, and state legislative district boundaries since the 2018 election. For more information, contact your state's redistricting body or check out All About Redistricting https://redistricting.lls.edu/

##Redistricting Data Hub (RDH) Retrieval Date
08/12/21

##Sources
2020 PL 94-171 Legacy Format Data was retrieved from the US Census FTP site https://www2.census.gov/programs-surveys/decennial/2020/data/01-Redistricting_File--PL_94-171/

##Fields
For a full list of fields and descriptions, review the 2020 Census State Redistricting Data (Public Law 94-171) Summary File Technical Documentation at https://www2.census.gov/programs-surveys/decennial/2020/technical-documentation/complete-tech-docs/summary-file/2020Census_PL94_171Redistricting_StatesTechDoc_English.pdf

##Processing
The legacy format 2020 PL 94-171 Data was downloaded from the Census. 
The legacy format data is provided in one zip file per state. Each zip file contains four files: 3 “segments” containing the data for 1 or more standard redistricting tables, and 1 “geographic header” file. 
The first segment contains the data for Tables P1 (Race) and P2 (Hispanic or Latino, and Not Hispanic or Latino by Race). The second segment contains data for Tables P3 (Race for the Population 18 Years and Over), P4 (Hispanic or Latino, and Not Hispanic or Latino by Race for the Population 18 Years and Over), and H1 (Occupancy Status). The third segment contains Table P5 (Group Quarters Population by Major Group Quarters Type), which was not part of the 2010 PL 94-171 data release.
The files were imported into Python as pipe-delimited data frames and the columns renamed. The segments were joined together and then to the geo file, using the logical record number, or LOGRECNO.
The 10 different summary levels that the RDH processed were queried from the data, each corresponding to a particular unit of geography: state, county, tract, block group, block, congressional district, state legislative district - lower, state legislative district - upper, minor civil division, and census place.
Then the corresponding geographies were merged with the 2020 Census Redistricting Data (P.L. 94-171) shapefiles based on Census GEOIDs. 
Finally, the tabulated data were exported in CSV and shapefile formats.
The RDH verified results internally and externally with partner organizations. 
For additional processing information, review our GitHub https://github.com/nonpartisan-redistricting-datahub

##Additional Notes
For more information about this data set, visit our PL 94-171 article included in https://redistrictingdatahub.org/data/about-our-data/. 

For additional questions, email info@redistrictingdatahub.org.