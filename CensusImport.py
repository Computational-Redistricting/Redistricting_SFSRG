# import the census api communication packages
import cenpy as cp
import pandas as pd
# from cenpy import products
from pathlib import Path
import geopandas as gp
import pathlib
import json
import os
import argparse
from Functions.products import ACS, Decennial2010
import zipfile


parser = argparse.ArgumentParser(description='Optional arguments for running python script')

# Optional argument
parser.add_argument('--parameter_file', type=str,
                    help='specify the name of the parameter file to be used')
# parse the arguments
args = parser.parse_args()

# use them in the global context
if args.parameter_file is not None:
    file_name = args.parameter_file
else:
    file_name = 'params.json'
parameters = json.load(open(file_name, 'rb'))

default_state_name = parameters['state'].strip()
default_state_name_abbr = parameters['state_abbr'].strip()
default_state_id = parameters['state_id']
default_level = 'tract'
# default_level = 'blockgroup'

acs_year = 2018

cvap_location = 'Data/CVAP_2014-2018_ACS_csv_files.zip'

def main(state_id=default_state_id, state_name=default_state_name,state_name_abbr=default_state_name_abbr, level=default_level):
    print('starting the main function')

    # check if the directory exists, otherwise create it for the given state
    # note: if the directory already exists, this should do nothing
    path = 'Data/' + state_name_abbr + "/"
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    path_cd = path + 'Congressional_Districts'

    if not os.path.exists(path_cd):
        # first, get the congressional districts from the 2010 census for this state
        conn = cp.remote.APIConnection('DECENNIALSF12010')

        # layer 62 is for the voting districts
        conn.set_mapservice('tigerWMS_Census2010')
        state_specifier = 'STATE = ' + str(state_id)
        cd = conn.mapservice.query(layer=62, where=state_specifier)

        cd = cd.to_crs({'init': 'epsg:4326'})
        cd.to_file(path_cd, index='False')


    # then, get the census tract information along with some ACS information for population and income
    acs_file_tracts = Path("Data/" + state_name_abbr + '/' + "Census_Tracts_" + str(2018))
    acs_file_blockgroups = Path("Data/" + state_name_abbr + '/' + "Census_Blockgroups_" + str(2018))

    print('only get the data if it does not already exist')

    # if not acs_file.is_file():
    if not os.path.exists(acs_file_tracts):
        # acs = products.ACS(2018)
        acs = ACS(acs_year)
        # get the acs and census data for the specified state from ACS and Census
        # for the ACS data, we restrict ourselves to certain population columns for now (can add income as well)

        # first, let's get the tracts. afterwards, we can get the blockgroups
        acs_state = acs.from_state(state_name, level='tract',
                                   variables=['B03002_001E', 'B03002_003E', 'B03002_004E', 'B03002_012E',
                                              'B07011_001E',
                                              'B29001_001E', 'B29004_001E'])
        # census_state = census.from_state(state, level='tract')

        # do some conversions to get the final table
        acs_state.rename(columns={'B03002_001E': 'total_population', 'B03002_003E': 'white_population',
                                  'B03002_004E': 'black_population', 'B03002_012E': 'latino_population',
                                  'B07011_001E': 'median_income', 'GEOID': 'GEOID10', 'NAME': 'geo_name',
                                  'B29001_001E': 'total_population_voting_age',
                                  'B29004_001E': 'median_income_voting_age'},
                         inplace=True)

        # calculate the other population based on existing values
        acs_state['other_population'] = acs_state['total_population'] - acs_state['white_population'] - acs_state[
            'black_population'] - acs_state['latino_population']

        # acs_state.to_csv(acs_file, index=False)
        # for index, row in acs_state.iterrows():
        #     geom = row['geometry']
        #     if len(geom.coords) <= 2:
        #         print("This row ", index, "has an invalid polygon geometry")
        #         # this is just one example of invalid geometries, there are also overlapping vertices, .
        # acs = gp.GeoDataFrame(acs)

        # acs_state['GEOID10'] = acs_state['geo_fips']

        acs_state = acs_state[["GEOID10", "geo_name", "geometry", "state", "county", "tract", "total_population",
                               "white_population", "black_population",
                               "latino_population", "other_population", "median_income", "total_population_voting_age", "median_income_voting_age"]]

        acs_shape = gp.GeoDataFrame(acs_state[["GEOID10", "geo_name", "geometry", "state", "county", "tract"]])

        # set the crs
        acs_shape.crs = {'init': 'epsg:3857'}
        acs_shape = acs_shape.to_crs({'init': 'epsg:4326'})
        acs_shape.to_file(acs_file_tracts)

        # also save the core information in another dataset
        acs_demographic = acs_state[['GEOID10','total_population','white_population','black_population','other_population','latino_population', 'median_income', 'total_population_voting_age', 'median_income_voting_age']]
        # acs_demographic = acs_state[['GEOID10','total_population','white_population','black_population','other_population','latino_population', 'median_income']]
        acs_demographic.rename(columns={'GEOID10': 'geo_fips'})
        acs_demographic_path = Path("Data/" + state_name_abbr + '/' +'Census_Tracts_Demographics.csv')
        acs_demographic.to_csv(acs_demographic_path)

        zf = zipfile.ZipFile(cvap_location)
        df_tracts = pd.read_csv(zf.open('Tract.csv'), encoding='latin1')
        # print(df_tracts.head())

        df_tracts_pivot = df_tracts.pivot(index='geoname', columns='lntitle', values='cvap_est')

        df_tracts_pivot = df_tracts_pivot.rename(columns={'Total': 'total_population', 'White Alone': 'white_population',
                                                          'Black or African American Alone': 'black_population',
                                                          'Hispanic or Latino': 'latino_population'})
        df_tracts_pivot['other_population'] = df_tracts_pivot['total_population'] - df_tracts_pivot['white_population'] - \
                                              df_tracts_pivot['black_population'] - df_tracts_pivot['latino_population']

        df_tracts_pivot = df_tracts_pivot.reset_index()

        df_tracts_pivot = df_tracts_pivot[
            ['geoname', 'total_population', 'white_population', 'black_population', 'latino_population',
             'other_population']]
        cvap_df = acs_state[['GEOID10','geo_name', 'median_income', 'median_income_voting_age', 'total_population_voting_age']].merge(df_tracts_pivot, left_on='geo_name', right_on="geoname")
        cvap_path = 'Data/' + state_name_abbr + '/Census_Tracts_Demographics_CVAP.csv'
        cvap_df.to_csv(cvap_path, index=False)

    # then, get the blockgroups
    if not os.path.exists(acs_file_blockgroups):
        # acs = products.ACS(2018)
        acs = ACS(acs_year)
        acs_state = acs.from_state(state_name, level='blockgroup',
                                   variables=['B03002_001E', 'B03002_003E', 'B03002_004E', 'B03002_012E',
                                              'B07011_001E',
                                              'B29001_001E', 'B29004_001E'])
        # census_state = census.from_state(state, level='tract')

        # do some conversions to get the final table
        acs_state.rename(columns={'B03002_001E': 'total_population', 'B03002_003E': 'white_population',
                                  'B03002_004E': 'black_population', 'B03002_012E': 'latino_population',
                                  'B07011_001E': 'median_income', 'GEOID': 'GEOID10', 'NAME': 'geo_name',
                                  'B29001_001E': 'total_population_voting_age',
                                  'B29004_001E':'median_income_voting_age'},
                         inplace=True)

        # calculate the other population based on existing values
        acs_state['other_population'] = acs_state['total_population'] - acs_state['white_population'] - acs_state[
            'black_population'] - acs_state['latino_population']

        # acs_state.to_csv(acs_file, index=False)
        # for index, row in acs_state.iterrows():
        #     geom = row['geometry']
        #     if len(geom.coords) <= 2:
        #         print("This row ", index, "has an invalid polygon geometry")
        #         # this is just one example of invalid geometries, there are also overlapping vertices, .
        # acs = gp.GeoDataFrame(acs)

        # acs_state['GEOID10'] = acs_state['geo_fips']

        acs_state = acs_state[["GEOID10", "geo_name", "geometry", "state", "county", "tract", "total_population",
                               "white_population", "black_population",
                               "latino_population", "other_population", "median_income", "total_population_voting_age", "median_income_voting_age"]]

        acs_shape = gp.GeoDataFrame(acs_state[["GEOID10", "geo_name", "geometry", "state", "county", "tract"]])

        # set the crs
        acs_shape.crs = {'init': 'epsg:3857'}
        acs_shape = acs_shape.to_crs({'init': 'epsg:4326'})
        acs_shape.to_file(acs_file_blockgroups)

        # also save the core information in another dataset
        acs_demographic = acs_state[['GEOID10','total_population','white_population','black_population','other_population','latino_population', 'median_income', 'total_population_voting_age', 'median_income_voting_age']]
        acs_demographic.rename(columns={'GEOID10': 'geo_fips'})
        acs_demographic_path = Path("Data/" + state_name_abbr + '/' +'Census_Blockgroups_Demographics.csv')
        acs_demographic.to_csv(acs_demographic_path)

        zf = zipfile.ZipFile(cvap_location)
        df_blockgroups = pd.read_csv(zf.open('BlockGr.csv'), encoding='latin1')
        # print(df_tracts.head())

        df_blockgroups_pivot = df_blockgroups.pivot(index='geoname', columns='lntitle', values='cvap_est')

        df_blockgroups_pivot = df_blockgroups_pivot.rename(columns={'Total': 'total_population', 'White Alone': 'white_population',
                                                          'Black or African American Alone': 'black_population',
                                                          'Hispanic or Latino': 'latino_population'})
        df_blockgroups_pivot['other_population'] = df_blockgroups_pivot['total_population'] - df_blockgroups_pivot['white_population'] - \
                                              df_blockgroups_pivot['black_population'] - df_blockgroups_pivot['latino_population']

        df_blockgroups_pivot = df_blockgroups_pivot.reset_index()

        df_blockgroups_pivot = df_blockgroups_pivot[
            ['geoname', 'total_population', 'white_population', 'black_population', 'latino_population',
             'other_population']]
        cvap_df = acs_state[['GEOID10','geo_name', 'median_income', 'median_income_voting_age', 'total_population_voting_age']].merge(df_blockgroups_pivot, left_on='geo_name', right_on="geoname")
        cvap_path = 'Data/' + state_name_abbr + '/Census_Blockgroups_Demographics_CVAP.csv'
        cvap_df.to_csv(cvap_path, index=False)




if __name__ == '__main__':
    print('Starting the main function of Census Import')
    main(state_id=default_state_id, state_name=default_state_name, state_name_abbr=default_state_name_abbr, level=default_level)
