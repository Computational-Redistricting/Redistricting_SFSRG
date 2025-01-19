"""

This file is a main file to run the Seed-Fill-Shift-Repair-Greedy (SFSR-G) Heuristic from one single function/file

Initially created on Fr 2019-11-15

"""
import pandas as pd
from Functions.sfsr import SFSR
from Functions.EvaluationFunctions import *

import os
import re
import json
import pathlib
import glob
import geopandas as gp
from functools import partial
import multiprocessing
from simpledbf import Dbf5
import argparse
import CensusImport

from Functions.sfsr import prepare_shapefile_data_table

script_dir = os.path.dirname(__file__)  # <-- absolute directory the script is in


def read_parameters():
    '''
    Reads the algorithm parameters from the specified json file
    :return: A set of parameters used to specify the SFSR run
    '''
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
    state_name = parameters['state'].strip()
    print("state name to be used: ", state_name)
    state_name_abbr = parameters['state_abbr'].strip()
    print("state name abbreviation to be used: ", state_name_abbr)
    state_id = parameters['state_id']
    print("state ID to be used: ", state_id)
    use_SFSR = parameters['use_SFSR']
    print("use_SFSR: ", use_SFSR)
    num_sfsr_runs = parameters['num_SFSR_runs']
    print("number of sfsr runs: ", num_sfsr_runs)
    numProcesses = parameters['num_processes']
    print("number of processes: ", numProcesses)
    ensure_contiguity = parameters['ensure_contiguity'].strip()
    print("ensure_contiguity: ", ensure_contiguity)
    ensure_population_tolerance_sfsr = parameters['ensure_population_tolerance_sfsr'].strip()
    print("ensure_population_tolerance_sfsr: ", ensure_population_tolerance_sfsr)
    save_pareto_solutions_in_folder = parameters['save_pareto_solutions_in_folder'].strip()
    print("save_pareto_solutions_in_folder: ", save_pareto_solutions_in_folder)
    population_tolerance = parameters['population_tolerance']
    print("population_tolerance: ", population_tolerance)
    use_num_districts_from_census = parameters['use_num_districts_from_census'].strip()
    print("use_num_districts_from_census: ", use_num_districts_from_census)
    num_districts = parameters['num_districts']
    print("num_districts: ", num_districts)
    use_cvap = parameters['use_cvap'].strip()
    print("use_cvap: ", use_cvap)
    level = parameters['level'].strip()
    print("level: ", level)
    if not (level == 'tract' or level == 'blockgroup'):
        print(f"level = {level} should be 'tract' or 'blockgroup'")
    shapefile_location = parameters['shapefile_location'].strip()
    print("shapefile_location: ", shapefile_location)
    data_location = parameters['data_location'].strip()
    print("data_location: ", data_location)
    shapefile_id_column = parameters['shapefile_id_column'].strip()
    print("shapefile_id_column: ", shapefile_id_column)
    data_id_column = parameters['data_id_column'].strip()
    print("data_id_column: ", data_id_column)
    manual_save_location = parameters['manual_save_location'].strip()
    print("manual_save_location: ", manual_save_location)
    population_column = parameters['population_column'].strip()
    if population_column == '':
        population_column = 'total_population'
    print("population_column: ", population_column)
    use_original_sfsr = parameters['use_original_sfsr'].strip()
    print("use original sfsr: ", use_original_sfsr)
    return state_name, state_name_abbr, state_id, use_SFSR, num_sfsr_runs, numProcesses, \
           ensure_contiguity, ensure_population_tolerance_sfsr, \
           save_pareto_solutions_in_folder, \
           population_tolerance, use_num_districts_from_census, \
           num_districts, use_cvap, level, shapefile_location, data_location, shapefile_id_column, data_id_column, \
           manual_save_location, population_column, use_original_sfsr


def chunks(lst, k):
    """Yield successive k-sized chunks from lst.
    Code taken from https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks"""
    for i in range(0, len(lst), k):
        yield lst[i:i + k]


def read_files_parallel(filename, data_table, shapefile_id_column,
                        number_postprocessing_runs_per_solution):
    """
    Read in files in parallel
    :param data_table: the geo dataframe with population and geometry information
    :param filename: the filename to read
    :param number_postprocessing_runs_per_solution: number of different postprocessing heuristics that were run
    :return: a list of merged plans and their filenames
    """
    if 'runlogs.csv' in filename:
        return
    print('next filename: ', filename)
    df_cols = pd.read_csv(filename, nrows=0)
    cols_used = [shapefile_id_column, 'District']
    # if 'GEOID' in df_cols:
    #     cols_used = ['GEOID', 'District']
    # elif 'GEOID10' in df_cols:
    #     cols_used = ['GEOID10', 'District']
    # else:
    #     ValueError('Either GEOID or GEOID10 column needs to be in the data')

    df = pd.read_csv(filename, usecols=cols_used)
    df.set_index(shapefile_id_column, inplace=True)
    df.index = df.index.astype(str)

    local_data_table = copy.deepcopy(data_table)
    local_data_table.index = local_data_table.index.astype(str)

    # lookup_dict = pd.Series(df.District.values, index=df.index).to_dict()
    # local_data_table['District'] = local_data_table.index.map(lookup_dict)
    local_data_table = local_data_table.merge(df, on=shapefile_id_column)

    temp_plans = []
    temp_filenames = []
    for i in range(number_postprocessing_runs_per_solution):
        if local_data_table is not None:
            temp_plans.append(local_data_table)
        if filename is not None:
            temp_filenames.append([filename, i])

    if temp_plans is not None and temp_filenames is not None:
        return temp_plans, temp_filenames


def getDuplicateColumns(df):
    '''
    This code was taken from https://thispointer.com/how-to-find-drop-duplicate-columns-in-a-dataframe-python-pandas/
    Get a list of duplicate columns.
    It will iterate over all the columns in dataframe and find the columns whose contents are duplicate.
    :param df: Dataframe object
    :return: List of columns whose contents are duplicates.
    '''
    duplicateColumnNames = set()
    # Iterate over all the columns in dataframe
    for x in range(df.shape[1]):
        # Select column at xth index.
        col = df.iloc[:, x]
        # Iterate over all the columns in DataFrame from (x+1)th index till end
        for y in range(x + 1, df.shape[1]):
            # Select column at yth index.
            otherCol = df.iloc[:, y]
            # Check if two columns at x 7 y index are equal
            if col.equals(otherCol):
                duplicateColumnNames.add(df.columns.values[y])
    return list(duplicateColumnNames)


def remove_duplicates_from_path(solution_directory, id_column=None):
    '''
    Basic function to remove duplicate SFSR solutions from the respective folder.
    :param solution_directory: Directory that should be checked for duplicate solutions
    '''
    all_files = glob.glob(solution_directory + "/*.csv")

    # delete the runlog and runtimes files
    all_files = [item for item in all_files if "runlogs.csv" not in item]
    all_files = [item for item in all_files if "Runtimes.csv" not in item]
    all_files = [item for item in all_files if "Plan_Metrics.csv" not in item]
    all_files = [item for item in all_files if "Costs.csv" not in item]

    if len(all_files) == 0:
        return

    filename = all_files[0]
    cols_used = [id_column]

    solutions_df = pd.read_csv(filename, usecols=cols_used)

    # if 'GEOID10' in df_cols:
    #     solutions_df = solutions_df.rename(columns={"GEOID10": "GEOID"})
    solutions_df = solutions_df.rename(columns={id_column: "GEOID"})

    # get initial GEOIDs for the first solution to start the df

    for i, filename in enumerate(all_files):

        df = pd.read_csv(filename, usecols=[id_column,'District'])
        df = df.rename(columns={"District": str("District" + str(i))})
        # if 'GEOID10' in df_cols:
        #     df = df.rename(columns={"GEOID10": "GEOID"})
        df = df.rename(columns={id_column: "GEOID"})

        solutions_df = solutions_df.merge(df, on=('GEOID'))

    # detect potential duplicates
    duplicated_solutions = getDuplicateColumns(solutions_df)
    df = solutions_df.drop(columns=duplicated_solutions)

    # for each duplicated column, delete the corresponding solution from the directory
    for duplicate in duplicated_solutions:
        ix = df.columns.get_loc(duplicate) - 1 # the first index would be the GEOID
        os.remove(all_files[ix])

def calculate_and_save_metrics(chunk_size, congress_districts, path, type, data_table,
                               include_minority_majority, include_num_county_split, level_string, manual_save_location,
                               num_districts, pool, population_column, population_tolerance, state_path, use_cvap,
                               shapefile_id_column):

    existing_cost_filename = None
    # now, calculate the metrics for the county optimized solutions
    cost_columns = ['solution_id', 'compactness',
                    'compactness_radial', 'compactness_polsby_popper', 'compactness_schwartzberg',
                    'compactness_reock', 'compactness_convex_hull_ratio', 'compactness_moment_of_inertia',
                    'population_deviation', 'lower_tolerance', 'upper_tolerance']
    if include_num_county_split == "True":
        cost_columns.append('num_county_split')
    if include_minority_majority == 'True':
        cost_columns.extend(['w75', 'w50', 'b50', 'b40', 'nw50', 'nw40', 'l50', 'l40'])
    cost_df = pd.DataFrame(columns=cost_columns)
    if manual_save_location == "":
        all_files = glob.glob(path + "/*.csv")
    else:
        manual_files = glob.glob(manual_save_location + "/*.csv")
        all_files = manual_files
    # delete the runlog and runtimes files
    all_files = [item for item in all_files if "runlogs.csv" not in item]
    all_files = [item for item in all_files if "Runtimes.csv" not in item]
    all_files = [item for item in all_files if "Plan_Metrics.csv" not in item]
    # # calculate the weight object and re-use it to speed up calculations
    # files = []
    # for file in all_files[0]:
    #     files.extend(read_files_parallel(file, data_table, 1))
    # solution_candidate = files[0][0]
    # weightsDF = ps.lib.weights.Rook.from_dataframe(solution_candidate, ids=solution_candidate.index.tolist())
    # read in the existing data and costs
    if manual_save_location == "":
        cvap_str = "Total_NumDist" + str(num_districts) + "_PopTol" + str(population_tolerance)
        if use_cvap == "True":
            cvap_str = "CVAP_NumDist" + str(num_districts) + "_PopTol" + str(population_tolerance)
        existing_cost_filename = state_path + '' + type + '_Costs_' + level_string + '_' + cvap_str + '.csv'
        print('existing cost file name: ', existing_cost_filename)
        print('does the cost file exist already? ', os.path.exists(existing_cost_filename))
        if os.path.exists(existing_cost_filename):
            cost_df = pd.read_csv(existing_cost_filename, usecols=cost_columns)

            # then, get the names of the solutions which are not already in the csv
            existing_solutions = cost_df['solution_id'].to_list()
            existing_solutions = '\t'.join(existing_solutions)

            all_files_new = []
            for solution in all_files:
                if solution not in existing_solutions:
                    all_files_new.append(solution)

            all_files = all_files_new
    # if we have new solutions, run the subsequent evaluation functions to calculate various metrics
    if len(all_files) != 0:

        chunk_size = min(len(all_files), chunk_size)

        all_files_chunks = chunks(all_files, chunk_size)

        for files in all_files_chunks:

            # read in all plans in parallel
            results = pool.map(partial(read_files_parallel, data_table=data_table,
                                       shapefile_id_column = shapefile_id_column,
                                       number_postprocessing_runs_per_solution=1), files)
            plans = [i[0][0] for i in results]
            plan_names = [i[1][0] for i in results]

            # calculate the costs in parallel
            costs = pool.map(partial(calculate_metrics,
                                     congress_districts=congress_districts, population_column=population_column,
                                     include_num_county_split=include_num_county_split,
                                     include_minority_majority=include_minority_majority), plans)

            for i, solution in enumerate(plans):
                # print("plan name: ", plan_names[i][0])
                next_row = [plan_names[i][0]]

                next_row.extend(costs[i])

                print('next row', next_row)
                cost_df.loc[len(cost_df)] = next_row

        # delete potential duplicates
        cost_df.drop_duplicates()
        if manual_save_location == "":
            cost_df.to_csv(existing_cost_filename)
        else:
            cost_df.to_csv((manual_save_location + 'Costs.csv'))


def main_run(pool):
    '''
    Main function to run SFSR and potential post-processing heuristics
    First, shapefiles are downloaded via the Census API
    Then, based on provided parameters, SFSR iterations and postprocessing heuristics are run.
    :param pool: A multiprocessing pool to parallelize the calculations.
    '''

    print('## Load parameter values from param.json ##')
    state_name, state_name_abbr, state_id, use_SFSR, num_sfsr_runs, numProcesses, \
        ensure_contiguity, ensure_population_tolerance_sfsr, \
        save_pareto_solutions_in_folder, population_tolerance, \
        use_num_districts_from_census, num_districts, use_cvap, level, shapefile_location, data_location,\
        shapefile_id_column, data_id_column, manual_save_location, population_column, use_original_sfsr\
        = read_parameters()

    print('If necessary, create the folder structure')
    state_path = 'Solutions/' + state_name_abbr + '/'

    # specify the number of districts that should be used
    if use_num_districts_from_census == "True" and os.path.isfile('Data/' + state_name_abbr +
                                                                  '/Congressional_Districts/Congressional_Districts.dbf'):
        num_districts_file = Dbf5('Data/' + state_name_abbr + '/Congressional_Districts/Congressional_Districts.dbf')

        num_districts_df = num_districts_file.to_dataframe()

        numDistricts = num_districts_df.shape[0]  # Looking to the 2010 census
    elif num_districts == '':
        ValueError("Need to specify the number of districts in the parameter file.")
    else:
        numDistricts = num_districts

    level_string = "Tracts"
    if level == 'blockgroup':
        level_string = "Blockgroups"
    elif level == 'block':
        level_string = 'Blocks'

    if shapefile_location == "":
        print('## Download Census Data if necessary ##')
        CensusImport.main(state_id=state_id, state_name=state_name, state_name_abbr=state_name_abbr)

        shapefile_location = 'Data/' + state_name_abbr + '/Census_' + level_string+'_2018/Census_' + \
                             level_string+'_2018.shp'
        if level =='block':
            shapefile_location = 'Data/' + state_name_abbr + '/Census_' + level_string + '_2010/Census_' + \
                                 level_string + '_2010.shp'

    cvap_str = "Total/NumDist" + str(num_districts) + "/PopTol"+str(population_tolerance)
    population = 'total'
    if use_cvap == "True":
        population = 'cvap'
        cvap_str = "CVAP/NumDist" + str(num_districts) + "/PopTol"+str(population_tolerance)
        if data_location == "":
            data_location = 'Data/' + state_name_abbr + '/Census_'+level_string+'_Demographics_CVAP.csv'
    else:
        if data_location == "":
            data_location = 'Data/' + state_name_abbr + '/Census_'+level_string+'_Demographics.csv'
    if manual_save_location == "":
        sfsr_path = 'Solutions/' + state_name_abbr + '/SFSR/'+level_string+'/' + cvap_str
    else:
        sfsr_path = manual_save_location
    pareto_solutions_path = 'Solutions/' + state_name_abbr + '/Pareto/'+level_string+'/' + cvap_str

    congress_districts_location = 'Data/' + state_name_abbr + '/Congressional_Districts/Congressional_Districts.shp'
    if os.path.isfile(congress_districts_location):
        congress_districts = gp.read_file(congress_districts_location, engine="pyogrio")
        # project it to the 3857 (Mercator) CRS for more accurate boundary calculations
        congress_districts = congress_districts.to_crs(epsg=3857)
    else:
        congress_districts = None
    pathlib.Path(sfsr_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(pareto_solutions_path).mkdir(parents=True, exist_ok=True)

    # due to the large number of files and potential memory restrictions, we subset the files into chunks
    chunk_size = 100

    if shapefile_id_column is None or shapefile_id_column == "":
        shapefile_id_column = 'GEOID10'
    if data_id_column is None or data_id_column == "":
        data_id_column = 'GEOID10'

    # Step 0
    # read in basic geometry for the current State
    data_table = prepare_shapefile_data_table(shapefile_location=shapefile_location,
                                              data_location=data_location,
                                              shapefile_id_column=shapefile_id_column,
                                              data_id_column=data_id_column,
                                              population_column=population_column)

    # Step 1: Start SFSR
    print('## Start SFSR ##')
    if use_SFSR == "True":
        # start the SFSR code

         # # if we look at block level data we also need to provide blockgroup data
        shapefile_location_blockgroups = None
        data_location_blockgroups = None
        if level == 'block':
            shapefile_location_blockgroups = 'Data/' + state_name_abbr + \
                                             '/Census_Blockgroups_2018/Census_Blockgroups_2018.shp'
            data_location_blockgroups = 'Data/' + state_name_abbr + '/Census_Blockgroups_Demographics.csv'

        # then, start the main SFSR function
        sfsr = SFSR(numRuns=num_sfsr_runs, pathToData='Data/' + state_name_abbr + '/', nDistricts=numDistricts,
                    tolerance=population_tolerance,
                    neighborsType='rook', pathToOutputData=sfsr_path,
                    cvap=use_cvap, level=level, data_table=data_table,
                    shapefile_id_column=shapefile_id_column,
                    data_id_column=data_id_column,
                    population_column=population_column, use_original_sfsr=use_original_sfsr)

        sfsr.basicgo_parallel(num_runs=num_sfsr_runs, pool=pool)

    else:
        print("Skip SFSR")

    print("## Redistricting completed ##")

    # Step 2: Drop potential duplicate solutions
    if manual_save_location == "":
        remove_duplicates_from_path(sfsr_path, shapefile_id_column)
    else:
        remove_duplicates_from_path(manual_save_location, shapefile_id_column)

    # Step 3: Calculate the solution properties
    print('read in previously generated solutions')

    # here, calculate the population-based statistics for both Total and CVAP population if available.
    # otherwise, if CVAP is not available only calculate the regular population-based statistics

    try:
        path = sfsr_path
        type = "SFSR"
        # note: only include the following two aspects if minority/majority data is available and
        # if there is a county column in the dataset
        include_minority_majority = "False"
        include_num_county_split = "False"

        calculate_and_save_metrics(chunk_size, congress_districts, path, type, data_table,
                                   include_minority_majority, include_num_county_split, level_string,
                                   manual_save_location, numDistricts, pool, population_column, population_tolerance,
                                   state_path, use_cvap, shapefile_id_column)

    except Exception as e:
        print('Exception during calculating solution costs/characteristics')
        print(e)


if __name__ == '__main__':
    print('Starting the main function')

    parser = argparse.ArgumentParser(description='Optional arguments for running python script')

    # Optional argument
    parser.add_argument('--parameter_file', type=str,
                        help='specify the name of the parameter file to be used')
    # parse the arguments
    args = parser.parse_args()
    if args.parameter_file is not None:
        file_name = args.parameter_file
    else:
        file_name = 'params.json'
    parameters = json.load(open(file_name, 'rb'))
    numProcesses = parameters['num_processes']
    print("number of processes: ", numProcesses)

    # set the number of SFSR runs and the number of cores/processes to use for parallelization
    maxProcesses = multiprocessing.cpu_count() # this is the maximum number of cores in your system
    numProcesses = min(numProcesses, maxProcesses) # make sure you're using not more than the available number of cores

    print("number processes: ", numProcesses)
    pool = multiprocessing.Pool(processes=numProcesses)
    print('created multiprocessing pool')

    main_run(pool=pool)

