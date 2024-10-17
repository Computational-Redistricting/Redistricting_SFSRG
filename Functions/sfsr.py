#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the new class for the SFSR algorithm exlusively using the pysal library to calculate and use neighbors
"""
import random
import pandas as pd
import numpy as np
import copy
import os
import datetime
import re
import libpysal as lp
import warnings
from functools import partial
import math
import geopandas as gp
from scipy.sparse.csgraph import connected_components
import time
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from Functions.EvaluationFunctions import calculate_population_equality, calculate_population_tolerance

time_seed = 0
time_fill = 0
time_shift = 0
time_repair = 0

count_seed = 0
count_fill = 0
count_shift = 0
count_repair = 0


def integrate_disconnected_components(data_table, population_column="total_population"):
    '''
    This function identifies disconnected components and island, if they exist
    It then calculates the pairwise distances between the centroid(s) of these units and the main component (mainland)
    units/precincts
    Finally, it adds the population, geometry, and other information to the closest precinct
    :param data_table: The geodataframe, potentially including islands and disconnected components
    :param population_column: the name of the main population column to integrate
    :return: a geo dataframe without disconnected components or islands
    '''
    weightsDF = lp.weights.Rook.from_dataframe(data_table, ids=data_table.index.tolist())
    u, count = np.unique(weightsDF.component_labels, return_counts=True)

    count_sort_ind = np.argsort(-count)
    # now, the first count_sort_ind has the component with the largest number of units. we can use the rest to be added
    # to the repair list
    toReturn = []
    largest_component_label = u[count_sort_ind][0]
    smaller_component_index = np.where(weightsDF.component_labels != largest_component_label)
    for i, ix in enumerate(smaller_component_index[0]):
        # print((aDistrictDF.index[i], aDistrictDF.iloc[ix]['GEOID']))
        toReturn.append(weightsDF.id_order[ix])
    data_table_large_component = data_table.loc[~data_table.index.isin(toReturn)]
    disconnected_components = []
    for ix in toReturn:
        # print(ix)
        # print(data_table.loc[ix,'geometry'].centroid)
        # print(data_table_large_component['geometry'])
        local_centroid = data_table.loc[ix, 'geometry'].centroid
        distances = [local_centroid.distance(centroid2) for centroid2 in
                     data_table_large_component['geometry'].centroid]
        min_distance_index = distances.index(min(distances))
        min_distance_index = data_table_large_component.index.values[min_distance_index]

        # add this information to save which GEOID was added to a GEOID in the connected component
        disconnected_components.append([ix, min_distance_index])

        # print(min_distance_index)
        data_table_large_component.loc[min_distance_index]['geometry'] = \
            data_table_large_component.loc[min_distance_index, 'geometry'].union(
            data_table.loc[ix, 'geometry'])

        data_table_large_component.loc[min_distance_index, population_column] = \
            data_table_large_component.loc[min_distance_index, population_column] + data_table.loc[
            ix, population_column]

        if ("total_population" in data_table.columns) and (population_column != "total_population"):
            data_table_large_component.loc[min_distance_index, 'total_population'] = \
                data_table_large_component.loc[min_distance_index, 'total_population'] + data_table.loc[
                ix, 'total_population']

        if "total_population_voting_age" in data_table.columns:
            data_table_large_component.loc[min_distance_index, 'total_population_voting_age'] = \
                data_table_large_component.loc[min_distance_index, 'total_population_voting_age'] + \
                         data_table.loc[ix, 'total_population_voting_age']

        if "median_income" in data_table.columns:
            data_table_large_component.loc[min_distance_index, 'median_income'] = \
                data_table_large_component.loc[min_distance_index, 'median_income'] + data_table.loc[
                ix, 'median_income']

        if "median_income_voting_age" in data_table.columns:
            data_table_large_component.loc[min_distance_index, 'median_income_voting_age'] = \
                data_table_large_component.loc[min_distance_index, 'median_income_voting_age'] + \
                         data_table.loc[ix, 'median_income_voting_age']

        if "white_population" in data_table.columns:
            data_table_large_component.loc[min_distance_index, 'white_population'] = \
                data_table_large_component.loc[min_distance_index, 'white_population'] + data_table.loc[
                ix, 'white_population']

        if "black_population" in data_table.columns:
            data_table_large_component.loc[min_distance_index, 'black_population'] = \
                data_table_large_component.loc[min_distance_index, 'black_population'] + data_table.loc[
                ix, 'black_population']

        if "latino_population" in data_table.columns:
            data_table_large_component.loc[min_distance_index, 'latino_population'] = \
                data_table_large_component.loc[min_distance_index, 'latino_population'] + data_table.loc[
                ix, 'latino_population']

        if "other_population" in data_table.columns:
            data_table_large_component.loc[min_distance_index, 'other_population'] = \
                data_table_large_component.loc[min_distance_index, 'other_population'] + data_table.loc[
                ix, 'other_population']

    # convert disconnected component data to data frame
    disconnected_components_df = pd.DataFrame(disconnected_components, columns=['Initial_GEOID', 'Merged_GEOID'])

    return data_table_large_component, disconnected_components_df


def prepare_shapefile_data_table(shapefile_location, data_location, shapefile_id_column=None,
                                 data_id_column=None, population_column='total_population'):
    '''
    This function joins the shapefile with the demographic data based on standard or custom id columns
    :param shapefile_location: the location of the shapefile
    :param data_location: the location of the demographic data file
    :param shapefile_id_column: the name of the column with the precinct ids. if none, GEOID10 is used
    :param data_id_column: the name of the column with the precinct ids. if none, GEOID10 is used
    :param population_column: the name of the population column in case population data need to be aggregated
    :return: a combined GeoDataFrame based on the shapefile and demographic data file
    '''
    # read in the shapefile

    data_table = None

    try:
        data_table = gp.read_file(shapefile_location, encoding='utf-8')

        # project it to the 3857 (Mercator) CRS for more accurate boundary calculations
        data_table = data_table.to_crs(epsg=3857)

        # new: convert to multipolygon for easier calculation of compactness measures
        multipolygons = []
        for geo in data_table.geometry:
            if isinstance(geo, Polygon):
                multipolygons.append(MultiPolygon([geo]))
            elif isinstance(geo, MultiPolygon):
                multipolygons.append(geo)
            else:
                ValueError('geometry should be either Polygon or Multipolygon')

        #         print(multipolygons)
        data_table['geometry'] = gp.GeoSeries(multipolygons, index=data_table.geometry.index)

        # merge with additional census data
        data_table[shapefile_id_column] = data_table[shapefile_id_column].astype(str)

        data_table.sort_values(by=shapefile_id_column)
    except FileNotFoundError as fe:
        print(fe)
        print('Shapefile not found at location: ', shapefile_location)
    except KeyError as ke:
        print(ke)
        print('column ', shapefile_id_column, ' not found in shapefile')

    try:
        if data_location is not None and data_location != "":
            census = pd.read_csv(data_location)
            census[data_id_column] = census[data_id_column].astype(str)

            if shapefile_id_column not in census.columns:
                census[shapefile_id_column] = census[data_id_column]

            # identify common columns to avoid column duplication
            common_col_names = np.intersect1d(data_table.columns, census.columns).tolist()

            # make sure the columns have the same type
            for name in common_col_names:
                data_table[name] = data_table[name].astype(census[name].dtype)

            # data_table = pd.merge(data_table, census, left_on=left_names, right_on=right_names)
            data_table = pd.merge(data_table, census, left_on=common_col_names, right_on=common_col_names)
            data_table.set_index(shapefile_id_column, inplace=True)
    except FileNotFoundError as fe:
        print(fe)
        print('Demographic data cannot be read from location: ', data_location)

    data_table, disconnected_components = integrate_disconnected_components(data_table,
                                                                            population_column=population_column)

    return data_table


def first_n_digits(num, n):
    """
    This function simply returns the first n digits from a given number
    :param num: the number from which digits need to be retrieved
    :param n: how many digits should be retrieved
    :return: the first n digits of num
    """
    return num // 10 ** (int(math.log(num, 10)) - n + 1)


def getDistrictPopulation(theDistrictNum,assignDF, population_column='total_population'):
    '''
    Function to return the popualtion for a given district number
    :param theDistrictNum: the district number
    :param assignDF: the current GeoDataFrame with assignments of precincts to districts
    :param population_column: the population column name
    :return: the population for the specified district
    '''
    districtRecords = assignDF[assignDF['District'] == theDistrictNum]
    districtPop = districtRecords[population_column].sum()
    return districtPop


def getDistrictPops(assignDF,districtList, population_column='total_population'):
    '''
    Returns the population for all districts as a list
    :param assignDF: the current GeoDataFrame with assignments of precints to districts
    :param districtList: the list of distinct district IDs
    :param population_column: the population column name
    :return: a list of district populations
    '''
    temp_list = assignDF[assignDF['District'].isin(districtList)]
    thePops = temp_list[['District',population_column]].groupby('District')[population_column].sum().tolist()
    return thePops


def notTooBig(populationList,idealPop,topMaxPerCent):
    '''
    Function to check if the district with maximum population is outside the allowed population range
    :param populationList: the list of district populations
    :param idealPop: the ideal population of a district (overall population divided by number of districts)
    :param topMaxPerCent: the maximum deviation from the ideal population
    :return: True if all district populations are smaller than the maximum positive deviation from the ideal population
    '''
    theMax = max(populationList)
    if theMax < idealPop * (1 + topMaxPerCent):
        return True
    else:
        return False


def notTooSmall(populationList,idealPop,bottomMinPerCent):
    '''
    Function to check if the district with minimum population is outside the allowed population range
    :param populationList: the list of district populations
    :param idealPop: the ideal population of a district (overall population divided by number of districts)
    :param bottomMinPerCent: the minimum (bottom) deviation from the ideal population
    :return: True if all district populations are larger than the minimum deviation from the ideal population
    '''
    theMin = min(populationList)
    if theMin > idealPop * (1 - bottomMinPerCent):
        return True
    else:
        return False


def popSizeOK(anAssignmentDF,idealPop,topMaxPerCent,bottomMinPerCent,districtList,population_column):
    '''
    Function to check if the district populations are in the allowed population range
    :param anAssignmentDF: the current precinct-district assignment GeoDataFrame
    :param idealPop: the ideal population of a disrict
    :param topMaxPerCent: the maximum positive deviation from the ideal population
    :param bottomMinPerCent: the maximum negative deviation from the ideal population
    :param districtList: the list of districts
    :param population_column: the name of the column with the population information per precinct
    :return: True if all district populations are within the allowed population range
    '''
    populationList = getDistrictPops(anAssignmentDF,districtList, population_column)
    if notTooBig(populationList,idealPop,topMaxPerCent) and \
    notTooSmall(populationList,idealPop,bottomMinPerCent):
        return True
    else:
        return False


def getGroupPopulation(theDistrictNum, assignDF, population_column='total_population'):
    '''
    Gets the population of a specific district
    :param theDistrictNum:
    :param assignDF:
    :param population_column:
    :return:
    '''
    districtRecords = assignDF[assignDF['Group'] == theDistrictNum]
    districtPop = districtRecords[population_column].sum()
    return districtPop


def getGroupPops(assignDF, districtList, population_column='total_population'):
    '''
    Returns the list of populations for a list of districts
    :param assignDF:
    :param districtList:
    :param population_column:
    :return:
    '''
    temp_list = assignDF[assignDF['Group'].isin(districtList)]
    # print(assignDF.head)
    thePops = temp_list[['Group', population_column]].groupby('Group')[population_column].sum().to_dict()
    return thePops


def getDistrictPopsDict(assignDF, districtList, population_column='total_population'):
    '''
    Returns a dictionary of populations for a specified list of districts.
    :param assignDF:
    :param districtList:
    :param population_column:
    :return:
    '''
    temp_list = assignDF[assignDF['District'].isin(districtList)]
    thePops = temp_list[['District', population_column]].groupby('District')[population_column].sum().to_dict()
    return thePops


def fastContiguityCheckForChangedDistrict(weights, listOfPrecincts):
    '''
    This is a fast method of checking if a list of precincts is contiguous by using the sparse matrix from weights
    The list of precincts needs to be a subset of the id's of the weights object
    :param weights:
    :param listOfPrecincts:
    :return:
    '''
    ids = weights.id_order
    idxs = np.searchsorted(ids, listOfPrecincts)

    num_components = connected_components(weights.sparse[idxs][:, idxs], directed=False)[0]
    if num_components > 1:
        return False
    else:
        return True


def seed_fill_single_blockgroup(initialDF, numberSubgroupsPerBlockgroup):
    '''
    Runs the initial seed and fill step of the SFSR procedure
    First, random precincts are chosen and seeded with different district IDs
    Then, Fill assigns previously unassigned precincts to an existing district
    Notes:
        - the initial seed and fill routine creates a contiguous assignment where each precinct belongs to a district
    :param initialDF: the GeoDataFrame with precinct-district assignments
    :param numberSubgroupsPerBlockgroup: the number of subgroups that each district should be subdivided in
    :return: a GeoDataFrame with precinct-district mappings that represents a contiguous assignment
    '''
    assignmentsDF = initialDF.copy()
    w_df = lp.weights.Rook.from_dataframe(assignmentsDF, ids=assignmentsDF.index.tolist())
    unassigned = []
    assignedDict = {}
    precinctList = list(assignmentsDF.index)

    assignmentsDF['Group'] = -1

    # only divide if we have enough blocks in the blockgroup
    numSubgroupsPerBlockgroup = min(numberSubgroupsPerBlockgroup, initialDF.shape[0])
    # if we have more than one component, meaning that the original area is not contiguous, we don't split
    if w_df.n_components > 1:
        numSubgroupsPerBlockgroup = 1

    groupList = list(range(0, numSubgroupsPerBlockgroup))

    # Shuffle the list. shuffle works in place.
    random.shuffle(precinctList)
    for group in range(0, numSubgroupsPerBlockgroup):
        assignmentsDF.loc[precinctList[group], 'Group'] = group

    for item in assignmentsDF.index:
        if assignmentsDF.loc[item, 'Group'] < 0:
            unassigned.append(item)
        else:
            assignedDict[item] = assignmentsDF.loc[item, 'Group']

    print('number of precincts: ', len(precinctList), ' number of assigned: ', len(assignedDict),
          ' number of unassigned: ', len(unassigned))

    if numSubgroupsPerBlockgroup == 1:
        # get the initial group/district and assign this to each sub-blockgroup
        group = initialDF.District.unique()
        for precinct in precinctList:
            assignedDict[precinct] = group
            assignmentsDF.loc[precinct, 'Group'] = group
    #             unassigned.remove(precinct)
    else:

        while len(unassigned) > 0:  # We have assignments to make

            # pick a random precinct that has been already allocated
            # precinct = random.choice(list(assignedDict.keys()))
            # group = assignedDict[precinct]

            # alternatively, we can add to the smallest of the groups to aim for a more balanced distribution
            groupPops = getGroupPops(assignmentsDF, groupList)

            groupPopsSorted = {k: v for k, v in sorted(groupPops.items(), key=lambda item: item[1])}

            group_done = False
            for key in groupPopsSorted:
                group = key
                precincts = [key for key, value in assignedDict.items() if value == group]
                done = False
                for precinct in precincts:

                    # Get the neighbors of precinct
                    neighbors = list(w_df[precinct].keys())

                    # check if at least one neighbor is yet unassigned. if yes, assign
                    if set(unassigned) & set(neighbors):
                        # Shuffle, randomize the order of, the neighbors.
                        # This is done in place.
                        random.shuffle(neighbors)
                        # Go through the neighbors as shuffled. If we find
                        # an unassigned neighbor, assign it to the district
                        # of the current precinct, precinct.
                        for item in neighbors:
                            if item in unassigned:  # so if it's not there, not assigned
                                assignedDict[item] = group
                                assignmentsDF.loc[item, 'Group'] = group
                                unassigned.remove(item)
                                done = True
                                break
                    if done:
                        group_done = True
                        break
                if group_done:
                    break
    assignmentsDF['GEOID'] = assignmentsDF.index.tolist()
    assignmentsDF['GEOID'] = assignmentsDF.apply(lambda x: first_n_digits(x['GEOID'], 12), axis=1)
    assignmentsDF['GEOID'] = assignmentsDF['GEOID'] * 10 + assignmentsDF['Group']

    # then, dissolve the geometries while summing up some of the values
    constantDF = assignmentsDF[['STATEFP10', 'COUNTYFP10', 'TRACTCE10', 'BLOCKCE', 'BLOCKID10', 'PARTFLG', 'state',
                                'county', 'tract', 'block', 'TRACT', 'STATE', 'PLACE', 'BLKGRP', 'District']]
    newDF = assignmentsDF[['HOUSING10', 'POP10', 'geometry', 'total_population', 'white_population', 'black_population',
                           'latino_population', 'other_population', 'GEOID']].dissolve(by='GEOID', aggfunc='sum')
    newDF['GEOID'] = newDF.index.tolist()
    constantDF = constantDF.iloc[0:numSubgroupsPerBlockgroup]
    constantDF.index = newDF.index
    newDF = pd.concat([newDF, constantDF], axis=1)
    newDF.index = newDF['GEOID']

    return newDF


def seed_fill(assignmentsDF, w_df, numDistricts, use_original_sfsr, population_column):
    '''
    Runs the initial seed and fill step of the SFSR procedure
    First, random precincts are chosen and seeded with different district IDs
    Then, Fill assigns previously unassigned precincts to an existing district
    Notes:
        - the initial seed and fill routine creates a contiguous assignment where each precinct belongs to a district
    :param assignmentsDF: the current GeoDataFrame with precinct-district assignments (initialized)
    :param w_df: the pysal weights object with neighbor information
    :param numDistricts: the number of districts to create
    :param use_original_sfsr: should the original sfsr be used ("True") or the greedy version ("False")
    :return: a GeoDataFrame with precinct-district mappings that represents a contiguous assignment
    '''

    tempTime = datetime.datetime.now()

    unassigned = []  # This will be a list of 2-tuples (idx, -1)
    assignedDict = {}
    # randomly assign a bunch of precincts

    precinctList = list(assignmentsDF.index)

    assignmentsDF['District'] = -1

    districtList = list(range(1, numDistricts+1))

    # Shuffle the list. shuffle works in place.
    random.shuffle(precinctList)
    for district in districtList:
        assignmentsDF.loc[precinctList[district], 'District'] = district

    for item in assignmentsDF.index:
        if assignmentsDF.loc[item, 'District'] < 0:
            unassigned.append((item, assignmentsDF.loc[item, 'District']))
        else:
            assignedDict[item] = assignmentsDF.loc[item, 'District']

    global time_seed, count_seed
    time_seed += (datetime.datetime.now() - tempTime).total_seconds()
    count_seed += 1

    tempTime = datetime.datetime.now()

    while len(unassigned) > 0:  # We have assignments to make
        # Pick a random precinct already assigned.

        if use_original_sfsr == "True":
            precinct = random.choice(list(assignedDict.keys()))
            district = assignedDict[precinct]

            # Get the neighbors of precinct
            neighbors = list(w_df[precinct].keys())

            unassigned_precincts = [item[0] for item in unassigned]

            # check if at least one neighbor is yet unassigned. if yes, assign
            if set(unassigned_precincts) & set(neighbors):
                # Shuffle, randomize the order of, the neighbors.
                # This is done in place.
                random.shuffle(neighbors)
                # Go through the neighbors as shuffled. If we find
                # an unassigned neighbor, assign it to the district
                # of the current precinct, precinct.
                for item in neighbors:
                    if assignedDict.get(item, -1) < 0:  # so if it's not there, not assigned
                        assignedDict[item] = district
                        assignmentsDF.loc[item, 'District'] = district
                        unassigned.remove((item, -1))
                        break

        else:
            # alternatively, we can add to the smallest of the groups to aim for a more balanced distribution
            districtPops = getDistrictPopsDict(assignmentsDF, districtList, population_column=population_column)

            districtPopsSorted = {k: v for k, v in sorted(districtPops.items(), key=lambda item: item[1])}

            district_done = False
            for key in districtPopsSorted:
                district = key
            # district = min(range(len(districtPops)), key=districtPops.__getitem__) + 1

                # print('district pops: ', districtPops)
                # print('pick district with min population: ', district)
                precincts = [key for key, value in assignedDict.items() if value == district]
                # print('precincts: ', precincts)
                done = False
                for precinct in precincts:

                    # Get the neighbors of precinct
                    neighbors = list(w_df[precinct].keys())

                    unassigned_precincts = [item[0] for item in unassigned]

                    # check if at least one neighbor is yet unassigned. if yes, assign
                    if set(unassigned_precincts) & set(neighbors):
                        # Shuffle, randomize the order of, the neighbors.
                        # This is done in place.
                        random.shuffle(neighbors)
                        # Go through the neighbors as shuffled. If we find
                        # an unassigned neighbor, assign it to the district
                        # of the current precinct, precinct.
                        for item in neighbors:
                            if assignedDict.get(item, -1) < 0:  # so if it's not there, not assigned
                                assignedDict[item] = district
                                assignmentsDF.loc[item, 'District'] = district
                                unassigned.remove((item, -1))
                                break
                    if done:
                        district_done = True
                        break
                if district_done:
                    break

    global time_fill, count_fill
    time_fill += (datetime.datetime.now() - tempTime).total_seconds()
    count_fill += 1

    assignmentsDF['GEOID'] = assignmentsDF.index.tolist()

    # print(assignmentsDF.head())

    return assignmentsDF


def get_border_AAU_list(assignmentDF, adjMatrix, adjMatrixIDs, district):
    '''
    This function calculates the border AAUs for the specified district.
    It returns a list of GEO-IDs corresponding to the border AAUs
    :param assignmentDF: the current assignment DF
    :param adjMatrix: the adjacency matrix that has the neighbor information
    :param adjMatrixIDs: the ids for the AAUs used to indicate which row/column an AAU corresponds to
    :param district: the district ID/number
    :return: a list of GEO-IDs of border AAUs for the given district
    '''
    # get list of AAUs
    AAUList = assignmentDF[assignmentDF['District'] == district].index.tolist()

    # get precincts with neighbors in different districts
    # get indices of AAUs
    idxs = np.searchsorted(adjMatrixIDs, AAUList)
    # indexes of the precincts:
    # print('idxs', idxs)

    # rowadj is the rows of the adj matrix that are in the district
    rowadj = adjMatrix[idxs, :]

    # Note that this also gives us district:
    # distr = adj[idxs, :][:, idxs]
    # print('distr', distr)

    # idea: we change the matrix such that we delete the adjacency information of all AAUs in the district among themselves
    # the only remaining adjacencies will be neighbors in a different district
    # hence, if a row (AAU) has at least one adjacency in a different district, it's a border AAU

    # set the district columns in rowadj to 0
    rowadj[:, idxs] = 0
    # print('zeroed rowadj', rowadj)
    rowsums = rowadj.sum(1)
    # print('rowsums', rowsums)
    # print('Get the list of adj indices for district that are on the border')

    indices = idxs[rowsums > 0]
    border = [adjMatrixIDs[i] for i in indices]

    return border


def get_border_AAU_List_matrix(weightsDF, assignmentDF, district):

    # get list of AAUs
    AAUList = assignmentDF[assignmentDF['District'] == district].index.tolist()

    # get indices of AAUs
    idxs = np.searchsorted(weightsDF.id_order, AAUList)

    mat = weightsDF.sparse.tolil()

    # mat = weightsDF.sparse
    rowadj = mat[idxs, :]

    rowadj[:, idxs] = 0

    flat_list = [item for sublist in rowadj[:, :].sum(axis=1).tolist() for item in sublist]

    indices = np.nonzero(flat_list)
    indices = indices[0].tolist()

    # get the actual indices
    list_of_border_AAUs = [AAUList[i] for i in indices]

    return list_of_border_AAUs


def shift(assignedDF, idealPop, minPop, maxPop, weightsDF, Wmatrix, ids, population_column, only_nonideal_population_districts=False,
          only_contiguous_shifts=False, only_opposite_population_shifts=True, only_border_districts=True):
    '''
    At this point we have assignedDF, which records a district assignment
    for each precinct, as GEOID, as well as the populations of the
    precincts. In addition, the overall assignment should be contiguous.
    Now we need to shift things to create population balance.
    Uses: assignmentDF, which records a district assignment
    for each precinct.
    :param assignedDF: the current GeoDataFrame with precint-district assignments
    :param idealPop: the ideal population per district
    :param weightsDF: the pysal weights object with neighbor information
    :param population_column: the name of the column with population information
    :return: a new GeoDataFrame with changed precinct-district assignments
    '''
    # Make a copy and return it. This keeps the original intact.
    assignedDFTemp = assignedDF.copy()
    assignedDFTemp.index = assignedDF.index

    # # determine a random district, and iterate through its precincts until a shift can be executed
    if only_border_districts:

        # optionally, only consider the precincts in Districts that are outside the population tolerance
        if only_nonideal_population_districts:
            district_candidates = list()
            for district in assignedDFTemp.District.unique():
                temp_pop = getDistrictPopulation(district, assignedDFTemp, population_column)

                if temp_pop < minPop or temp_pop > maxPop:
                    district_candidates.append(district)

        else:
            district_candidates = [random.choice(range(1, len(assignedDFTemp.District.unique()) + 1))]

        # note: if the list is empty, we pick a random district
        if not district_candidates:
            district_candidates = [random.choice(range(1, len(assignedDFTemp.District.unique()) + 1))]

        district = random.choice(district_candidates)

        precinctList = get_border_AAU_list(assignmentDF=assignedDFTemp, adjMatrix=Wmatrix, adjMatrixIDs=ids,
                                           district=district)
        # precinctList = get_border_AAU_List_matrix(weightsDF=weightsDF, assignmentDF=assignedDFTemp,
        #                                           district=district)
    else:
        precinctList = assignedDFTemp.index.tolist()

    random.shuffle(precinctList)

    # iterate over the randomly sorted precincts until a potential shift is found. then break
    for precinct in precinctList:
        aPrecinctIdx = precinct
        aPrecinctID = assignedDFTemp.loc[aPrecinctIdx, 'GEOID']

        district = assignedDFTemp.loc[aPrecinctID, 'District']
        precinctDistrictPop = getDistrictPopulation(district, assignedDFTemp, population_column)

        theNeighbors = list(weightsDF[aPrecinctID].keys())
        # Christian: added: only consider the subset of the neighbors with a different District
        neighborsInDifferentDistrict = assignedDFTemp[
            assignedDFTemp['GEOID'].isin(theNeighbors)]
        neighborsInDifferentDistrict = \
            neighborsInDifferentDistrict[neighborsInDifferentDistrict['District'] != district]['GEOID'].tolist()
        if neighborsInDifferentDistrict:
            # Now shuffle them to create a random ordering
            # randomize the rows

            random.shuffle(neighborsInDifferentDistrict)
            # Now we go through them and switch at most one.
            for neighbor in neighborsInDifferentDistrict:
                neighborDistrict = assignedDFTemp.loc[neighbor]['District']
                neighborDistrictPop = getDistrictPopulation(neighborDistrict, assignedDFTemp, population_column)

                precinctSmaller = None
                if only_opposite_population_shifts:
                    if precinctDistrictPop <= idealPop and neighborDistrictPop >= idealPop:
                        precinctSmaller = True
                    elif precinctDistrictPop >= idealPop and neighborDistrictPop <= idealPop:
                        precinctSmaller = False
                else:
                    if precinctDistrictPop <= neighborDistrictPop :
                        precinctSmaller = True
                    elif precinctDistrictPop >= neighborDistrictPop:
                        precinctSmaller = False

                if precinctSmaller and precinctSmaller is not None:
                    # Put the neighborDistrict into the precinct's district
                    if only_contiguous_shifts:
                        # check if the new district is contiguous
                        precincts = assignedDFTemp[assignedDFTemp['District'] == neighborDistrict].GEOID.unique().tolist()
                        precincts.remove(neighbor)

                        contiguous = fastContiguityCheckForChangedDistrict(weights=weightsDF,
                                                                           listOfPrecincts=precincts)

                        if contiguous:
                            assignedDFTemp.loc[neighbor, 'District'] = district
                        else:
                            continue
                    else:
                        assignedDFTemp.loc[neighbor, 'District'] = district
                    return assignedDFTemp

                # Mirrows the above if...
                if not precinctSmaller and precinctSmaller is not None:

                    # Put the precinct's district into the neighborDistrict
                    if only_contiguous_shifts:
                        # check if the new district is contiguous
                        precincts = assignedDFTemp[assignedDFTemp['District'] == district].GEOID.unique().tolist()
                        precincts.remove(aPrecinctIdx)

                        contiguous = fastContiguityCheckForChangedDistrict(weights=weightsDF,
                                                                           listOfPrecincts=precincts)

                        if contiguous:
                            assignedDFTemp.loc[aPrecinctIdx, 'District'] = neighborDistrict

                        else:
                            continue
                    else:
                        #                     assignedDFTemp.loc[aPrecinctIdx, 'District'] = neighborDistrictRow['District']
                        assignedDFTemp.loc[aPrecinctIdx, 'District'] = neighborDistrict
                    #                     print('regular shift done')
                    return assignedDFTemp
    return assignedDFTemp


def shift_original(assignedDF, idealPop, weightsDF, population_column):
    '''
    At this point we have assignedDF, which records a district assignment
    for each precinct, as GEOID, as well as the populations of the
    precincts. In addition, the overall assignment should be contiguous.
    Now we need to shift things to create population balance.
    Uses: assignmentDF, which records a district assignment
    for each precinct.
    :param assignedDF: the current GeoDataFrame with precint-district assignments
    :param idealPop: the ideal population per district
    :param weightsDF: the pysal weights object with neighbor information
    :param population_column: the name of the column with population information
    :return: a new GeoDataFrame with changed precinct-district assignments
    '''
    # Make a copy and return it. This keeps the original intact.
    assignedDFTemp = assignedDF.copy()
    assignedDFTemp.index = assignedDF.index

    precinctList = assignedDFTemp.index.tolist()

    random.shuffle(precinctList)

    # iterate over the randomly sorted precincts until a potential shift is found. then break
    for precinct in precinctList:

        aPrecinctIdx = precinct
        aPrecinctID = assignedDFTemp.loc[aPrecinctIdx, 'GEOID']

        district = assignedDFTemp.loc[aPrecinctID, 'District']
        precinctDistrictPop = getDistrictPopulation(district, assignedDFTemp,population_column=population_column)

        theNeighbors = list(weightsDF[aPrecinctID].keys())

        neighborsInDifferentDistrict = assignedDFTemp[
            assignedDFTemp['GEOID'].isin(theNeighbors)]
        neighborsInDifferentDistrict = \
            neighborsInDifferentDistrict[neighborsInDifferentDistrict['District'] != district]['GEOID'].tolist()

        if neighborsInDifferentDistrict:
            # Now shuffle them to create a random ordering
            # randomize the rows

            random.shuffle(neighborsInDifferentDistrict)
            # Now we go through them and switch at most one.
            for neighbor in neighborsInDifferentDistrict:
                #             neighborDistrictPop = getDistrictPopulation(neighborDistrictRow['District'], assignedDFTemp)
                neighborDistrict = assignedDFTemp.loc[neighbor]['District']
                neighborDistrictPop = getDistrictPopulation(neighborDistrict, assignedDFTemp,population_column)

                precinctSmaller = None
                if precinctDistrictPop <= idealPop and neighborDistrictPop >= idealPop:
                    precinctSmaller = True
                elif precinctDistrictPop >= idealPop and neighborDistrictPop <= idealPop:
                    precinctSmaller = False

                if precinctSmaller and precinctSmaller is not None:
                    assignedDFTemp.loc[neighbor, 'District'] = district
                    return assignedDFTemp

                # Mirrows the above if...
                if not precinctSmaller and precinctSmaller is not None:
                    assignedDFTemp.loc[aPrecinctIdx, 'District'] = neighborDistrict
                    return assignedDFTemp
    return assignedDFTemp


def checkDistrictContiguity(districtNum,districtAssignmentDF):
    '''
    Given a district number and a district assignment GEODataFrame,
    reports back on the
    contiguity of the district.
    Sometimes given areal units, e.g., census tracts, are discontiguous.
    We are not checking for that.
    repairList form:
    [22003950300, 22005030900], False/True
    :param districtNum: the district number for which contiguity needs to be checked
    :param districtAssignmentDF: the current GeoDataFrame with precinct-district assignments and geographic information
    :return: either a tuple of [], True in case the district is contiguous, otherwise tuple [GEOIDs], False of
    a list of GEOIDs that cause the discontiguity which can be used for repair
    '''

    # note: adding a warning filter here as otherwise everytime pysal finds an island it prints a warning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Get a DataFrame for the district we wish to examine.
        aDistrictDF = districtAssignmentDF[districtAssignmentDF['District'] == districtNum]

        w_df = lp.weights.Rook.from_dataframe(aDistrictDF)

        if w_df.n_components == 1: # Every unit is reachable!!
            return ([],True)
        else: # Not every unit is reachable, so not contiguous; return the IDs

            # if we have more than one component, we need to add all smaller components to the repair list
            # first, identify the largest component

            u, count = np.unique(w_df.component_labels, return_counts=True)

            count_sort_ind = np.argsort(-count)

            # now, the first count_sort_ind has the component with the largest number of units. we can use the rest to be added to the repair list
            toReturn = []
            largest_component_label = u[count_sort_ind][0]
            smaller_component_index = np.where(w_df.component_labels != largest_component_label)

            # for each district in the non-connected components, add them to the repair list
            for i, ix in enumerate(smaller_component_index[0]):
                # print((aDistrictDF.index[i], aDistrictDF.iloc[ix]['GEOID']))
                toReturn.append(aDistrictDF.iloc[ix]['GEOID'])

            return (toReturn,False) #np.where(connectiv[0,:] ==0)[1] # list(connectiv[1])


def checkPlanContiguity(anAssignmentDF,districtList):
    '''
    Checks the contiguity of the entire plan, i.e., for all districts
    :param anAssignmentDF: a GeoDataFrame with the current precinct-district assignments
    :param districtList: the list of district IDs
    :return: True if all districts are contiguous, False otherwise
    '''
    for district in districtList:
        connectiv = checkDistrictContiguity(district,anAssignmentDF)
        if connectiv[1] is False:
            return False
    return True


def shiftWhile(toShiftDF,idealPop,maxPop,minPop, topMaxPerCent,bottomMinPerCent,districtList,weightsDF,verbose, population_column, remainingTime, use_original_sfsr):
    '''
    The main function that calls the shift() procedure until the district populations are within
    the population tolerance
    :param toShiftDF: the current GeoDataFrame with precinct-district assignments
    :param idealPop: the ideal population per district
    :param maxPop: the maximum allowed population per district
    :param minPop: the minimum allowed population per district
    :param topMaxPerCent: the maximum positive percent deviation from the ideal population
    :param bottomMinPerCent: the maximum minimum percent deviation from the ideal population
    :param districtList: the list of districts
    :param weightsDF: the pysal weights object with neighbor information
    :param verbose: argument if additional printouts should be provided
    :param population_column: the name of the population column
    :return: an updated assignment of precinct-district mapping in a GeoDataFrame
    '''
    tempTime = datetime.datetime.now()

    if verbose:
        print('Starting shiftWhile() at {}'.format(str(datetime.datetime.now())))
        print('District populations: {}'.format(getDistrictPops(toShiftDF,districtList, population_column)))
    count = 0

    Wmatrix, ids = weightsDF.full()

    startTime = time.time()

    while time.time() <= startTime + remainingTime:
    # while not popSizeOK(toShiftDF,idealPop,topMaxPerCent,bottomMinPerCent,districtList):
    #     print('current iteration: ', count)
        count += 1
        if verbose:
            if count % 300 == 0:
                print(count)
                balance = calculate_population_equality(toShiftDF, True, population_column)
                print('District populations: {}'.format(getDistrictPops(toShiftDF,districtList, population_column)))
                print('current population balance: ', balance)
                lower, upper = calculate_population_tolerance(toShiftDF, population_column)
                print('current upper threshold: ', upper, ' and lower threshold: ', lower)
                # if balance < 0.020:
                #     break
        if use_original_sfsr == "True":
            tempDF = shift(toShiftDF, idealPop, minPop, maxPop, weightsDF, Wmatrix, ids, population_column,
                       only_nonideal_population_districts=False,
                       only_contiguous_shifts=False, only_opposite_population_shifts=True, only_border_districts=False)
        else:
            # update: use random numbers instead of fixed sequences
            # rnd = random.random()

            # if rnd < 1/7:
            if count > 1000 and count % 7 == 0:
                tempDF = shift(toShiftDF,idealPop,minPop, maxPop, weightsDF, Wmatrix, ids, population_column,
                               only_nonideal_population_districts=True,
                               only_contiguous_shifts=True)
            elif count > 1000 and count % 3 == 0:
            # elif rnd < 1/3:
                tempDF = shift(toShiftDF,idealPop,minPop, maxPop, weightsDF, Wmatrix, ids, population_column,
                               only_nonideal_population_districts=True,
                               only_contiguous_shifts=True,
                               only_opposite_population_shifts=False)
            elif count > 5000 and count % 1000 == 0:
            # elif rnd > 0.999:
                local_count_shift = 0
                tempDF = toShiftDF
                # print('current population balance: ', calculate_population_equality(toShiftDF))
                while local_count_shift < 50:
                    tempDF = shift(tempDF,idealPop,minPop,maxPop,weightsDF, Wmatrix, ids, population_column, True, False, True)
                    local_count_shift += 1
                    count += 1

                contiguous = checkPlanContiguity(tempDF,
                                                 districtList)

                # If we get here, contiguous == False
                while not contiguous:
                    (tempDF, oldDF) = repairFor(tempDF, districtList, weightsDF, verbose, population_column)

                    contiguous = checkPlanContiguity(tempDF,districtList)
            else:
                tempDF = shift(toShiftDF, idealPop, minPop, maxPop, weightsDF, Wmatrix, ids, population_column,
                               only_nonideal_population_districts=False,
                               only_contiguous_shifts=True)

        toShiftDF = tempDF

        if popSizeOK(toShiftDF,idealPop,topMaxPerCent,bottomMinPerCent,districtList, population_column):
            break

    if verbose:
        print('District populations: {}'.format(getDistrictPops(toShiftDF,districtList, population_column)))
        print('idealPop = {}, maxPop = {}, minPop = {}'.format(idealPop,maxPop,minPop))

    contiguity=checkPlanContiguity(toShiftDF,districtList)

    if verbose:
        print("Contiguous? {} Done with shiftWhile()".format(contiguity))

    global time_shift, count_shift
    time_shift += (datetime.datetime.now() - tempTime).total_seconds()
    count_shift += 1

    return toShiftDF


def getQualifiedNeighbors(districtNum, neighbors, assignedDF):
    '''
    Gets neighbor precincts that would be able to repair the missing contiguity
    :param districtNum: the district number
    :param neighbors: the neighbor information from pysal
    '''

    neighboringQualifiedDistricts = []

    defaultGeoid = neighbors[0]

    theIndex = assignedDF[assignedDF['GEOID'] == defaultGeoid].index[0]
    defaultDistrict = assignedDF.loc[theIndex, 'District']
    for geoid in neighbors:
        theIndex = assignedDF[assignedDF['GEOID'] == geoid].index[0]
        district = assignedDF.loc[theIndex, 'District']
        if district != districtNum:
            neighboringQualifiedDistricts.append((geoid, district))
        if neighboringQualifiedDistricts == []:  # you are surrounded by your own
            neighboringQualifiedDistricts.append((defaultGeoid, defaultDistrict))

    return neighboringQualifiedDistricts


def repairOneDistrict(districtNum,assignedDF,weightsDF):
    '''
    Repairs one discontiguous district
    repairList comes from checkDistrictContiguity, e.g.,
    [22003950300, 22005030900]
    Start here with only 1 item
    :param districtNum: the district number
    :param assignedDF: the current GeoDataFrame with precint-district assignments
    :param weightsDF: the pysal weights object with neighbor information
    :return: a repaired GeoDataFrame with updated precicnt-district
    '''
    contiguity = checkDistrictContiguity(districtNum,assignedDF)
    repairList = contiguity[0]

    if repairList == []:
        return assignedDF
    workingAssignedDF = assignedDF.copy()

    for element in repairList: # later, relax this

        geoIDToRepair = element
        rowToRepair = assignedDF[assignedDF['GEOID']==geoIDToRepair].index[0]

        neighbors = list(weightsDF[element].keys())
        qualifiedNeighbors = getQualifiedNeighbors(districtNum,neighbors,workingAssignedDF)

        toChange = random.choice(qualifiedNeighbors)

        toChangeTo = toChange[1]

        workingAssignedDF.loc[rowToRepair,'District'] = toChangeTo

    return workingAssignedDF


def repairFor(toRepairDF,
              districtList, weightsDF, verbose, population_column):
    '''
    A function that loops through the different districts to repair each one of them (if necessary)
    :param toRepairDF: the current GeoDataFrame with precinct-district mappings, potentially discontiguous
    :param districtList: the list of districts
    :param weightsDF: the pysal weights object with neighbor information
    :param verbose: parameter if additional printouts should be provided
    :param population_column: the name of the population column
    :return: an updated GeoDataFrame with repaired districts
    '''
    tempTime = datetime.datetime.now()

    if verbose:
        print('**** Starting repairFor() ****')
        print('District populations: {}'.format(getDistrictPops(toRepairDF, districtList, population_column)))

    random.shuffle(districtList)

    for districtNum in districtList:
        toRepairDF = repairOneDistrict(districtNum, toRepairDF, weightsDF)
        # checkDistrictContiguity(districtNum, toRepairDF)

    global time_repair, count_repair
    time_repair += (datetime.datetime.now() - tempTime).total_seconds()
    count_repair += 1

    return (toRepairDF, toRepairDF.copy())  # the idea is (new,old)


def shiftRepair(toShiftRepairDF,
                idealPop, maxPop, minPoP,
                topMaxPerCent, bottomMinPerCent, districtList, weightsDF, use_original_sfsr, verbose, population_column):
    '''
    The main shift-repair loop until a contiguous solution within the population tolerance is found.
    repair() then shift() if not ok on population.
    Then repair() if not contiguous. Embedded in a
    while loop.
    :param toShiftRepairDF: the current GeoDataFrame with precinct-district information
    :param idealPop: the ideal population per district
    :param maxPop: the maximum population per district
    :param minPoP: the minimum population per district
    :param topMaxPerCent: the maximum positive deviation from the ideal population
    :param bottomMinPerCent: the maximum negative deviation from the ideal population
    :param districtList: the list of districts
    :param weightsDF: the pysal weights object with neighbor information
    :param verbose: parameter if additional details should be printed
    :param population_column: the name of the population column
    :return: a GeoDataFrame with contiguous districts within the population tolerance
    '''
    if verbose:
        print('Starting shiftRepair at {}'.format(str(datetime.datetime.now())))

    timeout = 60 * 60 * 48  # 48 hours time limit
    if use_original_sfsr == "True":
        timeout = 60 * 60 * 10000000  # no time limit

    orig_df = copy.deepcopy(toShiftRepairDF)

    done = False
    while not done:
        startTime = time.time()
        toShiftRepairDF = seed_fill(orig_df, weightsDF, len(districtList), use_original_sfsr, population_column)

        while time.time() < startTime + timeout:
            # print('Top of the while loop of shiftRepair. District populations: {}'.format(getDistrictPops(toShiftRepairDF,districtList)))
            remainingTime = timeout - (time.time() - startTime)
            toShiftRepairDF = shiftWhile(toShiftRepairDF, idealPop, maxPop, minPoP,
                                         topMaxPerCent,
                                         bottomMinPerCent, districtList, weightsDF, verbose, population_column, remainingTime, use_original_sfsr)
            contiguous = checkPlanContiguity(toShiftRepairDF, districtList)

            if contiguous and popSizeOK(toShiftRepairDF,idealPop,topMaxPerCent,bottomMinPerCent,districtList, population_column):  # Then we are done!!
                done = True  # Not needed
                return toShiftRepairDF
            # If we get here, contiguous == False
            while not contiguous:
                (toShiftRepairDF, oldDF) = repairFor(toShiftRepairDF, districtList, weightsDF, verbose, population_column)

                contiguous = checkPlanContiguity(toShiftRepairDF, districtList)
            # We are contiguous if we get here.

            if popSizeOK(toShiftRepairDF, idealPop, topMaxPerCent, bottomMinPerCent, districtList, population_column):
                done = True
                break
                # print('**********Success from repairFor(): Finishing repairFor() at {}'.format(str(datetime.datetime.now())))
            # print('Bottom of the while loop of shiftRepair. District populations: {}'.format(getDistrictPops(toShiftRepairDF,districtList)))

        # new: if the returned assignmentDF is not contiguous and /or not within population tolerance, repeat process
        if not popSizeOK(toShiftRepairDF, idealPop, topMaxPerCent, bottomMinPerCent, districtList, population_column) or \
            not checkPlanContiguity(toShiftRepairDF, districtList):
            toShiftRepairDF = seed_fill(orig_df, weightsDF, len(districtList), use_original_sfsr, population_column)

    # print('*****Success: Finishing shiftRepair() at {}'.format(str(datetime.datetime.now())))
    return toShiftRepairDF


def go(idealPop,minPop,maxPop,numDistricts,topMaxPerCent,bottomMinPerCent,
       districtList,assignmentsDF, use_original_sfsr,weightsDF,
       verbose,population_column, seed=None):
    '''
    The key procedure.  First we initialize, then
    seed, then fill, then shiftRepair until done.
    :param idealPop: the ideal population per district
    :param minPop: the minimum population allowed per district
    :param maxPop: the maximum population allowed per district
    :param numDistricts: the number of districts to be created
    :param topMaxPerCent: the maximum positive deviation from the ideal population
    :param bottomMinPerCent: the maximum negative deviation from the ideal population
    :param districtList: the list of districts
    :param assignmentsDF: a GeoDataFrame with precinct-district assignments
    :param verbose: parameter if additional information should be printed
    :param seed: an optional random number seed (currently not used)
    :param population_column: the name of the population column
    :return: a GeoDataFrame with precinct-district information where each district is contiguous and within the
    population tolerance
    '''

    if verbose:
        print('Welcome to go()')

    print('starting seed fill')
    # assignmentsDF = seed_fill(assignmentsDF_bg, assignmentsDF, weightsDF, numDistricts)
    # print('level: ', level)

    # weightsDF = lp.weights.Rook.from_dataframe(assignmentsDF, ids=assignmentsDF.index.tolist())

    # u, count = np.unique(weightsDF.component_labels, return_counts=True)
    #
    # count_sort_ind = np.argsort(-count)
    #
    # # now, the first count_sort_ind has the component with the largest number of units. we can use the rest to be added to the repair list
    # toReturn = []
    # largest_component_label = u[count_sort_ind][0]
    # smaller_component_index = np.where(weightsDF.component_labels != largest_component_label)
    #
    # # for each district in the non-connected components, add them to the repair list
    # for i, ix in enumerate(smaller_component_index[0]):
    #     # print((aDistrictDF.index[i], aDistrictDF.iloc[ix]['GEOID']))
    #     toReturn.append([assignmentsDF.iloc[ix]['GEOID'], weightsDF.component_labels[ix]])
    #
    # print('largest component label: ', largest_component_label , ' with how many units? ', len(np.where(weightsDF.component_labels == largest_component_label)))
    # print('smaller component: ', toReturn)

    # assignmentsDF = seed_fill(assignmentsDF, weightsDF, numDistricts, use_original_sfsr, population_column)

    ###################### shiftRepair #######################
    # The heart of the procedure.
    toReturnDF = shiftRepair(assignmentsDF, idealPop, maxPop, minPop, topMaxPerCent,
                             bottomMinPerCent,districtList,weightsDF,use_original_sfsr,verbose, population_column)
    return toReturnDF


class SFSR():

    def __init__(self,seed=None,numRuns=1,pathToData = 'Data/PA/',
                 nDistricts = 17,tolerance=0.05,
                 neighborsType='rook',
                 pathToOutputData='Solutions/PA/SFSR/Tracts/',
                 cvap='False', level='blockgroup', data_table=None, shapefile_id_column=None, data_id_column=None,
                 population_column='total_population', use_original_sfsr='True'
                 ):

       self.numRuns = numRuns
       self.seed = seed
       self.pathToData = pathToData
       self.numDistricts = nDistricts
       self.tolerance = tolerance
       self.topMaxPerCent = self.tolerance
       self.bottomMinPerCent = self.tolerance
       self.neighborsType = neighborsType
       self.pathToOutputData = pathToOutputData
       self.cvap = cvap
       self.level = level
       self.data_table = data_table

       if shapefile_id_column is None or shapefile_id_column == "":
           self.shapefile_id_column = 'GEOID10'
       else:
           self.shapefile_id_column = shapefile_id_column
       if data_id_column is None or data_id_column == "":
            self.data_id_column = 'GEOID10'
       else:
           self.data_id_column = data_id_column
       if population_column is None or population_column == "":
           self.population_column = 'total_population'
       else:
           self.population_column = population_column
       self.use_original_sfsr = use_original_sfsr
       self.version = 1.0

        # initialize
       self.setup()

       if not (self.neighborsType == 'rook' or self.neighborsType == 'queen'):
           print(f"neighborsType = {self.neighborsType} should be 'rook' or 'queen'")

    def getInits(self):
        toReturn = ''
        toReturn += f'version={self.version}\n'
        toReturn += f'numRuns={self.numRuns}\n'
        toReturn += f'seed={self.seed}\n'
        toReturn += f'pathToData={self.pathToData}\n'
        toReturn += f'numDistricts={self.numDistricts}\n'
        toReturn += f'tolerance={self.tolerance}\n'
        toReturn += f'topMaxPerCent={self.topMaxPerCent}\n'
        toReturn += f'bottomMinPerCent={self.bottomMinPerCent}\n'
        toReturn += f'neighborsType={self.neighborsType}\n'
        toReturn += f'pathToOutputData={self.pathToOutputData}\n'
        return toReturn

    def setup(self):
        '''
        Initialization of the SFSR procedure.
        Sets up files and runs the procedure.
        '''
        # prepare the initial data from the specified locations

        self.data_table['GEOID'] = self.data_table.index

        self.assignmentsDF = copy.deepcopy(self.data_table)

        self.numPrecincts = len(self.assignmentsDF.index)

        self.idealPop = self.assignmentsDF[self.population_column].sum()/self.numDistricts
        self.maxPop = self.idealPop * (1 + self.topMaxPerCent)
        self.minPop = self.idealPop * (1 - self.bottomMinPerCent)
        self.districtList = list(range(1,self.numDistricts+1))

        self.weightsDF = lp.weights.Rook.from_dataframe(self.data_table, ids=self.data_table.index.tolist())

        random.seed(self.seed)

        # print(self.assignmentsDF.index)

        print("SFSR Setup completed.")

    def getSettings(self):
        toReturn = ''
        toReturn += f'numPrecincts={self.numPrecincts}\n'
        toReturn += f'idealPop={self.idealPop}\n'
        toReturn += f'maxPop={self.maxPop}\n'
        toReturn += f'minPop={self.minPop}\n'
        toReturn += f'districtList={self.districtList}\n'
        return toReturn

    def basicgo_single(self, run, verbose=True):
        '''
        Runs a single iteration of the SFSR procedure
        :param run: the run ID
        :param verbose: parameter if additional information should be printed
        '''
        global time_seed, count_seed, time_fill, count_fill, time_shift, count_shift, time_repair, count_repair
        time_seed = 0
        time_fill = 0
        time_shift = 0
        time_repair = 0

        count_seed = 0
        count_fill = 0
        count_shift = 0
        count_repair = 0

        # self.setup()
        startTime = str(datetime.datetime.now())
        # Make copies for a single run of the various DataFrames read in
        print(f'Start time = {startTime}')

        newDF = go(self.idealPop, self.minPop, self.maxPop,
                   self.numDistricts, self.topMaxPerCent, self.bottomMinPerCent,
                   self.districtList, self.assignmentsDF,
                   self.use_original_sfsr, self.weightsDF,
                   verbose, self.population_column, seed=self.seed)

        stopTime = str(datetime.datetime.now())
        stopTimex = re.sub(r'(:|\.| )', 'x', stopTime)
        startTimex = re.sub(r'(:|\.| )', 'x', startTime)

        # columns_to_save = ['GEOID', 'District', self.population_column]
        columns_to_save = ['District', self.population_column, self.shapefile_id_column]

        newDF = newDF[columns_to_save]

        save_string = '-SFSROrig'
        if self.use_original_sfsr == "False":
            save_string = '-SFSRGreedy'

        if self.cvap == "True":
            save_string = save_string + '-cvap'

        newDF.to_csv(self.pathToOutputData + os.sep + 'df' + startTimex + '-' + stopTimex + 'pc' + str(
            self.topMaxPerCent) + save_string + '.csv')

        print(f'Stop time = {stopTime}')
        if verbose:
            print(startTime, 'startTime')
            print('minPop {}, maxPop {}, total pop1 {}, total pop2 {}'.format(self.minPop, self.maxPop,
                                                                              self.assignmentsDF[self.population_column].sum(), newDF[
                                                                                  self.population_column].sum()))
        # write the runtime information to a file
        runtime_path = self.pathToOutputData[:12] + '/Runtimes.csv'
        # if the file exists, append, otherwise create the file
        if os.path.isfile(runtime_path):
            runtimes = pd.read_csv(runtime_path)
        else:
            runtimes = pd.DataFrame(columns=['run_id', 'num_dist', 'time_seed', 'count_seed', 'time_fill', 'count_fill',
                                             'time_shift', 'count_shift', 'time_repair', 'count_repair'])
        try:
            runtimes.loc[len(runtimes)] = [self.pathToOutputData, self.numDistricts, time_seed, count_seed, time_fill, count_fill,
                                           time_shift, count_shift, time_repair, count_repair]
            runtimes.to_csv(runtime_path, index=False)
        except:
            print('Error during writing runtime results. Check if file layout is correct.')


    def basicgo_parallel(self, num_runs, pool=None, verbose=True):
        '''
        Procedure that runs independent SFSR procedures based on an available multiprocessing pool
        :param num_runs: the number of SFSR runs that should be run
        :param pool: the multiprocessing pool
        :param verbose: parameter if additional information should be printed
        '''
        pool.map(partial(self.basicgo_single, verbose=verbose), range(num_runs))
