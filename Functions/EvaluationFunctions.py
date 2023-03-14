"""

This file is a basic test file to create and check various functions for the redistricting problem

Created on Wed 2019-10-30

"""

import libpysal as lp
import numpy as np
from statistics import mean
import geocompactness.compactness_measures as cm
import geopandas as gp
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon
import copy
from decimal import DivisionByZero

def calculate_compactness_boundary_lengths(district_df_agg, state_perimeter):
    """
    This function calculates the compactness as defined in Bozkaya 2003 for a specific district
    :param district_df_agg: the dataframe that maps the census tract IDs with the district number
    :param state_perimeter: the perimeter of the entire state
    :return: the compactness score for the specified district
    """

    try:
        # note: we use km, not m here
        district_lengths = district_df_agg.length / 1000

        local_state_perimeter = state_perimeter / 1000

        if local_state_perimeter != 0:
            compactness = (sum(district_lengths) - local_state_perimeter) / \
                2 * local_state_perimeter
            return compactness
        else:
            raise DivisionByZero
    except DivisionByZero:
        print("Division by zero exception while calculating the boundary lengths: ")



def calculate_compactness_radial(district_df, district_df_agg):
    """
    This function calculates the compactness as defined in Bacao 2005 for a specific district
    :param district_df: the dataframe that maps the census tract IDs with the district number
    :return: the compactness score for the specified district
    """
    # calculate the centroids
    centroids = district_df_agg.centroid

    overall_distance = 0
    unique_districts = district_df.District.unique()
    # if we manually removed some precincts, e.g., islands, we need to make sure the unassigned districts don't cause a NaN error
    unique_districts = [x for x in unique_districts if str(x) != 'nan']

    for district_id in unique_districts:
        # the temp_centroid is the overall centroid for the current district
        temp_centroid = centroids[district_id]
        temp_df = district_df[district_df.District == district_id]

        # the temp_df centroid are the centroids for the individual tracts
        for district in temp_df.centroid:
            # note: we use km, not m here
            overall_distance += (temp_centroid.distance(district) / 1000)

    return overall_distance


def calculate_compactness_measures(district_df_agg):
    """
    This function callls the compactness_measures model to calculate various compactness measures
    :param district_df_agg: the dataframe that maps the census tract IDs with the district number
    :return: various compactness score for the specified district
    """

    compactness_polsby_popper = mean(cm.polsby_popper(district_df_agg))
    compactness_schwartzberg = mean(cm.schwartzberg(district_df_agg))
    compactness_reock = mean(cm.reock(district_df_agg))
    compactness_convex_hull = mean(cm.c_hull_ratio(district_df_agg))
    compactness_moment_of_inertia = mean(cm.moment_of_inertia(district_df_agg))

    return [compactness_polsby_popper, compactness_schwartzberg, compactness_reock, compactness_convex_hull,
            compactness_moment_of_inertia]


def check_contiguity_district(num_district, district_df):
    """
    This function checks if the specified district is contiguous or not
    :param num_district: the specific district for which compactness is calculated
    :param district_df: the dataframe that maps the census tract IDs with the district number
    :return: boolean true/false if the specified district is contiguous or not
    """

    w_d1 = lp.weights.Rook.from_dataframe(district_df[district_df.District == num_district])

    if w_d1.n_components == 1:
        return True
    else:
        return False


def check_contiguity(district_df):
    """
    This function checks if the specified district is contiguous or not
    :param district_df: the dataframe that maps the census tract IDs with the district number
    :return: boolean true/false if the specified district is contiguous or not
    """
    unique_districts = district_df.District.unique()
    for district_id in unique_districts:
        contiguous = check_contiguity_district(district_id, district_df)
        if not contiguous:
            return False

    return True


def calculate_socioeconomic_homogeneity(district_df):
    """
    This function calculates the socioeconomic homogeneity as defined in Bozkaya 2003 for a specific district
    :param district_df: the dataframe that maps the census tract IDs with the district number
    :return: the socioeconomic homogeneity score for the specified district
    """
    mean_income = district_df.median_income.mean()
    temp_df = district_df[['District', 'median_income']].groupby('District')['median_income'].agg(
        mean_median_income='mean', std_median_income='std')

    homogeneity = temp_df.std_median_income.sum() / mean_income

    return homogeneity


def calculate_population_equality(district_df, in_percent=True, population_column='total_population'):
    """
    This function calculates the socioeconomic homogeneity similar to Bozkaya 2003 for a specific district
    In contrast to Bozkaya, we record the actual population equality, but do not consider solutions within a certain
    limit as 0
    :param district_df: the dataframe that maps the census tract IDs with the district number
    :return: the socioeconomic homogeneity score for the specified district
    """

    population_per_district = []
    unique_districts = district_df.District.unique()
    for district_num in unique_districts:
        temp_df = district_df[district_df.District == district_num]
        population_per_district.append(temp_df[population_column].sum())

    average_population = mean(population_per_district)

    equality = 0

    for district_pop in population_per_district:
        equality += abs(district_pop - average_population)

    if in_percent:
        return equality / (average_population * len(district_df.District.unique()))
    else:
        return equality


def calculate_num_counties_split(district_df, in_percent=False):
    """
    This function calculates the snumber of counties that are split into several districts
    :param district_df: the dataframe that maps the census tract IDs with the district number
    :param in_percent: should the result be returned in percent or absolute numbers
    :return: the number or percent of counties that are split into several districts
    """

    # first, aggregate the district_df table by county ID
    temp_df = district_df[['county', 'District']].groupby('county')['District'].nunique()

    # then, count the countyIDs for which we have more than one district
    count = len(temp_df[temp_df > 1])

    # finally, return either the percentage or the actual number of counties with more than one district
    if in_percent:
        num_counties = len(district_df.county.unique())
        return count / num_counties
    else:
        return count


def calculate_districts_per_county(district_df, name_of_index='county'):
    """
    This function calculates the number of districts that a county is split into
    :param district_df: the dataframe that maps the census tract IDs with the district number
    # :return: a dictionary that gives the number of districts (value) for each county (key)
    :return: an aggregated dataframe with the number of districts per county
    """

    # first, aggregate the district_df table by county ID
    temp_df = district_df[[name_of_index, 'District']].groupby(name_of_index)['District'].nunique()

    return temp_df


def check_population_tolerance(district_df, tolerance=0.05, population_column='total_population'):
    """
    This function checks if the current solution is within a certain population tolerance limit for each district
    :param district_df: the solution candidate for the redistricting solution
    :param tolerance: the tolerance that each district is allowed to have compared to the average population per
    district
    :return: boolean yes/no if the candidate solution is within the tolerance limit
    """

    average_population = district_df[population_column].sum() / len(district_df.District.unique())

    upper_tolerance = average_population * (1 + tolerance)
    lower_tolerance = average_population * (1 - tolerance)

    district_populations = district_df[['District', population_column]].groupby('District')[
            population_column].sum()

    for pop in district_populations:
        if pop > upper_tolerance or pop < lower_tolerance:
            return False
    return True


def calculate_population_tolerance(district_df, population_column='total_population'):
    """
    This function calculates the lower and upper deviation to the average population
    :param district_df: the solution candidate for the redistricting solution
    :return: lower, upper actual tolerance
    """

    average_population = district_df[population_column].sum() / len(district_df.District.unique())
    district_populations = district_df[['District', population_column]].groupby('District')[population_column].sum()
    lower = 1 - min(district_populations) / average_population
    upper = 1 - max(district_populations) / average_population

    return lower, upper


def calculate_similarity_to_existing_districting(district_df_agg, congress_districts, state_area):
    """
    This function calculates a similarity score that represents how close the districts are to the existing solution.
    For this, each district is compared to the most overlapping district in the existing solution.
    For now, this follows the similarity score suggested by Bozkaya 2003
    :return:
    """
    # print('crs of district_df: ', district_df.crs)
    # print('crs of congress_district: ', congress_districts.crs)
    # print('number of congress districts: ', congress_districts.shape[0])

    # district_df_agg = district_df.dissolve(by='District')
    # print('crs of district_df_agg: ', district_df_agg.crs)
    overlaps = []
    max_overlaps = []

    if congress_districts is not None:

        # print('state area: ', state_area)

        congress_districts['area'] = congress_districts['geometry'].area

        for i, congress_district in congress_districts.iterrows():
            # print('i: ', i, ' and congress district: ', congress_district)
            temp_overlaps = []
            for j, district in district_df_agg.iterrows():
                # print(district['geometry'])
                # print(congress_district['geometry'])
                # note: sometimes geometries are not perfect and could lead to error messages due to self-intersection
                # hence, if we encounter an error we set the overlap to 0 for now
                try:
                    intersect = congress_district['geometry'].intersection(district['geometry'])
                    # print(intersect.area)
                    # temp_overlaps.append(intersect.area / congress_district['area'])
                    temp_overlaps.append(intersect.area / state_area)
                except:
                    temp_overlaps.append(0)
            overlaps.append(temp_overlaps)
            temp_max_overlap = max(temp_overlaps)
            max_overlaps.append(temp_max_overlap)
            # print('max overlap: ', temp_max_overlap)

        # print('similarity score: ', (1-sum(max_overlaps)))
        return 1 - sum(max_overlaps)

    else:
        return 0


# Faster than is_pareto_efficient_simple, but less readable.
# code from: https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
def is_pareto_efficient(costs, return_mask=True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index < len(costs):
        nondominated_point_mask = np.any(costs < costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index]) + 1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype=bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient


# Fairly fast for many datapoints, less fast for many costs, somewhat readable
# code from: https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
def is_pareto_efficient_simple(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype=bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient] < c, axis=1)  # Keep any point with a lower cost
            is_efficient[i] = True  # And keep self
    return is_efficient


def getMMstats(solution_candidate, column_suffix=None, population_column='total_population'):
    '''
    Returns a 4-tuple of district counts:
    (white percent >= 75, white percent >= 50,
    black percent >= 50, black percent >= 40).
    White excludes black, latino, and other.
    '''

    if column_suffix is None:
        column_suffix = ""
    w75 = 0
    w50 = 0
    b50 = 0
    b40 = 0
    nw50 = 0
    nw40 = 0
    l50 = 0
    l40 = 0
    # print(list(solution_candidate))
    for row in solution_candidate.index:
        whitepop = solution_candidate.loc[row, 'white_population'+column_suffix]
        blackpop = solution_candidate.loc[row, 'black_population'+column_suffix]
        latinopop = solution_candidate.loc[row, 'latino_population'+column_suffix]
        otherpop = solution_candidate.loc[row, 'other_population'+column_suffix]
        totalpop = solution_candidate.loc[row, population_column+column_suffix]
        bpercent = blackpop / totalpop
        wpercent = whitepop / totalpop
        lpercent = latinopop / totalpop
        nwpercent = (blackpop + latinopop + otherpop) / totalpop
        if bpercent >= 0.4:
            b40 += 1
        if bpercent >= 0.5:
            b50 += 1
        if wpercent >= 0.75:
            w75 += 1
        if wpercent >= 0.5:
            w50 += 1
        if nwpercent >= 0.5:
            nw50 += 1
        if nwpercent >= 0.4:
            nw40 += 1
        if lpercent >= 0.5:
            l50 += 1
        if lpercent >= 0.4:
            l40 += 1

    return w75, w50, b50, b40, nw50, nw40, l50, l40


def calculate_minority_majority_districts(solution_candidate, column_suffix=None, population_column='total_population'):
    districtSums = solution_candidate.groupby(['District']).sum()

    w75, w50, b50, b40, nw50, nw40, l50, l40 = getMMstats(districtSums, column_suffix, population_column)

    return w75, w50, b50, b40, nw50, nw40, l50, l40


def calculate_metrics(solution_candidate, congress_districts, population_column,
                                     include_num_county_split='False', include_minority_majority='False'):
    """
    Calculates various solution metrics and adds them to the existing costs list
    :param solution_candidate: the solution candidate for which costs should be calculated
    # :param costs: a list of costs
    :return: the updated list of costs
    """

    # due to an update in geopandas, we now use a different crs to calculate compactness
    # solution_candidate = solution_candidate.to_crs(epsg=3035)

    # first, we aggregate the solution candidate by district
    district_df_agg = solution_candidate.dissolve(by='District', aggfunc='sum')
    # state_df_agg = solution_candidate.dissolve(by='state', aggfunc='sum')
    # state_area = state_df_agg['geometry'].area[0]

    # also calculate the state perimeter
    territory_id = [1] * solution_candidate.shape[0]

    # calculate the state perimeter
    data_table_perimeter = copy.deepcopy(solution_candidate)
    data_table_perimeter['territory'] = territory_id
    data_table_agg = data_table_perimeter.dissolve(by='territory')
    state_perimeter = data_table_agg.length.iloc[0]

    # new: convert to multipolygon for easier calculation of compactness measures
    multipolygons = []
    for geo in district_df_agg.geometry:
        if isinstance(geo, Polygon):
            multipolygons.append(MultiPolygon([geo]))
        elif isinstance(geo, MultiPolygon):
            multipolygons.append(geo)
        else:
            ValueError('geometry should be either Polygon or Multipolygon')

    #         print(multipolygons)
    district_df_agg['geometry'] = gp.GeoSeries(multipolygons, index=district_df_agg.geometry.index)

    # then, we can leverage this dissolved dataframe in the subsequent calculations without need for re-calculation

    print('calculate compactness')
    temp_compactness = calculate_compactness_boundary_lengths(district_df_agg=district_df_agg,
                                                              state_perimeter=state_perimeter)

    print('calculate compactness radial')
    temp_compactness_radial = calculate_compactness_radial(district_df=solution_candidate,
                                                           district_df_agg=district_df_agg)

    print('calculate compactness measures')
    temp_compactness_measures = calculate_compactness_measures(district_df_agg=district_df_agg)

    print('calculate population equality')
    temp_equality = calculate_population_equality(district_df=solution_candidate, population_column=population_column)

    # print('calculate socioeconomic homogeneity')
    # temp_homogeneity = calculate_socioeconomic_homogeneity(district_df=solution_candidate)

    # print('calculate similarity')
    # temp_similarity = calculate_similarity_to_existing_districting(district_df_agg=district_df_agg,
    #                                                                congress_districts=congress_districts,
    #                                                                state_area=state_area)

    print('calculate lower and upper tolerance')
    lower, upper = calculate_population_tolerance(district_df=solution_candidate, population_column=population_column)

    result = [temp_compactness, temp_compactness_radial]
    result.extend(temp_compactness_measures)
    result.extend([temp_equality, lower, upper])

    if include_num_county_split == "True":
        print('calculate num county split')
        temp_num_county_split = calculate_num_counties_split(district_df=solution_candidate)
        result.append(temp_num_county_split)
    if include_minority_majority == "True":
        print('calculate majority minority districts')
        result.extend(calculate_minority_majority_districts(solution_candidate,column_suffix=None,
                                                            population_column=population_column))

    return result

