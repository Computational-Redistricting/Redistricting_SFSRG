# import the census api communication packages
from cenpy import products

acs = products.ACS(2017)
census = products.Decennial2010()


# specify the state
state = "Nebraska"

# get the acs and census data for the specified state from ACS and Census
# for the ACS data, we restrict ourselves to certain population columns for now (can add income as well)
acs_state = acs.from_state(state, level='tract', variables=['B03002_001E','B03002_003E','B03002_004E','B03002_012E','B07011_001E'])