from pyesgf.search import SearchConnection
import os
import pandas as pd
import requests
from tqdm import tqdm

conn = SearchConnection('https://esgf-node.llnl.gov/esg-search', distrib=True)

query = conn.new_context(
    latest=True,
    facets='null',
    project='CMIP5',
    model='MRI-ESM1',
    experiment='rcp85',
    variable='co2mass',
    time_frequency="mon",
    realm='atmos',
    ensemble='r1i1p1')

results = query.search()

# We can use an anonymous “lambda” function to extract the filename and URL of each result.
hit = results[0].file_context().search()
files = map(lambda f : {'filename': f.filename, 'url': f.download_url}, hit)
files = list(files)
files = list(filter(lambda x: ('co2mass_' in x['filename'] or 'co2_' in x['filename'] or 'tas_' in x['filename']) and '200601' in x['filename'], files))
files


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import os



co2 = Dataset("co2_Amon_MRI-ESM1_rcp85_r1i1p1_200601-201512.nc",mode = "r")
np_co2 = np.array(co2["co2"])
np_co2[0].mean(axis=0)


tas = Dataset("tas_Amon_MRI-ESM1_rcp85_r1i1p1_200601-210012.nc",mode="r")
np_tas = np.array(tas["tas"])
