from pyesgf.search import SearchConnection
import os
import pandas as pd
import requests
from tqdm import tqdm



cmip_url = ["https://esgf-node.llnl.gov/esg-search",
            "https://esgf-node.ipsl.upmc.fr/esg-search",
            "https://esgf-data.dkrz.de/esg-search",
            "https://esgf-index1.ceda.ac.uk/esg-search"]
cmip_url = ["https://esgf-node.llnl.gov/esg-search"]


gcm_models = ['GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','MIROC-ESM','MIROC-ESM-CHEM','MIROC5',
                'MPI-ESM-LR','MPI-ESM-MR','MRI-ESM1','NorESM1-M','NorESM1-ME']
rcps = ['historical','rcp26','rcp45','rcp85']

dic_connections = {}
for cmip in cmip_url:
    print("Búsqueda en la liga {}".format(cmip))
    dic_connections[cmip]= SearchConnection(cmip, distrib=True)

    for gcm in gcm_models:
        for rcp in rcps:
            results = dic_connections[cmip].new_context(
                latest=True,
                facets='null',
                project='CMIP5',
                model=gcm,
                experiment=rcp,
                time_frequency="mon",
                realm='atmos',
                ensemble='r1i1p1').search()


            if len(results) >0:
                # We can use an anonymous “lambda” function to extract the filename and URL of each result.
                hit = results[0].file_context().search()
                files = map(lambda f : {'filename': f.filename, 'url': f.download_url}, hit)
                files = list(files)
                files = list(filter(lambda x: ('co2mass_' in x['filename'] or 'co2_' in x['filename'] or 'tas_' in x['filename']), files))

                verifica = []
                for f in files:
                    for var in ['co2mass_','co2_','tas_']:
                        if var in f["filename"]:
                            verifica.append(var)
                if len(set(verifica))==3:
                    print("Liga: {}\nModelo: {}\nRCP{}".format(cmip,gcm,rcp))
                    f = open("{}_{}.txt".format(gcm,rcp),"x")
                    for file in files:
                        f.write("{}\n".format(file["url"]))
                        print(file["url"])
                    f.close()
####################################
### %%%%%%%%%%%%%% Descarga CMIP y xarray
# https://claut.gitlab.io/man_ccia/lab2.html
# https://esgf-pyclient.readthedocs.io/en/stable/
# https://www.indierocks.mx/musica/noticias/unete-a-fuerzazopi-para-apoyar-a-gerardo-pimentel/
# http://xarray.pydata.org/en/stable/user-guide/dask.html
# https://esgf-pyclient.readthedocs.io/en/latest/notebooks/demo/subset-cmip6.html
# https://esgf-pyclient.readthedocs.io/en/latest/notebooks/demo/subset-cmip5.html
# http://gallery.pangeo.io/repos/pangeo-data/pangeo-tutorial-gallery/xarray.html
# http://gallery.pangeo.io/repos/pangeo-data/pangeo-tutorial-gallery/geopandas.html
# http://gallery.pangeo.io/repos/pangeo-gallery/cmip6/index.html
# http://gallery.pangeo.io/repos/pangeo-gallery/cmip6/basic_search_and_load.html
# https://kpegion.github.io/Pangeo-at-AOES/examples/advanced-analysis.html
# https://readthedocs.org/projects/esgf-pyclient/downloads/pdf/latest/
# https://nci-data-training.readthedocs.io/en/latest/_notebook/climate/1_04_Xarray_calculate_metrics_CMIP6.html
# https://nordicesmhub.github.io/Norway_Sweden_training/pangeo/CMIP6_example.html
# https://github.com/ESGF/sproket
# https://esgf-node.llnl.gov/projects/cmip6/
# http://xarray.pydata.org/en/stable/user-guide/indexing.html
# http://xarray.pydata.org/en/stable/user-guide/plotting.html
# https://scitools.org.uk/cartopy/docs/v0.13/matplotlib/gridliner.html
# https://jeffpollock9.github.io/maximum-likelihood-estimation-with-tensorflow-probability-and-pystan/
#
#
