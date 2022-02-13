
using NCDatasets, DataStructures

gcms = ["CESM2-WACCM","GFDL-ESM4"]

paths = ["/home/milo/PCIC/Maestría/2doSemestre/seminario/github/data/cmip6/CESM2-WACCM/historical/co2/",
        "/home/milo/PCIC/Maestría/2doSemestre/seminario/github/data/cmip6/GFDL-ESM4/historical/co2/"]

netcdfs = ["co2_Amon_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412.nc","co2_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_185001-194912.nc"]

cd(paths[2])
printstyled("Procesando datos de CO2 del GCM"*gcms[2]*"\n"; color = :yellow)

ds_historical = NCDataset(paths[2]*netcdfs[2],"r");
inicio = 1
fin = 12

for i in 1:100

   filename = netcdfs[2][1:43]*string(1849+i)*".nc"

   println("Procesando NetCFD: "*filename)

   ds = NCDataset(filename,"c", attrib = OrderedDict(
       "external_variables"        => "areacella",
       "history"                   => "File was processed by fremetar (GFDL analog of CMOR). TripleID: [exper_id_MFLg3OOf97,realiz_id_6UiFuoEKMa,run_id_PhuSv75why]",
       "table_id"                  => "Amon",
       "activity_id"               => "CMIP",
       "branch_method"             => "standard",
       "branch_time_in_child"      => 0.0,
       "branch_time_in_parent"     => 36500.0,
       "comment"                   => "<null ref>",
       "contact"                   => "gfdl.climate.model.info@noaa.gov",
       "Conventions"               => "CF-1.7 CMIP-6.0 UGRID-1.0",
       "creation_date"             => "2019-08-06T01:24:53Z",
       "data_specs_version"        => "01.00.27",
       "experiment"                => "all-forcing simulation of the recent past",
       "experiment_id"             => "historical",
       "forcing_index"             => Int32(1),
       "frequency"                 => "mon",
       "further_info_url"          => "https://furtherinfo.es-doc.org/CMIP6.NOAA-GFDL.GFDL-ESM4.historical.none.r1i1p1f1",
       "grid"                      => "atmos data regridded from Cubed-sphere (c96) to 180,288; interpolation method: conserve_order2",
       "grid_label"                => "gr1",
       "initialization_index"      => Int32(1),
       "institution"               => "National Oceanic and Atmospheric Administration, Geophysical Fluid Dynamics Laboratory, Princeton, NJ 08540, USA",
       "institution_id"            => "NOAA-GFDL",
       "license"                   => "CMIP6 model data produced by NOAA-GFDL is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License (https://creativecommons.org/licenses/). Consult https://pcmdi.llnl.gov/CMIP6/TermsOfUse for terms of use governing CMIP6 output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further_info_url (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law.",
       "mip_era"                   => "CMIP6",
       "nominal_resolution"        => "100 km",
       "parent_activity_id"        => "CMIP",
       "parent_experiment_id"      => "piControl",
       "parent_mip_era"            => "CMIP6",
       "parent_source_id"          => "GFDL-ESM4",
       "parent_time_units"         => "days since 0001-1-1",
       "parent_variant_label"      => "r1i1p1f1",
       "physics_index"             => Int32(1),
       "product"                   => "model-output",
       "realization_index"         => Int32(1),
       "realm"                     => "atmos",
       "source"                    => "GFDL-ESM4 (2018):
   atmos: GFDL-AM4.1 (Cubed-sphere (c96) - 1 degree nominal horizontal resolution; 360 x 180 longitude/latitude; 49 levels; top level 1 Pa)
   ocean: GFDL-OM4p5 (GFDL-MOM6, tripolar - nominal 0.5 deg; 720 x 576 longitude/latitude; 75 levels; top grid cell 0-2 m)
   seaIce: GFDL-SIM4p5 (GFDL-SIS2.0, tripolar - nominal 0.5 deg; 720 x 576 longitude/latitude; 5 layers; 5 thickness categories)
   land: GFDL-LM4.1
   aerosol: interactive
   atmosChem: GFDL-ATMCHEM4.1 (full atmospheric chemistry)
   ocnBgchem: GFDL-COBALTv2
   landIce: GFDL-LM4.1
   (GFDL ID: 2019_0353)",
       "source_id"                 => "GFDL-ESM4",
       "source_type"               => "AOGCM AER CHEM BGC",
       "sub_experiment"            => "none",
       "sub_experiment_id"         => "none",
       "title"                     => "NOAA GFDL GFDL-ESM4 model output prepared for CMIP6 all-forcing simulation of the recent past",
       "tracking_id"               => "hdl:21.14100/e2089252-64eb-4588-969c-a506012e4224",
       "variable_id"               => "co2",
       "variant_info"              => "N/A",
       "references"                => "see further_info_url attribute",
       "variant_label"             => "r1i1p1f1",
   ))

   # Dimensions

   ds.dim["bnds"] = 2
   ds.dim["time"] = Inf # unlimited dimension
   ds.dim["plev"] = 19
   ds.dim["lat"] = 180
   ds.dim["lon"] = 288

   # Declare variables

   ncbnds = defVar(ds,"bnds", Float64, ("bnds",), attrib = OrderedDict(
       "long_name"                 => "vertex number",
   ))

   ncco2 = defVar(ds,"co2", Float32, ("lon", "lat", "plev", "time"), attrib = OrderedDict(
       "long_name"                 => "Mole Fraction of CO2",
       "units"                     => "mol mol-1",
       "missing_value"             => Float32(1.0e20),
       "_FillValue"                => Float32(1.0e20),
       "cell_methods"              => "time: mean",
       "cell_measures"             => "area: areacella",
       "standard_name"             => "mole_fraction_of_carbon_dioxide_in_air",
       "interp_method"             => "conserve_order2",
       "original_name"             => "co2",
   ))

   nclat = defVar(ds,"lat", Float64, ("lat",), attrib = OrderedDict(
       "long_name"                 => "latitude",
       "units"                     => "degrees_north",
       "axis"                      => "Y",
       "bounds"                    => "lat_bnds",
       "standard_name"             => "latitude",
       "cell_methods"              => "time: point",
   ))

   nclat_bnds = defVar(ds,"lat_bnds", Float64, ("bnds", "lat"), attrib = OrderedDict(
       "long_name"                 => "latitude bounds",
       "units"                     => "degrees_north",
       "axis"                      => "Y",
   ))

   nclon = defVar(ds,"lon", Float64, ("lon",), attrib = OrderedDict(
       "long_name"                 => "longitude",
       "units"                     => "degrees_east",
       "axis"                      => "X",
       "bounds"                    => "lon_bnds",
       "standard_name"             => "longitude",
       "cell_methods"              => "time: point",
   ))

   nclon_bnds = defVar(ds,"lon_bnds", Float64, ("bnds", "lon"), attrib = OrderedDict(
       "long_name"                 => "longitude bounds",
       "units"                     => "degrees_east",
       "axis"                      => "X",
   ))

   ncplev = defVar(ds,"plev", Float64, ("plev",), attrib = OrderedDict(
       "long_name"                 => "pressure",
       "units"                     => "Pa",
       "axis"                      => "Z",
       "positive"                  => "down",
       "standard_name"             => "air_pressure",
       "description"               => "There are 19  requested levels. If the upper levels are above the model top they should be omitted.",
   ))

   nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
       "long_name"                 => "time",
       "units"                     => "days since 1850-01-01 00:00:00",
       "axis"                      => "T",
       "calendar_type"             => "noleap",
       "calendar"                  => "noleap",
       "bounds"                    => "time_bnds",
       "standard_name"             => "time",
       "description"               => "Temporal mean",
   ))

   nctime_bnds = defVar(ds,"time_bnds", Float64, ("bnds", "time"), attrib = OrderedDict(
       "long_name"                 => "time axis boundaries",
       "units"                     => "days since 1850-01-01 00:00:00",
   ))


    # Define variables

    ncco2[:] = ds_historical["co2"][1:288,1:180,1:19,inicio:fin];

    nclat[:] = ds_historical["lat"][:];
    nclon[:] = ds_historical["lon"][:];
    ncplev[:] = ds_historical["plev"][:];
    nctime[:] = ds_historical["time"][:][inicio:fin];
    nctime_bnds[:] = ds_historical["time_bnds"][:][1:2,inicio:fin];
    nclat_bnds[:] = ds_historical["lat_bnds"][:];
    nclon_bnds[:] = ds_historical["lon_bnds"][:];

    close(ds)

    inicio = inicio + 12
    fin = fin + 12

end


# https://alexander-barth.github.io/NCDatasets.jl/stable/
# https://github.com/Alexander-Barth/NCDatasets.jl
# https://docs.juliahub.com/NCDatasets/lxvtD/0.10.0/
