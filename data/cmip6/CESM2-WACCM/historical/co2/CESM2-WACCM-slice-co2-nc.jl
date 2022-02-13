
using NCDatasets, DataStructures

gcms = ["CESM2-WACCM","GFDL-ESM4"]

paths = ["/home/milo/PCIC/Maestría/2doSemestre/seminario/github/data/cmip6/CESM2-WACCM/historical/co2/",
        "/home/milo/PCIC/Maestría/2doSemestre/seminario/github/data/cmip6/GFDL-ESM4/historical/co2/"]

netcdfs = ["co2_Amon_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412.nc","co2_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_185001-194912.nc"]

cd(paths[1])
printstyled("Procesando datos de CO2 del GCM"*gcms[1]*"\n"; color = :yellow)

ds_historical = NCDataset(paths[1]*netcdfs[1],"r");
inicio = 1
fin = 12

for i in 1:165

   filename = netcdfs[1][1:44]*string(1849+i)*".nc"

   println("Procesando NetCFD: "*filename)

    ds = NCDataset(filename,"c", attrib = OrderedDict(
        "Conventions"               => "CF-1.7 CMIP-6.2",
        "activity_id"               => "CMIP",
        "case_id"                   => "4",
        "cesm_casename"             => "b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.001",
        "contact"                   => "cesm_cmip6@ucar.edu",
        "creation_date"             => "2019-01-30T21:27:29Z",
        "data_specs_version"        => "01.00.29",
        "experiment"                => "all-forcing simulation of the recent past",
        "experiment_id"             => "historical",
        "external_variables"        => "areacella",
        "forcing_index"             => Int64(1),
        "frequency"                 => "mon",
        "further_info_url"          => "https://furtherinfo.es-doc.org/CMIP6.NCAR.CESM2-WACCM.historical.none.r1i1p1f1",
        "grid"                      => "native 0.9x1.25 finite volume grid (192x288 latxlon)",
        "grid_label"                => "gn",
        "initialization_index"      => Int64(1),
        "institution"               => "National Center for Atmospheric Research, Climate and Global Dynamics Laboratory, 1850 Table Mesa Drive, Boulder, CO 80305, USA",
        "institution_id"            => "NCAR",
        "license"                   => "CMIP6 model data produced by <The National Center for Atmospheric Research> is licensed under a Creative Commons Attribution-[]ShareAlike 4.0 International License (https://creativecommons.org/licenses/). Consult https://pcmdi.llnl.gov/CMIP6/TermsOfUse for terms of use governing CMIP6 output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further_info_url (recorded as a global attribute in this file)[]. The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law.",
        "mip_era"                   => "CMIP6",
        "model_doi_url"             => "https://doi.org/10.5065/D67H1H0V",
        "nominal_resolution"        => "100 km",
        "parent_activity_id"        => "CMIP",
        "parent_experiment_id"      => "piControl",
        "parent_mip_era"            => "CMIP6",
        "parent_source_id"          => "CESM2-WACCM",
        "parent_time_units"         => "days since 0001-01-01 00:00:00",
        "parent_variant_label"      => "r1i1p1f1",
        "physics_index"             => Int64(1),
        "product"                   => "model-output",
        "realization_index"         => Int64(1),
        "realm"                     => "atmos",
        "source"                    => "CESM2 (2017): atmosphere: CAM6 (0.9x1.25 finite volume grid; 288 x 192 longitude/latitude; 70 levels; top level 4.5e-6 mb); ocean: POP2 (320x384 longitude/latitude; 60 levels; top grid cell 0-10 m); sea_ice: CICE5.1 (same grid as ocean); land: CLM5 0.9x1.25 finite volume grid; 288 x 192 longitude/latitude; 70 levels; top level 4.5e-6 mb); aerosol: MAM4 (0.9x1.25 finite volume grid; 288 x 192 longitude/latitude; 70 levels; top level 4.5e-6 mb); atmosChem: WACCM (0.9x1.25 finite volume grid; 288 x 192 longitude/latitude; 70 levels; top level 4.5e-6 mb; landIce: CISM2.1; ocnBgchem: MARBL (320x384 longitude/latitude; 60 levels; top grid cell 0-10 m)",
        "source_id"                 => "CESM2-WACCM",
        "source_type"               => "AOGCM BGC CHEM AER",
        "sub_experiment"            => "none",
        "sub_experiment_id"         => "none",
        "table_id"                  => "Amon",
        "tracking_id"               => "hdl:21.14100/6dfdee4d-d113-4701-9525-2835dab596f8",
        "variable_id"               => "co2",
        "variant_info"              => "CMIP6 CESM2 hindcast (1850-2014) with high-top atmosphere (WACCM6) with interactive chemistry (TSMLT1), interactive land (CLM5), coupled ocean (POP2) with biogeochemistry (MARBL), interactive sea ice (CICE5.1), and non-evolving land ice (CISM2.1)",
        "variant_label"             => "r1i1p1f1",
        "branch_time_in_parent"     => 20075.0,
        "branch_time_in_child"      => 674885.0,
        "branch_method"             => "standard",
    ))

    # Dimensions

    ds.dim["lat"] = 192
    ds.dim["plev"] = 19
    ds.dim["nbnd"] = 2
    ds.dim["lon"] = 288
    ds.dim["time"] = Inf # unlimited dimension

    # Declare variables

    ncco2 = defVar(ds,"co2", Float32, ("lon", "lat", "plev", "time"), attrib = OrderedDict(
        "_FillValue"                => Float32(1.0e20),
        "cell_measures"             => "area: areacella",
        "cell_methods"              => "time: mean",
        "comment"                   => "Mole fraction is used in the construction mole_fraction_of_X_in_Y, where X is a material constituent of Y.",
        "coordinates"               => "time plev lat lon",
        "description"               => "Mole fraction is used in the construction mole_fraction_of_X_in_Y, where X is a material constituent of Y.",
        "frequency"                 => "mon",
        "id"                        => "co2",
        "long_name"                 => "Mole Fraction of CO2",
        "mipTable"                  => "Amon",
        "missing_value"             => 1.0e20,
        "out_name"                  => "co2",
        "prov"                      => "Amon ((isd.003))",
        "realm"                     => "atmos",
        "standard_name"             => "mole_fraction_of_carbon_dioxide_in_air",
        "time"                      => "time",
        "time_label"                => "time-mean",
        "time_title"                => "Temporal mean",
        "title"                     => "Mole Fraction of CO2",
        "type"                      => "real",
        "units"                     => "mol mol-1",
        "variable_id"               => "co2",
    ))

    nclat = defVar(ds,"lat", Float64, ("lat",), attrib = OrderedDict(
        "axis"                      => "Y",
        "bounds"                    => "lat_bnds",
        "standard_name"             => "latitude",
        "title"                     => "Latitude",
        "type"                      => "double",
        "units"                     => "degrees_north",
        "valid_max"                 => 90.0,
        "valid_min"                 => -90.0,
    ))

    nclon = defVar(ds,"lon", Float64, ("lon",), attrib = OrderedDict(
        "axis"                      => "X",
        "bounds"                    => "lon_bnds",
        "standard_name"             => "longitude",
        "title"                     => "Longitude",
        "type"                      => "double",
        "units"                     => "degrees_east",
        "valid_max"                 => 360.0,
        "valid_min"                 => 0.0,
    ))

    ncplev = defVar(ds,"plev", Float64, ("plev",), attrib = OrderedDict(
        "axis"                      => "Z",
        "positive"                  => "down",
        "requested"                 => "100000. 92500. 85000. 70000. 60000. 50000. 40000. 30000. 25000. 20000. 15000. 10000. 7000. 5000. 3000. 2000. 1000. 500. 100.",
        "standard_name"             => "air_pressure",
        "title"                     => "pressure",
        "type"                      => "double",
        "units"                     => "Pa",
    ))

    nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
        "axis"                      => "T",
        "bounds"                    => "time_bnds",
        "standard_name"             => "time",
        "title"                     => "time",
        "type"                      => "double",
        "units"                     => "days since 0001-01-01 00:00:00",
        "calendar"                  => "noleap",
    ))

    nctime_bnds = defVar(ds,"time_bnds", Float64, ("nbnd", "time"), attrib = OrderedDict(
        "calendar"                  => "noleap",
        "units"                     => "days since 0001-01-01 00:00:00",
    ))

    nclat_bnds = defVar(ds,"lat_bnds", Float32, ("nbnd", "lat"), attrib = OrderedDict(
        "units"                     => "degrees_north",
    ))

    nclon_bnds = defVar(ds,"lon_bnds", Float32, ("nbnd", "lon"), attrib = OrderedDict(
        "units"                     => "degrees_east",
    ))


    # Define variables

    ncco2[:] = ds_historical["co2"][1:288,1:192,1:19,inicio:fin];

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
