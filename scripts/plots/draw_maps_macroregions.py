##### Concat historical an rcp TAS files
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

historical_name = 'tas_Amon_MRI-ESM1_historical_r1i1p1_185101-200512.nc'
rcp_name = 'tas_Amon_MRI-ESM1_rcp85_r1i1p1_200601-210012.nc'
historical = xr.open_dataset(historical_name, chunks={'time': 120})
rcp = xr.open_dataset(rcp_name, chunks={'time': 120})

historical_tas = historical['tas']
rcp_tas = rcp['tas']
total_tas = xr.concat([historical_tas,rcp_tas],dim='addHistorical')
total_tas_mean = total_tas.mean(dim='addHistorical')

## Group by year
total_tas_mean_c = total_tas_mean - 273.15

total_tas_mean_c_year = total_tas_mean_c.groupby('time.year').mean()


# Global plot
fig = plt.figure(1, figsize=[30,13])
# Set the projection to use for plotting
ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.coastlines()
gridlines = ax.gridlines(draw_labels=False, linewidth=2, color='black', alpha=1.0, linestyle='--')
gridlines.xlocator = mticker.FixedLocator([-30, 60])
gridlines.ylines = False
gridlines.xlabels_top = True
gridlines.xlabels_bottom = True
# Pass ax as an argument when plotting. Here we assume data is in the same coordinate reference system than the projection chosen for plotting
# isel allows to select by indices instead of the time values
total_tas_mean_c_year.mean(dim="year").plot.pcolormesh(ax=ax, cmap='coolwarm')
plt.title('Macro-Regiones')
plt.show()


# America
america = total_tas_mean_c_year.mean(dim="year").isel(lon=slice(170, 300))
fig = plt.figure(1, figsize=[30,13])
# Set the projection to use for plotting
ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.coastlines()
america.plot.pcolormesh(ax=ax, cmap='coolwarm')
plt.show()

# Europa-Africa
eurafrica =  total_tas_mean_c_year.mean(dim="year").where((total_tas_mean_c_year.lon <55) | (total_tas_mean_c_year.lon >330))
fig = plt.figure(1, figsize=[30,13])
# Set the projection to use for plotting
ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.coastlines()
eurafrica_uno.plot.pcolormesh(ax=ax, cmap='coolwarm')
plt.show()

# Asia
asia = total_tas_mean_c_year.mean(dim="year").isel(lon=slice(55, 170))
fig = plt.figure(1, figsize=[30,13])
# Set the projection to use for plotting
ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.coastlines()
asia.plot.pcolormesh(ax=ax, cmap='coolwarm')
plt.show()
