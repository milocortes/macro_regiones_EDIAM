from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import matplotlib.animation as animation
import re
import os


# Cambiamos de directorio
os.chdir("../../data/WC_nc")
# Generamos una lista de los directorios dentro de WC_nc
dir_data = os.listdir()
map_dir_data = {dir:var for dir,var in zip(dir_data,['tmax','prec','tmin'])}
dict_dir = {}

anioInicio = 1970
anioFinal = 2000

rp = r"[0-9]{4,7}-[0-9]{2,2}"

## Guardamos las matrices de los netcdf en el diccionario dict_dir
for dir in dir_data:
    print("Leyendo netCDF del directorio {}".format(dir))
    os.chdir(dir)
    files = os.listdir()
    files = sorted(files)
    files = [file for file in files if int(file[16:len(file)-6])>= anioInicio and int(file[16:len(file)-6])<=anioFinal]
    dict_dir[dir]=[]
    for file in files:
        fecha = re.findall(rp,file).pop()
        print(fecha)
        data = Dataset(file,mode = 'r')
        dict_dir[dir].append(data[map_dir_data[dir]][:])
        data.close()
    os.chdir("../")

dict_dir['WC_Tmax'] = np.array(dict_dir['WC_Tmax'])
dict_dir['WC_Tmin'] = np.array(dict_dir['WC_Tmin'])
dict_dir['WC_Prec'] = np.array(dict_dir['WC_Prec'])

## Calculamos las anomalías
# Valores medios
mean_tmax = np.array([np.nanmean(x) for x in dict_dir['WC_Tmax']])
mean_tmin = np.array([np.nanmean(x) for x in dict_dir['WC_Tmin']])
mean_prec = np.array([np.nanmean(x) for x in dict_dir['WC_Prec']])

# Temperatura media
tmean = (mean_tmax + mean_tmin)/2

# Temperaturas medias por mes
#periodos = ["{}".format(x)+"-{0:0=2d}".format(y) for x in range(1970,2001) for y in range(1,13)]
anios = np.array([anio for anio in range(1970,2001)])
meses = np.array([mes for mes in range(1,13)] * len(anios))

mes_mean_tmax = []
mes_mean_tmin = []
mes_mean_prec = []
mes_mean_tmean = []

for mes in range(1,13):
    indices = np.argwhere(meses==mes).flatten()
    mes_mean_tmax.append(np.nanmean(dict_dir['WC_Tmax'][indices]))
    mes_mean_tmin.append(np.nanmean(dict_dir['WC_Tmin'][indices]))
    mes_mean_prec.append(np.nanmean(dict_dir['WC_Prec'][indices]))
    mes_mean_tmean.append(np.nanmean(tmean[indices]))

# Anomalías
size = (anioFinal-anioInicio+1)*12
anoma_tmax = np.empty(size)
anoma_tmin = np.empty(size)
anoma_prec = np.empty(size)
anoma_tmean = np.empty(size)

for mes in range(1,13):
    indices = np.argwhere(meses==mes).flatten()
    anoma_tmax[indices] = mean_tmax[indices] - mes_mean_tmax[mes-1]
    anoma_tmin[indices] = mean_tmin[indices] - mes_mean_tmin[mes-1]
    anoma_prec[indices] = mean_prec[indices] - mes_mean_prec[mes-1]
    anoma_tmean[indices] = tmean[indices] - mes_mean_tmean[mes-1]

# Calculamos las matrices climáticas
Tmean = (dict_dir['WC_Tmax'] + dict_dir['WC_Tmin']) /2

fig = plt.figure()
ax = fig.add_subplot(111)
ims=[]

for i in range(1,13):
    indices  = np.argwhere(meses==i).flatten()
    ttl = plt.text(0.5, 1.01, "Temperatura media climatológica 1970-200 {0:0=2d}".format(i), horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes)
    ims.append([plt.imshow(Tmean[indices].mean(0), cmap='jet', vmin=-10, vmax=25,interpolation='none',origin='lower', animated=True),ttl])
    #plt.cla()


ani = animation.ArtistAnimation(fig, ims, interval=500, blit=False,repeat_delay=0)
plt.colorbar()
plt.show()
