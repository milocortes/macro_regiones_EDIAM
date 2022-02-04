#!/usr/bin/env bash
filename='gcm-models.txt'

# Declaramos un arreglo de los rcps
declare -a rcps=("historical" "ssp126" "ssp245" "ssp370" "ssp585")

# Declaramos un arreglo de las variables climáticas
declare -a climvariables=("co2mass" "co2" "tas")

# Creamos los directorios necesarios
echo "CREAMOS LOS DIRECTORIOS NECESARIOS"
while read model; do

  if [ -d "$model" ]; then
    rm -rf "$model"
  fi

  mkdir "$model"
  cd "$model"
  for rcp in ${rcps[@]}; do

    mkdir "$rcp"
    cd "$rcp"
    for climvar in ${climvariables[@]}; do
      mkdir "$climvar"
    done
    cd ..
  done
  cd ..

done < $filename

while read model; do

  for rcp in ${rcps[@]}; do
    model_download=$model"_"$rcp".txt"

    if [ -f "$model_download" ]; then
    echo "$model_download exists."
    wget -i $model_download
    echo "***********************************"
    echo "MOVEMOS LOS NCDF AL DIRECTORIO"
    echo "./$model/$rcp/"
    echo "***********************************"

    mv *.nc "./$model/$rcp/"

    fi

  done

done < $filename
