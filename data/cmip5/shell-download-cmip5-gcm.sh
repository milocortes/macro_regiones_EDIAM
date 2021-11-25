#!/usr/bin/env bash
filename='gcm-models.txt'

# Declaramos un arreglo de los rcps
declare -a rcps=("historical" "rcp26" "rcp45" "rcp85")

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
