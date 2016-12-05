#!/bin/bash

#matrix_N=(500 1000 5000 10000 50000)
#process=(2 4 8 12 16) 
#begin_omp=(500 1000)

matrix_N=(500)
process=(2) 
begin_omp=(500)
mkdir result
cd result/

for i in "${matrix_N[@]}" 
do
  mkdir N${i}
  for j in "${process[@]}"
  do 
    mkdir N${i}/process${j}
    for k in "${begin_omp[@]}"
    do
    mkdir N${i}/process${j}/beginomp${j}

    cd N${i}/process${j}/beginomp${j}
    cp ../../../../* ./
    make clean
    sed -i '6s/5000/'"${i}"'/' Makefile
    sed -i '21s/12/'"${j}"'/' Makefile
    sed -i '9s/1000/'"${k}"'/' Makefile
    make
    qsub run_hybrid.simple
    cd ../../../
    done
  done
done

