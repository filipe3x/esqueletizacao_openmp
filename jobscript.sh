#!/bin/bash
#
kernel=(ske5.30static) 
#kernel=(ske5.30static ske5.30dynamic ske6.50static ske6.50dynamic)

input=(horse 256Kcircle) 
#input=(horse 256Kcircle 512Kcircle 1Mcircle 4Mcircle 16Mcircle 32Mcircle)

output=skeleton.pgm

folder=ppmimages

errorlog=joberrors.log

maxthreads=$1

RESULTS=()


function calc_kbest() {
  array=("$@")

  a=${array[0]}
  b=${array[1]}
  c=${array[2]}

  sum=$((a + b + c))

  echo $(($sum / 3))
}

function calc_3_lower() {
  array=("$@")
  echo ${array[@]} | tr ' ' '\n' | sort -n | head -n3 | tr '\n' ' '
}

for k in "${kernel[@]}"
do
  echo "** $k kernel **"
  for i in "${input[@]}"
  do
    i=$folder/$i.pgm
    echo "input -> $i"
    echo "threads = 1"
    RESULTS=()
    for r in {1..3} ## number of runs for sequential
    do
      comm="./$k $i $folder/$output 1"
      echo try= $r comm= $comm
      time=$($comm 2>>$errorlog)
      echo - $time -
      RESULTS+=($time)
      sleep 2
    done
    KBEST=($(calc_3_lower "${RESULTS[@]}"))
    echo k-best score = $(calc_kbest "${KBEST[@]}")
    for t in {2..4} ## number of total threads here
    do
      echo "threads = $t"
      RESULTS=()
      for r in {1..3} ## number of runs for each parallel configuration
      do
        comm="./$k $i $folder/$output 0 $t"
        echo try= $r comm= $comm
        time=$($comm 2>>$errorlog)
        echo - $time -
        RESULTS+=($time)
        sleep 2
      done
      KBEST=($(calc_3_lower "${RESULTS[@]}"))
      echo k-best score = $(calc_kbest "${KBEST[@]}")
    done
  done
done

exit 0;

