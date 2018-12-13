#!/bin/sh
#

kernel=(skeletonize) 
#kernel=(ske5.30static ske5.30dynamic ske6.50static ske6.50dynamic)

input=(horse 256Kcircle) 
#input=(horse 256Kcircle 1Mcircle 2Mcircle 4Mcircle 16Mcircle 32Mcircle)

commseq="./$k $i $folder/$output 3"
commpar="./$k $i $folder/$output 2 $t"
commmpi="mpirun -np $t ./$k $i $folder/$output 999"
commmpi="mpirun -bynode -mca btl ^openib -np $t ./$k $i $folder/$output 999"
commmpi="mpirun --report-bindings --map-by node -mca btl ^openib -np $t ./$k $i $folder/$output 999"
commmpi="mpirun --map-by node -mca btl ^openib,mx,sm,self -np $t ./$k $i $folder/$output 999"

output=skeleton.pgm

folder=ppmimages

errorlog=results/joberrors.log

OUTPUTGRAPH=results/resultsgraph$(date +%d-%m-%Y-%H%M%S).log

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

calc() {
  awk "BEGIN { print "$*" }"
}

for k in "${kernel[@]}"
do
#  if [ $k == "skeletonize" ]; then module load gcc/5.3.0 2>>$errorlog; fi ## if we need specific modules for this kernel
  echo "** $k kernel | $(date)" | tee -a $OUTPUTGRAPH
  for i in "${input[@]}"
  do
    i=$folder/$i.pgm
    echo "input -> $i" | tee -a $OUTPUTGRAPH
    echo "threads = 1"
    RESULTS=()
    for r in {1..10} ## number of runs for sequential
    do
      comm="./$k $i $folder/$output 3"
      echo try= $r comm= $comm
      time=$($comm 2>>$errorlog)
      echo - $time -
      RESULTS+=($time)
      sleep 2
    done
    THREEBEST=($(calc_3_lower "${RESULTS[@]}"))
    KBESTSEQ=$(calc_kbest "${THREEBEST[@]}")
    echo k-best score = $KBESTSEQ
    echo "(1, 1)" >> $OUTPUTGRAPH
    for t in {2..8} ## number of total threads here
    do
      echo "threads = $t"
      RESULTS=()
      for r in {1..10} ## number of runs for each parallel configuration
      do
        #comm="./$k $i $folder/$output 2 $t"
        comm="mpirun -bynode -mca btl ^openib -np $t ./$k $i $folder/$output 999"
        echo try= $r comm= $comm
        time=$($comm 2>>$errorlog)
        echo - $time -
        RESULTS+=($time)
        sleep 2
      done
      THREEBEST=($(calc_3_lower "${RESULTS[@]}"))
      KBESTPAR=$(calc_kbest "${THREEBEST[@]}")
      echo k-best score = $KBESTPAR
      echo speedup = $(calc $KBESTSEQ/$KBESTPAR)
      echo "($t, $(calc $KBESTSEQ/$KBESTPAR))" >> $OUTPUTGRAPH
    done
  done
done

echo "That's all."

exit 0;

