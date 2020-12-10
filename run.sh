#!/bin/bash

helpFunction() {
   echo "Usage: $0 -n -p"
   echo -e "\t-n <number> number of times that the programs will run"
   echo -e "\t-t <n° of threads> sets the number of processes"
}

if [ $# -eq 0 ] 
    then
    helpFunction
    exit 1
fi

while getopts "n:t:h" opt
do
   case "$opt" in
      n ) n_repeticoes="$OPTARG" ;;
      t ) n_threads="$OPTARG" ;;
      h ) helpFunction; exit 1 ;; # Print helpFunction in case parameter is non-existent
      * ) exit 1
   esac
done

par_executable="lu_par.out"

sizes=(1000 2000 4000 8000)

seq_time=()
par_time=()

echo "Execução para ${n_threads} threads:"
echo
echo

for size in ${sizes[@]}; do
    echo "Size ${size}"
    for i in $(eval echo "{1..$n_repeticoes}"); do
        par_time+=($(upcrun -shared-heap 2GB -n ${n_threads} ${par_executable} -n ${size} | grep "Time elapsed" | grep -Eo "?([0-9]*[.])?[0-9]+"))
        seq_time+=($(upcrun -shared-heap 2GB -n 1 ${par_executable} -n ${size} | grep "Time elapsed" | grep -Eo "?([0-9]*[.])?[0-9]+"))
    done

    python3 speedup.py $size $n_threads ${n_repeticoes} ${seq_time[@]} ${par_time[@]}

    seq_time=()
    par_time=()
done 