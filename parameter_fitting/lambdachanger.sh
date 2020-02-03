#!/bin/bash

#filename="$(date '+%y%m%d%k%M')"
filename=L2_9ten_1e-5_gosanasi
g++-8.2 -O2 -lm main.cpp bregsplitlogistic.cpp -o logistic
mkdir $filename
cd $filename
#seq -6 1 -3 | awk '{print exp($1*log(10))}' | xargs -L 1 -P 4 -I{} .././logistic {}
seq 1 1 100 | awk '{print $1}' | xargs -L 1 -P 8 -I{} .././logistic {}
#for i in `seq -4 1 0` 
#do
#data no kazu to parameter no kazu de kaeru koto
#paste calc1_$i.000000.dat calc_$i.000000.dat | awk '{i+=($4-$2)*($4-$2)};NR==10{print sqrt(i)}' >> L1c.dat
#cat p1_$i.000000.dat | awk '{i+=$3*$3};NR==5{print sqrt(i)}' >> L1p.dat
#paste L1c.dat L1p.dat | awk '{print $1" "$2}' > L1.dat
#done
#gnuplot -e "
#set term png;
#set xlabel 'Cost' ;
#set ylabel 'Parameter norm' ;
#set out 'L1.png' ;
#plot 'L1.dat' w l lw 4 ;
#"
for i in `seq 1 1 100` ; do paste calc1_${i}.dat calcc.dat | awk '{i+=($2-$4)*($2-$4);if(NR % 9 == 0) print sqrt(i);if(NR % 9 == 0) i=0}' ; done > Lc.dat
for i in `seq 1 1 100` ; do paste p1a_${i}.dat pbase.dat | awk '{i+=($3-$6)*($3-$6);if(NR % 6 == 0) print sqrt(i);if(NR % 6 == 0) i=0}' ; done > Lp.dat

