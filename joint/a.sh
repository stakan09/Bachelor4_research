make
para=1
rm -f log
seq 10000 | xargs -I{} -P 8 sh -c './bregsplit {} 1 | grep -v "p\[" >> log'
sort -n -k 2 log | grep -v 'nan' |grep -v 'inf'| head -n 10 > random_result
cat random_result | awk "NR==${para}{print \$1}" |xargs -I{} ./bregsplit {} 0
gnuplot -e "
set term png;
set xrange[-1:1];
set out 'plot.png';
plot 'test.txt' w l lw 4 lc 3 title 'Calculation', x*x w l  lw 4 lc 7;
"

