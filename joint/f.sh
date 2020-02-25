#!/bin/sh
jikken_data=3spalt.txt
calcdata=plotdataS3sp
para=plot
paste ${jikken_data1} ${calcdata} | awk ' {print $1" "($6-$2)/$2" "($7-$3)/$3" "($8-$4)/$4} ' > ${para}_3spsoutaigosanew.txt
paste ${jikken_data1} ${calcdata} | awk ' {print $1" "($6-$2)" "($7-$3)" "($8-$4)} ' > ${para}_3spzettaigosanew.txt
gnuplot -e "
set term png;
set ylabel 'Processing Speed';
set xlabel 'Condition';
set mxtics 2;
set key outside;
set format x '';
set xtics('1' 0.5,'2' 1.5,'3' 2.5,'4' 3.5,'5' 4.5,'6' 5.5,'7' 6.5,'8' 7.5,'9' 8.5);
set title 'SCN';
set out '${para}_3sp_data0.png';
plot '${calcdata}' w l lw 4 title 'Calculation', '${jikken_data}' w l lw 4 lc 3 title 'Experiment';
set title 'Phenol';
set out '${para}_3sp_data1.png';
plot '${calcdata}' u 1:3  w l lw 4 title 'Calculation', '${jikken_data}' u 1:3 w l lw 4 lc 3 title 'Experiment';
set title 'S2O3';
set out '${para}_3sp_data2.png';
plot '${calcdata}' u 1:4  w l lw 4 title 'Calculation', '${jikken_data}' u 1:4 w l lw 4 lc 3 title 'Experiment';
set ylabel 'difference';
set xlabel 'condition';
set title 'SCN';
set out '${para}_3spSCNgosa.png';
plot '${para}_3spsoutaigosanew.txt' u 1:2 w l lw 4 title 'difference';
set title 'Phenol';
set out '${para}_3spPhegosa.png';
plot '${para}_3spsoutaigosanew.txt' u 1:3 w l lw 4 title 'difference';
set title 'S2O3';
set out '${para}_3spS2O3gosa.png';
plot '${para}_3spsoutaigosanew.txt' u 1:4 w l lw 4 title 'difference';
"

jikken_data=plot_SCN
calcdata=plotdataSSCN
paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)/$2" "($7-$3)/$3" "($8-$4)/$4} ' > ${para}_SCNsoutaigosa.txt
paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)" "($7-$3)" "($8-$4)} ' > ${para}_SCNzettaigosa.txt
gnuplot -e "
set term png;
set ylabel 'Processing Speed';
set xlabel 'Condition';
set mxtics 2;
set key outside;
set format x '';
set xtics('1' 0.5,'2' 1.5,'3' 2.5,'4' 3.5,'5' 4.5,'6' 5.5,'7' 6.5,'8' 7.5,'9' 8.5);
set xrange [0:6.99];
set title 'SCN';
set out '${para}_SCN_data.png';
plot '${calcdata}' w l lw 4 title 'Calculation', '${jikken_data}' w l lw 4 lc 3 title 'Experiment';
set ylabel 'difference';
set xlabel 'condition';
set title 'SCN';
set out '${para}_SCNgosa.png';
plot '${para}_SCNsoutaigosa.txt' u 1:2 w l lw 4 title 'difference';
"

