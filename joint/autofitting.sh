#!/usr/bin/bash
set +o posix
make
echo //////////////////Searching Starting Point/////////////////////////
#seq 1 10 | xargs -L 1 -P 8 bash seed.sh 
seq 1 1000 | xargs -L 1 -P 8 ./plot 2>  >(tee log$1) >/dev/null
#bash -c "ls |grep -P '^log\d\$'|xargs cat" >log 
echo //////////////////Sorting data/////////////////////////////////////
sort -n -k 2 log | grep -v 'nan' |grep -v 'inf'| head -n 10 > random_result
echo //////////////////Searching Near by Linear Search//////////////////
for para in `seq 1 5`
do
     jikken_data=3spalt.txt
     calcdata=optdataS3sp
   cat random_result | awk "NR==${para}{print \$1}" |xargs ./plot | awk 'NR<=30{print $2}'| ./psearch 1>${para}_parameter.txt 
   #cat parameter | awk 'NR<=30{print $2}'| ./psearch 1>${para}_parameter.txt 
     paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)/$2" "($7-$3)/$3" "($8-$4)/$4} ' > ${para}_3spsoutaigosa.txt
     paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)" "($7-$3)" "($8-$4)} ' > ${para}_3spzettaigosa.txt
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
     plot '${calcdata}' w l lw 4 title 'Calculation', '${jikken_data}' w l lw 4 lc 3 title 'Experimenti','plot_3sp' w l lw 4 lc 7 title 'Experimentiold';
     set title 'Phenol';
     set out '${para}_3sp_data1.png';
     plot '${calcdata}' u 1:3  w l lw 4 title 'Calculation', '${jikken_data}' u 1:3 w l lw 4 lc 3 title 'Experiment','plot_3sp' u 1:3 w l lw 4 lc 7 title 'Experimentiold';
     set title 'S2O3';
     set out '${para}_3sp_data2.png';
     plot '${calcdata}' u 1:4  w l lw 4 title 'Calculation', '${jikken_data}' u 1:4 w l lw 4 lc 3 title 'Experiment','plot_3sp' u 1:4 w l lw 4 lc 7 title 'Experimentiold';
     set ylabel 'difference';
     set xlabel 'condition';
     set title 'SCN';
     set out '${para}_SCNgosa.png';
     plot '${para}_3spsoutaigosa.txt' u 1:2 w l lw 4 title 'difference';
     set title 'Phenol';
     set out '${para}_Phegosa.png';
     plot '${para}_3spsoutaigosa.txt' u 1:3 w l lw 4 title 'difference';
     set title 'S2O3';
     set out '${para}_S2O3gosa.png';
     plot '${para}_3spsoutaigosa.txt' u 1:4 w l lw 4 title 'difference';
     "

     jikken_data=plot_SCN
     calcdata=optdataSSCN
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

     jikken_data=plot_Phe
     calcdata=optdataSPhe
     paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)/$2" "($7-$3)/$3" "($8-$4)/$4} ' > ${para}_Phesoutaigosa.txt
     paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)" "($7-$3)" "($8-$4)} ' > ${para}_Phezettaigosa.txt
     gnuplot -e "
     set term png;
     set ylabel 'Processing Speed';
     set xlabel 'Condition';
     set mxtics 2;
     set key outside;
     set format x '';
     set xtics('1' 0.5,'2' 1.5,'3' 2.5,'4' 3.5,'5' 4.5,'6' 5.5,'7' 6.5,'8' 7.5,'9' 8.5);
     set xrange [0:9];
     set title 'Phenol';
     set out '${para}_Phe_data.png';
     plot '${calcdata}' u 1:3  w l lw 4 title 'Calculation', '${jikken_data}' u 1:3 w l lw 4 lc 3 title 'Experiment';
     set title 'Phenol';
     set out '${para}_Phegosa.png';
     plot '${para}_Phesoutaigosa.txt' u 1:3 w l lw 4 title 'difference';
     "
     jikken_data=plot_S2O3
     calcdata=optdataSS2O3
     paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)/$2" "($7-$3)/$3" "($8-$4)/$4} ' > ${para}_S2O3soutaigosa.txt
     paste ${jikken_data} ${calcdata} | awk ' {print $1" "($6-$2)" "($7-$3)" "($8-$4)} ' > ${para}_S2O3zettaigosa.txt
     gnuplot -e "
     set term png;
     set ylabel 'Processing Speed';
     set xlabel 'Condition';
     set mxtics 2;
     set key outside;
     set format x '';
     set xtics('1' 0.5,'2' 1.5,'3' 2.5,'4' 3.5,'5' 4.5,'6' 5.5,'7' 6.5,'8' 7.5,'9' 8.5,'10' 9.5);
     set xrange [0:10];
     set title 'S2O3';
     set out '${para}_S2O3_data.png';
     plot '${calcdata}' u 1:4  w l lw 4 title 'Calculation', '${jikken_data}' u 1:4 w l lw 4 lc 3 title 'Experiment';
     set title 'S2O3';
     set out '${para}_S2O3gosa.png';
     plot '${para}_S2O3soutaigosa.txt' u 1:4 w l lw 4 title 'difference';
     "

     mv optdataX3sp ${para}_X3sp.txt
     mv optdataXSCN ${para}_XSCN.txt
     mv optdataXPhe ${para}_XPhe.txt
     mv optdataXS2O3 ${para}_XS2O3.txt
     mkdir ${para}opt
     mv ${para}_* ${para}opt

done
echo Calculation done | mail -s "Calculation" rnc9ndr.yg@gmail.com  
