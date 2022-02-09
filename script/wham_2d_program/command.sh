#for i in $(seq 1 50); do k=$(echo "1-0.02*${i}" | bc | awk '{printf "%.2f\n", $0}'); echo "${i}-wham.dat ${k} 0 200 0" >> metadata; done
#./wham-2d Px=0 0.0 1.0 20 Py=0 120.0 175.0 25 0.00001 300 0 metadata pmf-mask-t5.dat 1 > wham-mask-t5.out
#./wham-2d Px=0 0.0 1.0 20 Py=0 120.0 175.0 25 0.000001 300 0 metadata pmf-mask-t6.dat 1 > wham-mask-t6.out
./wham-2d Px=0 0.4 1.2 25 Py=0 30.0 40.0 25 0.00001 300 0 metadata pmf-mask-t5.dat 1 > wham-mask-t5.out
./wham-2d Px=0 0.4 1.2 25 Py=0 30.0 40.0 25 0.000001 300 0 metadata pmf-mask-t6.dat 1 > wham-mask-t6.out
python draw.py pmf-mask-t5.dat;mv pmf-p.png pmf-t5.png
python draw.py pmf-mask-t6.dat;mv pmf-p.png pmf-t6.png
