./wham-2d Px=0 0.0 1.0 20 Py=0 120.0 175.0 25 0.00001 300 0 metadata pmf-mask-t5.dat 1 > wham-mask-t5.out
./wham-2d Px=0 0.0 1.0 20 Py=0 120.0 175.0 25 0.000001 300 0 metadata pmf-mask-t6.dat 1 > wham-mask-t6.out
python draw.py pmf-mask-t5.dat;mv pmf-p.png pmf-t5.png
python draw.py pmf-mask-t6.dat;mv pmf-p.png pmf-t6.png
