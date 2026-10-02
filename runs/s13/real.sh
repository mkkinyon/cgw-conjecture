#!/bin/sh
./jm 50 600 150000 101 | python3 s13_real.py 50 --inst 4 --pairs 16 --nlevel 16 --seed 2 --csv runs/s13/real50.csv > runs/s13/real50.log 2>&1
./jm 100 300 1000000 102 | python3 s13_real.py 100 --inst 4 --pairs 16 --nlevel 16 --seed 3 --csv runs/s13/real100.csv > runs/s13/real100.log 2>&1
./jm 200 150 8000000 103 | python3 s13_real.py 200 --inst 4 --pairs 16 --nlevel 16 --seed 4 --csv runs/s13/real200.csv > runs/s13/real200.log 2>&1
