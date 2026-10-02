#!/bin/sh
./jm 50 600 150000 201 | python3 s13_real.py 50 --inst 4 --pairs 16 --nlevel 16 --seed 12 --csv runs/s13/real50.csv > runs/s13/real50.log 2>&1
./jm 30 300 30000 202 | python3 s13_trade.py 30 --inst 3 --steps 3000 --seed 14 > runs/s13/trade30.log 2>&1
./jm 100 300 1000000 203 | python3 s13_real.py 100 --inst 4 --pairs 16 --nlevel 16 --seed 13 --csv runs/s13/real100.csv > runs/s13/real100.log 2>&1
./jm 50 200 150000 204 | python3 s13_trade.py 50 --inst 3 --steps 4000 --seed 15 > runs/s13/trade50.log 2>&1
./jm 100 100 1000000 205 | python3 s13_trade.py 100 --inst 3 --steps 6000 --seed 16 > runs/s13/trade100.log 2>&1
