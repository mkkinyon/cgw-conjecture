#!/bin/sh
./jm 30 200 30000 78 | python3 s13_trade.py 30 --inst 3 --steps 3000 --seed 4 > runs/s13/trade30.log 2>&1
./jm 50 200 150000 79 | python3 s13_trade.py 50 --inst 3 --steps 4000 --seed 5 > runs/s13/trade50.log 2>&1
./jm 100 100 1000000 80 | python3 s13_trade.py 100 --inst 3 --steps 6000 --seed 6 > runs/s13/trade100.log 2>&1
