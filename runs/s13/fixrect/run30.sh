#!/bin/bash
# atypical lambda battery at n=30 on completions of a fixed rectangle
cd /home/claude/cgw
run() { spec=$1; cls=$2; tag=$3; init=$4; thin=$5
  python3 s13_fixrect.py init 30 $spec $init > runs/s13/fixrect/init30_$tag.bin
  ./jmfix 30 60000 $thin $((init+40)) < runs/s13/fixrect/init30_$tag.bin | python3 s13_fixstats.py 30 --marks 6 --classes $cls > runs/s13/fixrect/n30_$tag.log 2>&1
}
run cyc            30:2,30:5,30:15        cyc_s1    1 2000 &
run cyc            30:2,30:5,30:15        cyc_s2    2 500 &
run cyc:15,15      15:2,15:7              15-15     1 2000 &
run cyc:28,2       28:2,28:14             28-2      1 2000 &
run cyc:4,2,2,2,2,2,2,2,2,2,2,2,2,2  4:2  4-2s      1 2000 &
run cyc:6,3,3,3,3,3,3,3,3  6:2,6:3        6-3s      1 2000 &
run cyc:10,10,10   10:2,10:5              10s       1 2000 &
wait
echo done
