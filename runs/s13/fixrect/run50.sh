#!/bin/bash
cd /home/claude/cgw
run() { n=$1; spec=$2; cls=$3; tag=$4; nsq=$5; thin=$6; marks=$7
  python3 s13_fixrect.py init $n $spec 1 > runs/s13/fixrect/init${n}_$tag.bin
  ./jmfix $n $nsq $thin 77 < runs/s13/fixrect/init${n}_$tag.bin | python3 s13_fixstats.py $n --marks $marks --classes $cls --cand 8 > runs/s13/fixrect/n${n}_$tag.log 2>&1
}
( run 50 cyc 50:2,50:25 cyc 40000 2500 4
  run 50 cyc:4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 4:2 4-2s 40000 2500 4
  run 50 cyc:48,2 48:2,48:24 48-2 40000 2500 4 ) &
( run 50 cyc:25,25 25:2,25:12 25-25 40000 2500 4
  run 50 cyc:6,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2 6:2,6:3 6-3s 40000 2500 4
  run 100 cyc 100:2,100:50 cyc 12000 10000 3
  run 100 cyc:4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 4:2 4-2s 12000 10000 3 ) &
wait; echo done
