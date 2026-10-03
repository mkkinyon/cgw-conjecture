#!/bin/bash
cd /home/claude/cgw
run() { n=$1; spec=$2; cls=$3; tag=$4; nsq=$5; thin=$6; marks=$7; seed=$8
  python3 s13_fixrect.py init $n $spec 1 > runs/s13/fixrect/winit${n}_$tag.bin
  ./jmfix $n $nsq $thin $seed < runs/s13/fixrect/winit${n}_$tag.bin | python3 s13_wtail.py $n --marks $marks --classes $cls > runs/s13/fixrect/w${n}_$tag.log 2>&1
}
( run 30 cyc 30:2,30:15 cyc 120000 1000 6 301
  run 50 cyc 50:2,50:25 cyc 60000 2500 6 501
  run 100 cyc 100:2,100:50 cyc 20000 10000 4 1001 ) &
( run 50 cyc:4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 4:2 4-2s 60000 2500 4 502
  run 50 cyc:25,25 25:2,25:12 25-25 60000 2500 6 503
  run 100 cyc:4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 4:2 4-2s 20000 10000 4 1002 ) &
wait; echo done
