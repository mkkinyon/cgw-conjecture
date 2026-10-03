#!/bin/bash
cd /home/claude/cgw
python3 s13_fixrect.py init 150 cyc 1 > runs/s13/fixrect/init150_cyc.bin
python3 s13_fixrect.py init 150 cyc:4,$(python3 -c "print(','.join(['2']*73))") 1 > runs/s13/fixrect/init150_4-2s.bin
( ./jmfix 150 6000 22500 151 < runs/s13/fixrect/init150_cyc.bin | python3 s13_firstpair_id.py 150 --marks 3 --classes 150:2,150:75 > runs/s13/fixrect/fp150_cyc.log 2>&1
  ./jmfix 50 40000 2500 52 < runs/s13/fixrect/init50_25-25.bin | python3 s13_firstpair_id.py 50 --marks 4 --classes 25:2,25:12 > runs/s13/fixrect/fp50_25-25.log 2>&1 ) &
( ./jmfix 150 6000 22500 152 < runs/s13/fixrect/init150_4-2s.bin | python3 s13_firstpair_id.py 150 --marks 3 --classes 4:2 > runs/s13/fixrect/fp150_4-2s.log 2>&1
  ./jmfix 50 40000 2500 53 < runs/s13/fixrect/init50_4-2s.bin | python3 s13_firstpair_id.py 50 --marks 4 --classes 4:2 > runs/s13/fixrect/fp50_4-2s.log 2>&1
  ./jmfix 100 12000 10000 54 < runs/s13/fixrect/init100_4-2s.bin | python3 s13_firstpair_id.py 100 --marks 3 --classes 4:2 > runs/s13/fixrect/fp100_4-2s.log 2>&1 ) &
wait; echo done
