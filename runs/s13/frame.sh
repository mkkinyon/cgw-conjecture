#!/bin/sh
for m in 1 2; do for k in 8 16 32 64 128 256; do ./s13_frame $m $k 1000000 0 4096 5; done; done
for N in 50 100 200; do for k in 8 16 24; do ./s13_frame 2 $k 1000000 0 $N 9; done; done
