"""Python mirror of jm.c logic with invariant auditing, n=4."""
import random
n = 4
rng = random.Random(123)
T = [[[0]*n for _ in range(n)] for _ in range(n)]
rs = [[[] for _ in range(n)] for _ in range(n)]   # rs[r][s] = cols
cs = [[[] for _ in range(n)] for _ in range(n)]   # cs[c][s] = rows
rc = [[[] for _ in range(n)] for _ in range(n)]   # rc[r][c] = syms
neg = None
for r in range(n):
    for c in range(n):
        s = (r+c) % n
        T[r][c][s] = 1
        rs[r][s].append(c); cs[c][s].append(r); rc[r][c].append(s)

def inc(r,c,s):
    global neg
    T[r][c][s] += 1
    if T[r][c][s] == 1:
        rs[r][s].append(c); cs[c][s].append(r); rc[r][c].append(s)

def dec(r,c,s):
    global neg
    T[r][c][s] -= 1
    if T[r][c][s] == 0:
        rs[r][s].remove(c); cs[c][s].remove(r); rc[r][c].remove(s)
    else:
        neg = (r,c,s)

def move():
    global neg
    if neg is None:
        r = rng.randrange(n); c = rng.randrange(n)
        s0 = rc[r][c][0]
        while True:
            s = rng.randrange(n)
            if s != s0: break
        c0 = rs[r][s][0]
        r0 = cs[c][s][0]
    else:
        r, c, s = neg
        r0 = cs[c][s][rng.randrange(2)] if len(cs[c][s])>1 else cs[c][s][0]
        c0 = rs[r][s][rng.randrange(2)] if len(rs[r][s])>1 else rs[r][s][0]
        s0 = rc[r][c][rng.randrange(2)] if len(rc[r][c])>1 else rc[r][c][0]
    neg = None
    inc(r,c,s); dec(r,c0,s); dec(r0,c,s); dec(r,c,s0)
    inc(r,c0,s0); inc(r0,c,s0); inc(r0,c0,s)
    dec(r0,c0,s0)

def audit():
    # line sums and list consistency
    negcount = 0
    for r in range(n):
        for c in range(n):
            for s in range(n):
                assert T[r][c][s] in (-1,0,1)
                if T[r][c][s] == -1: negcount += 1
    assert negcount == (0 if neg is None else 1)
    for r in range(n):
        for s in range(n):
            assert sorted(rs[r][s]) == sorted(c for c in range(n) if T[r][c][s]==1)
            assert sum(T[r][c][s] for c in range(n)) == 1
    for c in range(n):
        for s in range(n):
            assert sorted(cs[c][s]) == sorted(r for r in range(n) if T[r][c][s]==1)
    for r in range(n):
        for c in range(n):
            assert sorted(rc[r][c]) == sorted(s for s in range(n) if T[r][c][s]==1)
    # improper: lines through neg have 2 entries
    if neg is not None:
        r,c,s = neg
        assert len(rs[r][s])==2 and len(cs[c][s])==2 and len(rc[r][c])==2

from collections import Counter
cnt = Counter()
for i in range(200000):
    move()
    if i % 997 == 0: audit()
    if neg is None and i > 5000:
        key = tuple(tuple(rc[r][c][0] for c in range(n)) for r in range(n))
        cnt[key] += 1
tot = sum(cnt.values())
print("distinct squares visited:", len(cnt), "(should be 576)")
import statistics
freqs = [v/tot for v in cnt.values()]
print(f"freq min {min(freqs):.6f} max {max(freqs):.6f} mean {1/len(cnt):.6f}")
# type distribution
tc = Counter()
for key, v in cnt.items():
    L = key
    sig=[0]*n
    for j in range(n): sig[L[0][j]]=L[1][j]
    seen=[False]*n; lens=[]
    for x0 in range(n):
        if not seen[x0]:
            l=0; x=x0
            while not seen[x]: seen[x]=True; x=sig[x]; l+=1
            lens.append(l)
    tc[tuple(sorted(lens))]+=v
for k,v in sorted(tc.items()): print(k, v/tot)
print("audits passed")
