"""s13_fixrect.py -- completions of a fixed 2 x n rectangle: initial completion + driver for ./jmfix.

  python3 s13_fixrect.py init n sigma_spec > init.bin        write one completion (rows 0,1 = (id, sigma))
  python3 s13_fixrect.py enum n sigma_spec                   enumerate all completions (small n), print count,
                                                              and write them to runs/s13/fixrect/all_n{n}_{spec}.bin
  python3 s13_fixrect.py validate n sigma_spec nsamp thin    chi-square of ./jmfix samples against the enumeration,
                                                              from two different initial completions
sigma_spec: 'cyc' (single n-cycle), 'cyc:a,b,c' (cycle type), or an explicit permutation 'perm:3,0,1,2,...'
Rows are 0-based; symbols 0..n-1; row 0 is the identity 0..n-1 and row 1 is sigma (sigma(L[0][j]) = L[1][j]).
"""
import sys, os, random, subprocess, math
from collections import Counter


def sigma_from_spec(n, spec):
    if spec == 'cyc':
        return [(i + 1) % n for i in range(n)]
    if spec.startswith('cyc:'):
        parts = [int(t) for t in spec[4:].split(',')]
        assert sum(parts) == n and all(p >= 2 for p in parts), spec
        sig = [0] * n; start = 0
        for p in parts:
            for i in range(p): sig[start + i] = start + (i + 1) % p
            start += p
        return sig
    if spec.startswith('perm:'):
        sig = [int(t) for t in spec[5:].split(',')]
        assert sorted(sig) == list(range(n)) and all(sig[i] != i for i in range(n)), spec
        return sig
    raise ValueError(spec)


def complete(rows, n, rng=None):
    """complete a Latin rectangle (list of rows) to a Latin square, row by row by bipartite matching
    (columns -> unused symbols); random order of trial if rng given.  Hall's theorem guarantees success."""
    rows = [r[:] for r in rows]
    used_col = [set(r[c] for r in rows) for c in range(n)]
    while len(rows) < n:
        avail = [[s for s in range(n) if s not in used_col[c]] for c in range(n)]
        if rng is not None:
            for a in avail: rng.shuffle(a)
        match_sym = [-1] * n   # symbol -> column
        def try_col(c, seen):
            for s in avail[c]:
                if s in seen: continue
                seen.add(s)
                if match_sym[s] < 0 or try_col(match_sym[s], seen):
                    match_sym[s] = c; return True
            return False
        order = list(range(n))
        if rng is not None: rng.shuffle(order)
        for c in order:
            ok = try_col(c, set()); assert ok
        row = [0] * n
        for s in range(n): row[match_sym[s]] = s
        rows.append(row)
        for c in range(n): used_col[c].add(row[c])
    return rows


def to_bytes(L):
    return bytes(x for row in L for x in row)


def enumerate_completions(rows, n, out):
    """all completions by backtracking cell by cell; writes squares to `out`, returns count."""
    L = [r[:] for r in rows] + [[-1] * n for _ in range(n - len(rows))]
    colused = [set(r[c] for r in rows) for c in range(n)]
    rowused = [set(r) for r in rows] + [set() for _ in range(n - len(rows))]
    count = 0
    r0 = len(rows)
    cells = [(r, c) for r in range(r0, n) for c in range(n)]
    def rec(i):
        nonlocal count
        if i == len(cells):
            count += 1; out.write(to_bytes(L)); return
        r, c = cells[i]
        for s in range(n):
            if s in rowused[r] or s in colused[c]: continue
            L[r][c] = s; rowused[r].add(s); colused[c].add(s)
            rec(i + 1)
            rowused[r].discard(s); colused[c].discard(s); L[r][c] = -1
    rec(0)
    return count


def run_jmfix(n, init, nsamp, thin, seed):
    p = subprocess.run(['./jmfix', str(n), str(nsamp), str(thin), str(seed)], input=to_bytes(init), capture_output=True)
    sys.stderr.write(p.stderr.decode())
    return p.stdout


def main():
    cmd = sys.argv[1]; n = int(sys.argv[2]); spec = sys.argv[3]
    sig = sigma_from_spec(n, spec)
    rows = [list(range(n)), sig[:]]
    if cmd == 'init':
        seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        L = complete(rows, n, random.Random(seed))
        sys.stdout.buffer.write(to_bytes(L)); return
    os.makedirs('runs/s13/fixrect', exist_ok=True)
    tag = spec.replace(':', '').replace(',', '-')
    allf = f'runs/s13/fixrect/all_n{n}_{tag}.bin'
    if cmd == 'enum':
        with open(allf, 'wb') as f: cnt = enumerate_completions(rows, n, f)
        print(f'n={n} sigma={spec}: {cnt} completions -> {allf}'); return
    if cmd == 'validate':
        nsamp = int(sys.argv[4]); thin = int(sys.argv[5])
        if not os.path.exists(allf):
            with open(allf, 'wb') as f: cnt = enumerate_completions(rows, n, f)
        data = open(allf, 'rb').read(); N = n * n; cnt = len(data) // N
        index = {data[i * N:(i + 1) * N]: i for i in range(cnt)}
        print(f'n={n} sigma={spec}: {cnt} completions; sampling {nsamp} x thin {thin} from two starts')
        for start_seed in (1, 2):
            init = complete(rows, n, random.Random(start_seed))
            buf = run_jmfix(n, init, nsamp, thin, 100 + start_seed)
            got = Counter()
            bad = 0
            for i in range(len(buf) // N):
                sq = buf[i * N:(i + 1) * N]
                if sq in index: got[index[sq]] += 1
                else: bad += 1
            hit = len(got); exp = nsamp / cnt
            chi2 = sum((got.get(i, 0) - exp) ** 2 / exp for i in range(cnt))
            df = cnt - 1
            z = (chi2 - df) / math.sqrt(2 * df)
            print(f'  start {start_seed}: not-a-completion={bad}  distinct states hit={hit}/{cnt}  chi2={chi2:.1f} df={df}  z=(chi2-df)/sqrt(2df)={z:+.2f}')
        return
    raise SystemExit('unknown command')


if __name__ == '__main__':
    main()
