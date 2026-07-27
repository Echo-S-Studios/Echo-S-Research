#!/usr/bin/env python3
"""run_all.py — the reproduction harness for the Exact Arithmetic evidence bundle.

    python3 run_all.py              # quick: full suite, census in smoke mode   (~90 s)
    python3 run_all.py --full       # everything, census over the whole box     (~3.5 min)
    python3 run_all.py --list       # environment + expected runtimes, no execution
    python3 run_all.py --verify     # checksums only

Exit status is 0 only if every checksum matches, every script exits 0, and every script's
SUMMARY line reports no failures. A source bundle is reproducible only when the environment and
the execution expectations travel with it, so both are printed before anything runs.
"""
import argparse, hashlib, re, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ── the environment this bundle was produced and last verified under ───────────────────────────
PINNED = {'python': '3.12.3', 'sympy': '1.14.0', 'mpmath': '1.3.0'}
# Known-compatible: CPython 3.12-3.13, SymPy 1.13-1.14, mpmath 1.3.
# SymPy < 1.13 changed the shape of Poly.intervals(all=True); lesson6_census.py tolerates both,
# but nothing older than 1.12 has been exercised.

#   name                        quick(s) full(s)  what it decides
SUITE = [
    ('lesson1_checks.py',            1,    1, 'Mahler measure, Lehmer by trace fold, the two walls'),
    ('lesson2_checks.py',            2,    2, 'semiring + Adams laws, tropical counterexamples, trace duality'),
    ('lesson3_checks.py',            2,    2, 'KL expansion, trace form, lambda = 2c, gate ladder, the flip'),
    ('lesson4_checks.py',            2,    2, 'capacity gate, charge groups, parity floor, the commutator door'),
    ('lesson5_checks.py',            1,    1, 'projector/capture, compositum, Cl(2,0), the phi-slack lexicon'),
    ('lesson6_checks.py',            3,    3, 'gauge/rigidity, the contact ledger, modulus pinning'),
    ('lesson6_census.py',           12,  130, 'the [-2,2]^5 quintic census; --full re-derives the D5 erratum'),
    ('lesson7_checks.py',            2,    2, 'seed ledger, gate ladder, the chart, q = i*tau'),
    ('lesson8_checks.py',           12,   12, 'closure guard, 27-subfield census, similarity, Emptiness IV-V'),
    ('lesson9_checks.py',            3,    3, 'torsion witnesses, ATR units, the orbit criterion'),
    ('lesson10_checks.py',           3,    3, 'Minkowski/index, the Phi retraction, the small-Salem box'),
    ('revision_checks.py',           1,    1, 'third-pass repairs, verified before they were written'),
    ('synthesis_checks.py',          1,    1, 'the Schinzel scoping and the grafted lessons'),
    ('certificates.py',             45,   45, 'Lessons 11-13: exact signs, exact rank, interval covariance'),
]
QUICK_ARGS = {'lesson6_census.py': ['--quick']}


def environment():
    rows = [('python', '.'.join(map(str, sys.version_info[:3])))]
    for mod in ('sympy', 'mpmath'):
        try:
            rows.append((mod, __import__(mod).__version__))
        except Exception:
            rows.append((mod, 'MISSING'))
    return rows


def show_env():
    print('environment')
    drift = False
    for name, got in environment():
        want = PINNED[name]
        flag = '' if got == want else f'   <- pinned {want}'
        if got != want:
            drift = True
        print(f'    {name:8s} {got}{flag}')
    if drift:
        print('    (version drift is not automatically a failure; the suite is known-compatible\n'
              '     across CPython 3.12-3.13, SymPy 1.13-1.14, mpmath 1.3)')
    return drift


def verify_sums():
    sums = HERE / 'SHA256SUMS'
    if not sums.exists():
        print('SHA256SUMS: not found — cannot verify provenance')
        return False
    ok = True
    for line in sums.read_text().splitlines():
        if not line.strip():
            continue
        want, name = line.split(None, 1)
        f = HERE / name.strip()
        if not f.exists():
            print(f'    MISSING  {name.strip()}')
            ok = False
            continue
        got = hashlib.sha256(f.read_bytes()).hexdigest()
        if got != want:
            print(f'    MISMATCH {name.strip()}\n             want {want}\n             got  {got}')
            ok = False
    print('checksums: ' + ('all match' if ok else 'FAILED'))
    return ok


def run(full):
    print(f'\n{"script":<24}{"time":>7}   result')
    print('-' * 78)
    failures, total = [], 0.0
    for name, tq, tf, _ in SUITE:
        args = [] if full else QUICK_ARGS.get(name, [])
        t0 = time.time()
        try:
            r = subprocess.run([sys.executable, str(HERE / name)] + args,
                               capture_output=True, text=True, timeout=3600)
            out = r.stdout
        except subprocess.TimeoutExpired:
            out, r = '', None
        dt = time.time() - t0
        total += dt
        line = ''
        for L in reversed(out.splitlines()):
            if 'SUMMARY' in L:
                line = L.strip()
                break
        bad = (r is None) or r.returncode != 0 or not line or 'FAILURES' in line
        if bad:
            failures.append(name)
            line = line or ('TIMEOUT' if r is None else f'no SUMMARY (exit {r.returncode})')
        tag = '' if not args else ' ' + ' '.join(args)
        print(f'{name + tag:<24}{dt:6.1f}s   {line}')
    print('-' * 78)
    print(f'{"total":<24}{total:6.1f}s   '
          + ('ALL GREEN' if not failures else f'FAILURES: {failures}'))
    return not failures


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--full', action='store_true', help='full census instead of the smoke slice')
    ap.add_argument('--list', action='store_true', help='print the plan and exit')
    ap.add_argument('--verify', action='store_true', help='checksums only')
    a = ap.parse_args()

    print('Exact Arithmetic — evidence bundle\n')
    show_env()
    print()

    if a.list:
        tq = sum(s[1] for s in SUITE)
        tf = sum(s[2] for s in SUITE)
        print(f'{"script":<24}{"quick":>7}{"full":>7}   decides')
        print('-' * 100)
        for name, q, f, what in SUITE:
            print(f'{name:<24}{q:6d}s{f:6d}s   {what}')
        print('-' * 100)
        print(f'{"expected total":<24}{tq:6d}s{tf:6d}s')
        print('\nCounts: 16 32 33 25 23 20 6 39 20 32 26 13 21 22.')
        print('354 numbered checks across the thirteen lessons; the census (6), the two revision')
        print('passes (13, 21) and the certificates (22 executable + 1 cited) are separate artifacts.')
        return 0

    ok_sums = verify_sums()
    if a.verify:
        return 0 if ok_sums else 1
    ok_run = run(a.full)
    if not a.full:
        print('\nnote: --quick ran lesson6_census.py on a slice of the box. The published census '
              'tallies\n      (625/50/638/1318/411 -> 83, and the D5 erratum) require --full.')
    return 0 if (ok_sums and ok_run) else 1


if __name__ == '__main__':
    sys.exit(main())
