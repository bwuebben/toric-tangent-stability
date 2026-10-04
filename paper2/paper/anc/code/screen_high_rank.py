"""Exact non-stability screen for the high-Picard-rank slices (dims 7-9).
For each variety: exact anticanonical degrees (accelerated hull, exact certification, as in
../src/toric_fast.py), then a search for a ray-spanned subspace V with slope >= mu(T_X) in
increasing dimension, stopping at the first witness (slope > mu: unstable; slope = mu: not
stable). Varieties with no witness of dimension <= KMAX get the full analysis (analyze_fast),
which decides stability exactly. Output: output/screen_<n>d.json + log."""
import sys, os, json, time
from fractions import Fraction as Fr
from itertools import combinations
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
import toric_fast as tf
from toric_stability import dot, rref, in_span, barycenter_and_volume_from_facets
import math

KMAX = 4

def degrees(rays):
    n = len(rays[0])
    dual = tf.dual_polytope_vertices_fast(rays)
    deg, tris = [], []
    for u in rays:
        fv = [m for m in dual if dot(m, u) == -1]
        assert len(fv) >= n
        v, t = tf.facet_volume_and_triangulation_fast(u, fv)
        deg.append(v); tris.append(t)
    total = sum(deg)
    bary, vol = barycenter_and_volume_from_facets(n, tris)
    assert total == math.factorial(n) * vol
    return deg, total, bary

def screen(rays, deg, mu, kmax):
    r = len(rays); seen = set(); best = None
    for k in range(1, kmax + 1):
        for comb in combinations(range(r), k):
            rk, canon = rref([rays[i] for i in comb])
            if rk < k or canon in seen: continue
            seen.add(canon)
            members = [i for i in range(r) if in_span(canon, rays[i])]
            slope = Fr(sum(deg[i] for i in members), k)
            if best is None or slope > best[0]:
                best = (slope, k, [rays[i] for i in members])
            if slope >= mu:
                return {'dim': k, 'rays_in_V': [rays[i] for i in members], 'slope': [slope.numerator, slope.denominator],
                        'strict': slope > mu}, best
    return None, best

def work(item):
    vid, rays = item['id'], [tuple(v) for v in item['vertices']]
    n = len(rays[0]); t0 = time.time()
    try:
        deg, total, bary = degrees(rays)
        mu = Fr(total, n)
        w, best = screen(rays, deg, mu, KMAX)
        rec = {'id': vid, 'dim': n, 'num_rays': len(rays), 'picard_rank': len(rays) - n, 'degrees': deg,
               'anticanonical_degree': total, 'mu_TX': [mu.numerator, mu.denominator],
               'kahler_einstein': all(c == 0 for c in bary)}
        if w is not None:
            rec['status'] = 'unstable' if w['strict'] else 'not-stable (equality witness)'
            rec['witness'] = w
        else:
            full = tf.analyze_fast(rays, name=vid)
            rec['status'] = full['verdict']; rec['witness'] = full['witness']
            rec['equality_witnesses'] = full['equality_witnesses']; rec['full_analysis'] = True
        rec['seconds'] = round(time.time() - t0, 1)
        return rec
    except Exception as e:
        return {'id': vid, 'status': 'ERROR', 'error': repr(e)}

def main(n, path, procs, min_rho=None):
    items = json.load(open(path))['polytopes']
    if min_rho is not None:
        items = [it for it in items if len(it['vertices']) - n >= min_rho]
    out_path = os.path.join(HERE, '..', 'output', f'screen_{n}d.json')
    done = {}
    if os.path.exists(out_path):
        done = {r['id']: r for r in json.load(open(out_path))}
    todo = [it for it in items if it['id'] not in done]
    print(f"dim {n}: {len(items)} varieties, {len(done)} done, {len(todo)} to do, {procs} procs", flush=True)
    results = list(done.values()); t0 = time.time()
    with Pool(procs) as pool:
        for i, rec in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            results.append(rec)
            if (i + 1) % 25 == 0 or rec.get('full_analysis') or rec['status'] in ('stable', 'ERROR'):
                print(f"  [{i+1}/{len(todo)} {time.time()-t0:.0f}s] {rec['id']} rho={rec.get('picard_rank')} {rec['status']}"
                      + (f" wdim={rec['witness']['dim']}" if rec.get('witness') else ''), flush=True)
            if (i + 1) % 100 == 0:
                json.dump(results, open(out_path, 'w'))
    json.dump(results, open(out_path, 'w'))
    from collections import Counter
    print("summary by (rho, status):", sorted(Counter((r.get('picard_rank'), r['status']) for r in results).items()), flush=True)

if __name__ == '__main__':
    n = int(sys.argv[1]); procs = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    min_rho = int(sys.argv[3]) if len(sys.argv) > 3 else None
    path = {7: 'smooth_toric_fano_7d_rho_ge_9.json', 8: 'smooth_toric_fano_8d_rho_ge_10.json', 9: 'smooth_toric_fano_9d_rho_ge_11.json'}[n]
    main(n, os.path.join(HERE, '..', 'data', path), procs, min_rho)
