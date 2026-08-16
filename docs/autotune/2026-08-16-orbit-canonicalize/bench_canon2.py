import sys, time, numpy as np
from ppvm._core import TranslationGroup
tag = sys.argv[1]
rng = np.random.default_rng(0)
out = {}
for L in (4, 8, 16, 32, 64, 96):
    g = TranslationGroup.ladder(L, 2); N = 2 * L
    words = []
    for _ in range(200):
        w = np.zeros(N, dtype=np.uint8)
        pos = rng.choice(N, size=4, replace=False)
        w[pos] = rng.integers(1, 4, size=4)
        words.append(w)
    g.canonicalize(words[0])
    reps = max(2, int(4000 // L))
    best = 1e9
    for _ in range(3):
        t0 = time.perf_counter()
        for _ in range(reps):
            for w in words:
                g.canonicalize(w)
        best = min(best, (time.perf_counter() - t0) / (reps * len(words)))
    out[L] = best * 1e6
    print(f"{tag} L={L:<3} |G|={L:<3} N={N:<4} {best*1e6:9.2f} us/call", flush=True)
np.save(f"/tmp/canon_{tag}.npy", np.array([[L, v] for L, v in out.items()]))
