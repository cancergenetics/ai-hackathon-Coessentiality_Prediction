"""H1/H1S design check from labels and covariates only (synthetic outcome; no dependency value): which events have
a leverage-1 / rank-deficient HC3 design. Writes results/part2/directions/H1S_silencing_GI/degenerate_designs.tsv.
Run from the project root."""
import importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path.cwd()
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("run_direction", ROOT / "scripts" / "44_run_direction.py")
R = importlib.util.module_from_spec(spec); sys.modules["run_direction"] = R; spec.loader.exec_module(R)
from part2 import assoc, covariates as covmod, stats, registry
env = R.real_env(dry=True)
for aid in ("H1_within_complex_GI", "H1S_silencing_GI"):
    entry, _ = R.entry_of(env, aid)
    events = R.load_events(aid, env.proj.root)
    loss = R.load_loss_A(env.proj.root, env.E_A.index)
    bad = []
    for e in events:
        ctx = loss[e.gene].rename(e.context)
        Z = covmod.assemble(covmod.covariates_for(entry, e.context)[0], env.cov)
        X, names, ok, dropped = assoc._design(ctx.reindex(env.E_A.index), Z.reindex(env.E_A.index), return_dropped=True)
        try:
            stats.HC3Design(X, 1); msg = ""
        except (np.linalg.LinAlgError, ValueError) as ex:
            msg = str(ex)
            h = None
        if msg:
            # which column causes leverage 1
            Q, Rr = np.linalg.qr(X); h = np.einsum("ij,ij->i", Q, Q)
            i1 = np.flatnonzero(h > 1 - 1e-10)
            cols = [names[j] for j in range(X.shape[1]) if X[i1, j].any() and X[:, j].sum() <= 3] if len(i1) else []
            bad.append((e.id, e.n_A, int(X[:, 1].sum()), X.shape, msg, list(env.E_A.index[ok][i1]), cols))
    print(aid, len(events), "events; design failures:", len(bad))
    for b in bad: print("  ", b)
    if aid.startswith("H1S"):
        import pandas as pd
        pd.DataFrame(bad, columns=["event","n_A","n_loss_in_design","design_shape","error","leverage1_models","sparse_columns_at_those_rows"]).to_csv("results/part2/directions/H1S_silencing_GI/degenerate_designs.tsv", sep="\t", index=False)
