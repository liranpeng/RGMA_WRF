#!/usr/bin/env python
"""Build online/diag_fn_ml_offline.npz for each case.

This is get_fn_ml() from plot_allcases_lost_removed.py lifted verbatim (same
deployed TorchScript file, same 200k batching, same clip to [0,1], same npz
keys) so the cache is bit-identical to what that script would have written --
without running its plotting side and overwriting its own figures.
"""
import os
import sys

import numpy as np
import torch

TEST = ("/scratch/07088/tg863871/Perlm_Backup/"
        "WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test")
PT = os.path.join(os.path.dirname(TEST), "ml_models", "mlp_wrapper.pt")
CASES = ["20170714_1aer_org", "20170714_1aer_mid",
         "20170714_3aer_org", "20170714_3aer_mid"]

torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "8")))
print("emulator: %s" % PT, flush=True)
m = torch.jit.load(PT, map_location="cpu").eval()

for case in CASES:
    on = os.path.join(TEST, case, "online")
    cache = os.path.join(on, "diag_fn_ml_offline.npz")
    X = np.load(os.path.join(on, "parcel_inputs.npz"))["X"]

    if os.path.exists(cache):
        c = np.load(cache)
        if c["fn_ml"].shape[0] == X.shape[0]:
            print("%-22s cache already valid (%d rows)" % (case, X.shape[0]), flush=True)
            continue
        print("%-22s cache shape mismatch -> recomputing" % case, flush=True)

    print("%-22s running emulator on %s cells ..." % (case, f"{X.shape[0]:,}"), flush=True)
    outs = []
    with torch.no_grad():
        for s in range(0, X.shape[0], 200000):
            t = torch.from_numpy(np.ascontiguousarray(X[s:s + 200000], dtype=np.float32))
            outs.append(m(t).numpy())
    ml = np.concatenate(outs, 0)
    fn_ml = np.clip(ml[:, :2], 0, 1)          # slot0 = Aitken, slot1 = Accum
    np.savez(cache, fn_ml=fn_ml, ml_raw=ml)
    print("%-22s wrote %s  (mean Aitken=%.4f Accum=%.4f)"
          % (case, os.path.basename(cache), fn_ml[:, 0].mean(), fn_ml[:, 1].mean()),
          flush=True)

print("done", flush=True)
