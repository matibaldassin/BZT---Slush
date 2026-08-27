#!/usr/bin/env python3
import argparse, numpy as np

def load_P(path):
    with open(path,'r') as f:
        nx, ny, nz = map(int, f.readline().split())
        _ = int(f.readline())  # Ncells, unused here
        data = []
        for line in f:
            s=line.strip()
            if not s: continue
            cols=s.split()
            if len(cols) < 4: continue
            step = int(float(cols[0]))
            px,py,pz = map(float, cols[1:4])  # units don't matter for rho(t)
            data.append((step, px, py, pz))
    arr = np.array(data, float)
    steps = arr[:,0].astype(int)
    P = arr[:,1:4]
    return steps, P  # P shape (N,3)

def autocorr_fft(x):
    """Return normalized autocorrelation rho(tau) for a 1D series x (length N)."""
    x = np.asarray(x, float)
    x -= x.mean()
    N = len(x)
    # Next power of 2 for zero-padding (speed + no circular wrap)
    nfft = 1 << (2*N-1).bit_length()
    f = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(f * np.conjugate(f), n=nfft)[:N]
    ac /= ac[0] if ac[0] != 0 else 1.0
    return ac  # rho[0]=1

def tau_int_sokal(rho, c=5.0):
    """
    Integrated autocorrelation time using Sokal's self-consistent window:
      tau = 0.5 + sum_{t=1}^{W} rho(t)
      W chosen as the largest t s.t. t <= c * tau
    """
    tau = 0.5
    W = 1
    N = len(rho)
    # Stop if rho becomes negative for stability
    while True:
        # extend W while within data and condition t <= c*tau holds
        while W < N and W <= c * tau and rho[W] > 0:
            tau += rho[W]
            W += 1
        # Check fixed point
        if W <= c * tau or W >= N or (W < N and rho[W] <= 0):
            break
        # Otherwise iterate again (rare)
    return tau, W

def main():
    ap = argparse.ArgumentParser(description="Estimate polarization autocorrelation time")
    ap.add_argument("file", help="Polarization file (nx ny nz / Ncells / step Px Py Pz)")
    ap.add_argument("--dt-sample-fs", type=float, default=None,
                    help="Time between saved frames (fs). If omitted, only 'frames' are reported.")
    ap.add_argument("--c-window", type=float, default=5.0,
                    help="Sokal window constant c (default 5)")
    ap.add_argument("--maxlen", type=int, default=None,
                    help="Optional: truncate to first maxlen frames to test stability")
    args = ap.parse_args()

    steps, P = load_P(args.file)
    if args.maxlen:
        P = P[:args.maxlen]
    N = len(P)
    print(f"[INFO] Loaded {N} frames.")

    labels = ["Px","Py","Pz"]
    taus_frames = []
    for i, lab in enumerate(labels):
        rho = autocorr_fft(P[:,i])
        tau, W = tau_int_sokal(rho, c=args.c_window)
        taus_frames.append(tau)
        if args.dt_sample_fs:
            tau_ps = tau * args.dt_sample_fs * 1e-3
            print(f"{lab}: tau_int = {tau:.2f} frames  (~{tau_ps:.3f} ps),  window W={W}")
        else:
            print(f"{lab}: tau_int = {tau:.2f} frames,  window W={W}")

    # Conservative block length suggestion:
    # choose max across components, multiply by 10–20 to ensure independence
    tau_max = max(taus_frames)
    block_len_10 = int(np.ceil(10 * tau_max))
    block_len_20 = int(np.ceil(20 * tau_max))
    print("\n[Suggestion] Use a block length between:")
    print(f"  --block-len {block_len_10}  (≈10× max tau_int)")
    print(f"  --block-len {block_len_20}  (≈20× max tau_int)")
    nblocks10 = N // block_len_10 if block_len_10>0 else 0
    nblocks20 = N // block_len_20 if block_len_20>0 else 0
    print(f"  That gives ~{nblocks10} or ~{nblocks20} blocks, respectively.")
    if min(nblocks10, nblocks20) < 5:
        print("  [Note] Consider longer runs or smaller sampling interval to get ≥5–10 blocks.")
    if args.dt_sample_fs:
        bl10_ps = block_len_10 * args.dt_sample_fs * 1e-3
        bl20_ps = block_len_20 * args.dt_sample_fs * 1e-3
        print(f"  (~{bl10_ps:.2f}–{bl20_ps:.2f} ps per block)")

if __name__ == "__main__":
    main()
