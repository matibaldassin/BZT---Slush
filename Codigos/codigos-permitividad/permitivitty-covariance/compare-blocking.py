#!/usr/bin/env python3
# compare_blocking.py
# Compare permittivity estimates with and without block averaging for ONE file.
# Your file format: line1 nx ny nz; line2 Ncells; then: step Px Py Pz (μC/cm^2)

import argparse, numpy as np

EPS0 = 8.8541878128e-12   # F/m
KB   = 1.380649e-23       # J/K

def read_P(path):
    with open(path,'r') as f:
        nx, ny, nz = map(int, f.readline().split())
        _ = int(f.readline())  # Ncells
        data=[]
        for line in f:
            s=line.strip()
            if not s: continue
            cols=s.split()
            if len(cols) < 4: continue
            step=int(float(cols[0]))
            px,py,pz = map(float, cols[1:4])  # μC/cm^2
            data.append((step,px,py,pz))
    arr = np.array(data, float)
    steps = arr[:,0].astype(int)
    P = arr[:,1:4] * 1e-2   # μC/cm^2 -> C/m^2
    return (nx,ny,nz), steps, P

def vol_from_lattice_file(lat_path):
    # 3 lines, each a lattice vector (Å). Volume = |det| in Å^3.
    M = np.loadtxt(lat_path, dtype=float).reshape(3,3)
    V_A3 = abs(np.linalg.det(M))
    return V_A3 * 1e-30   # -> m^3

def cov_no_block(P):
    # population covariance (bias=True)
    return np.cov(P.T, bias=True)

def cov_block(P, block_len):
    N=len(P)
    nB = N // block_len
    X  = P[:nB*block_len].reshape(nB, block_len, 3)
    Cs = np.stack([np.cov(X[b].T, bias=True) for b in range(nB)])
    Cmean = Cs.mean(0)
    # standard error of covariance estimate across blocks
    Cstderr = Cs.std(0, ddof=1)/np.sqrt(nB) if nB>1 else np.full((3,3), np.nan)
    return Cmean, Cstderr, nB

def eps_from_cov(C, T_K, V_m3):
    chi = (V_m3/(EPS0*KB*T_K)) * C
    eps = np.eye(3) + chi
    return eps

def main():
    ap = argparse.ArgumentParser(description="Compare ε with/without block averaging")
    ap.add_argument("pol_file", help="polarization file (one temperature)")
    ap.add_argument("--lat", required=True, help="lattice 3x3 file (Å) for same temperature")
    ap.add_argument("-T","--temperature", type=float, required=True, help="temperature in K")
    ap.add_argument("--block-len", type=int, default=6, help="frames per block (e.g., 6)")
    args = ap.parse_args()

    (nx,ny,nz), steps, P = read_P(args.pol_file)
    V = vol_from_lattice_file(args.lat)
    T = args.temperature

    # --- No blocking ---
    C0 = cov_no_block(P)
    eps0 = eps_from_cov(C0, T, V)

    # crude naive stderr: treat frames independent => Var(mean) ~ Var/ N.
    # propagate linearly to χ and ε; for a quick scalar comparison we show ε_iso only.
    N = len(P)
    varP = P.var(axis=0, ddof=0)             # per component variance
    # covariance of components' means ~ C0 / N
    C0_stderr = C0 / np.sqrt(N)              # rough (not exact for all entries)
    eps0_stderr = (V/(EPS0*KB*T)) * C0_stderr

    # --- With blocking ---
    if args.block_len <= 0 or args.block_len > len(P)//2:
        raise SystemExit("Choose a --block-len between 1 and N/2 frames.")
    Cb, Cb_se, nB = cov_block(P, args.block_len)
    epsb = eps_from_cov(Cb, T, V)
    epsb_se = (V/(EPS0*KB*T)) * Cb_se

    # helpers
    def iso(M): return np.trace(M)/3.0
    def fmt3x3(M):
        return f"[[{M[0,0]:8.2f} {M[0,1]:8.2f} {M[0,2]:8.2f}]\n [{M[1,0]:8.2f} {M[1,1]:8.2f} {M[1,2]:8.2f}]\n [{M[2,0]:8.2f} {M[2,1]:8.2f} {M[2,2]:8.2f}]]"

    print(f"\nFrames: {len(P)} | Blocks (len={args.block_len}): {nB}")
    print("\n--- NO BLOCKING ---")
    print("epsilon_r tensor:\n", fmt3x3(eps0))
    print(f"epsilon_iso (no block): {iso(eps0):.2f}")
    print("\n--- BLOCK AVERAGING ---")
    print("epsilon_r tensor (block mean):\n", fmt3x3(epsb))
    print("stderr on tensor (block):\n", fmt3x3(epsb_se))
    print(f"epsilon_iso (block): {iso(epsb):.2f}  ± {iso(epsb_se):.2f}")

    print("\nNote:")
    print("- If the two ε tensors are close but the block stderr is larger than the naive one,")
    print("  it means blocking is guarding you against underestimating uncertainty.")
    print("- Increase --block-len until it’s > ~10× your autocorrelation time (in frames);")
    print("  then increase trajectory length to raise the number of blocks (nB ≥ 20 is great).")

if __name__ == "__main__":
    main()
