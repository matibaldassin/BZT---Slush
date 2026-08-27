#!/usr/bin/env python3
import re, sys, argparse, numpy as np
from math import cos, sin, radians, sqrt
EPS0 = 8.8541878128e-12
KB   = 1.380649e-23

def volume_from_lattice(a,b,c,alpha=90.0,beta=90.0,gamma=90.0):
    ca, cb, cg = cos(radians(alpha)), cos(radians(beta)), cos(radians(gamma))
    return a*b*c*sqrt(1 + 2*ca*cb*cg - ca*ca - cb*cb - cg*cg)

# NEW: read 3×3 lattice (Å) and return |det| in Å^3
def volume_from_matrix_file(path):
    M = np.loadtxt(path, dtype=float).reshape(3,3)
    # Rows are lattice vectors A,B,C in Å (as in your example)
    V_A3 = abs(np.linalg.det(M))
    return V_A3

def read_P_file(path):
    with open(path,'r') as f:
        nx, ny, nz = [int(x) for x in f.readline().split()]
        Ncells = int(f.readline().strip())
        data=[]
        for line in f:
            s=line.strip()
            if not s: continue
            cols=s.split()
            if len(cols)<4: continue
            step=int(float(cols[0]))
            Px,Py,Pz=[float(c) for c in cols[1:4]]
            data.append((step,Px,Py,Pz))
    if not data: raise ValueError(f"No time series in {path}")
    arr=np.array(data,float)
    steps=arr[:,0].astype(int)
    P_mCcma2=arr[:,1:4]
    return (nx,ny,nz), Ncells, steps, P_mCcma2

def block_cov(P, block_len=None):
    X=np.asarray(P,float); N=len(X)
    if not block_len or block_len<=0 or block_len>N//2:
        C=np.cov(X.T, bias=True)
        return C, np.full((3,3), np.nan)
    nB=N//block_len
    X=X[:nB*block_len].reshape(nB,block_len,3)
    Cs=np.stack([np.cov(X[b].T, bias=True) for b in range(nB)])
    return Cs.mean(0), (Cs.std(0, ddof=1)/np.sqrt(nB) if nB>1 else np.full((3,3), np.nan))

def approx_eps_parallel_err(u, eps_err):
    """
    Approximate std error for eps_parallel = u^T eps u
    assuming elementwise independence: var = sum_{ij} (u_i u_j)^2 Var(eps_ij)
    eps_err contains std errors for each component.
    """
    if not np.isfinite(eps_err).any():
        return np.nan
    u = np.asarray(u, float)
    U = np.outer(u,u)  # U_ij = u_i u_j
    return np.sqrt(np.sum((U**2) * (eps_err**2)))

def main():
    ap=argparse.ArgumentParser(description="Permittivity from polarization fluctuations (with longitudinal/transverse projection)")
    ap.add_argument("files", nargs="+", help="Polarization files")
    ap.add_argument("-T","--temperature", type=float, help="Temperature in K (single-file mode)")
    ap.add_argument("--temp-from-name", type=str, help=r"Regex capturing T in K, e.g. 'T(\d+)K'")
    # ---- VOLUME OPTIONS (choose one family) ----
    ap.add_argument("--lat", action="append",
                    help="Lattice matrix file (Å) per polarization file; give once per input")
    ap.add_argument("--lat-template", type=str,
                    help="Template to build lattice filename from T, e.g. 'lat_T{T}K.dat'")
    ap.add_argument("--supercell-vol-A3", type=float, help="Supercell volume in Å^3 (constant)")
    ap.add_argument("--ucell-vol-A3", type=float, help="Unit-cell volume in Å^3 (multiplies nx*ny*nz)")
    ap.add_argument("--ucell-a", type=float); ap.add_argument("--ucell-b", type=float); ap.add_argument("--ucell-c", type=float)
    ap.add_argument("--alpha", type=float, default=90.0); ap.add_argument("--beta", type=float, default=90.0); ap.add_argument("--gamma", type=float, default=90.0)
    # --------------------------------------------
    ap.add_argument("--block-len", type=int, default=0, help="Frames per block for error bars (0=off)")
    ap.add_argument("--out", default="epsilon_vs_T.csv")
    args=ap.parse_args()

    # Enforce the requested filename
    args.out = "epsilon_vs_T.csv"

    temp_re = re.compile(args.temp_from_name) if args.temp_from_name else None
    if args.lat and len(args.lat)!=len(args.files):
        sys.exit("[ERROR] Number of --lat files must match number of polarization files.")

    rows=[]
    for i, path in enumerate(args.files):
        (nx,ny,nz), Ncells, steps, P_mCcma2 = read_P_file(path)

        # Temperature
        if temp_re:
            m=temp_re.search(path)
            if not m: sys.exit(f"[ERROR] Could not extract T from '{path}'")
            T_K=float(m.group(1))
        elif args.temperature is not None:
            T_K=args.temperature
        else:
            sys.exit(f"[ERROR] Provide -T or --temp-from-name for '{path}'")

        # Volume (Å^3) priority: per-T lattice -> constant supercell -> unit cell info
        if args.lat:
            V_A3 = volume_from_matrix_file(args.lat[i])
        elif args.lat_template:
            if not temp_re:
                sys.exit("[ERROR] --lat-template requires --temp-from-name")
            V_A3 = volume_from_matrix_file(args.lat_template.format(T=int(T_K)))
        elif args.supercell_vol_A3 is not None:
            V_A3 = args.supercell_vol_A3
        else:
            if args.ucell_vol_A3 is not None:
                Vcell_A3 = args.ucell_vol_A3
            else:
                if None in (args.ucell_a, args.ucell_b, args.ucell_c):
                    sys.exit("[ERROR] Provide volume via --lat/--lat-template OR --supercell-vol-A3 OR ucell options.")
                Vcell_A3 = volume_from_lattice(args.ucell_a,args.ucell_b,args.ucell_c,args.alpha,args.beta,args.gamma)
            V_A3 = Vcell_A3 * (nx*ny*nz)
        V_m3 = V_A3 * 1e-30

        # Units: mC/cm^2 -> C/m^2
        P_Cm2 = P_mCcma2 / 100.0

        # Covariance and permittivity
        C, Cerr = block_cov(P_Cm2, args.block_len if args.block_len>0 else None)
        factor = V_m3/(EPS0*KB*T_K)
        chi = factor*C
        eps = np.eye(3)+chi
        eps = 0.5*(eps + eps.T)  # enforce symmetry numerically
        eps_err = (factor*Cerr) if np.isfinite(Cerr).any() else np.full((3,3), np.nan)

        # --- NEW: Longitudinal/transverse projection wrt mean polarization ---
        muP = P_Cm2.mean(axis=0)  # <P> in C/m^2
        u_tol = 1e-9  # C/m^2, threshold to consider <P> ~ 0
        norm_mu = float(np.linalg.norm(muP))
        if norm_mu > u_tol:
            u = muP / norm_mu
        else:
            # Fallback: use principal eigenvector of eps (largest eigenvalue)
            w, V = np.linalg.eigh(eps)
            u = V[:, np.argmax(w)]
            u = u / np.linalg.norm(u)

        eps_parallel = float(u @ eps @ u)
        tr_eps = float(np.trace(eps))
        eps_perp = 0.5*(tr_eps - eps_parallel)

        # Approximate errors for the scalars (if available)
        eps_parallel_err = approx_eps_parallel_err(u, eps_err) if np.isfinite(eps_err).any() else np.nan
        if np.isfinite(eps_err).any():
            var_tr = float(np.sum(np.diag(eps_err)**2))
            var_par = float(eps_parallel_err**2)
            eps_perp_err = 0.5*np.sqrt(var_tr + var_par)  # assumes independence
        else:
            eps_perp_err = np.nan

        rows.append({
            "file": path, "T_K": T_K,
            "eps_xx": eps[0,0], "eps_yy": eps[1,1], "eps_zz": eps[2,2],
            "eps_xy": eps[0,1], "eps_xz": eps[0,2], "eps_yz": eps[1,2],
            "eps_iso": (eps[0,0]+eps[1,1]+eps[2,2])/3.0,
            "eps_xx_err": eps_err[0,0], "eps_yy_err": eps_err[1,1], "eps_zz_err": eps_err[2,2],
            "eps_xy_err": eps_err[0,1], "eps_xz_err": eps_err[0,2], "eps_yz_err": eps_err[1,2],
            "eps_iso_err": (eps_err[0,0]+eps_err[1,1]+eps_err[2,2])/3.0,
            # NEW columns:
            "eps_parallel": eps_parallel,
            "eps_perp": eps_perp,
            "eps_parallel_err": eps_parallel_err,
            "eps_perp_err": eps_perp_err,
            # Optional: store |<P>| for reference (not written unless you add to cols)
            # "P_mean_C_per_m2": norm_mu
        })

    # Write CSV
    import csv
    cols = [
        "file","T_K",
        "eps_xx","eps_yy","eps_zz","eps_xy","eps_xz","eps_yz","eps_iso",
        "eps_xx_err","eps_yy_err","eps_zz_err","eps_xy_err","eps_xz_err","eps_yz_err","eps_iso_err",
        # NEW scalar projections:
        "eps_parallel","eps_perp","eps_parallel_err","eps_perp_err"
    ]
    with open(args.out,"w",newline="") as f:
        w=csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in sorted(rows, key=lambda d:d["T_K"]):
            w.writerow({k:(f"{r[k]:.8g}" if isinstance(r[k],float) else r[k]) for k in cols})
    print(f"[OK] Wrote {args.out} ({len(rows)} temps).")

if __name__=="__main__":
    main()
