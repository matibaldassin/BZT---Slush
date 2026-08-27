#!/bin/bash
#SBATCH --job-name=0.40sol    # Job name
#SBATCH --nodes=1                 # Number of nodes requested
#SBATCH --ntasks-per-node=48      # Number of tasks (cores) per node
#SBATCH --time=48:00:00          # Time limit hrs:min:sec
#SBATCH --output=salida_%j.out    # Standard output file name
#SBATCH --error=salida_%j.err     # Standard error file name
#SBATCH --partition=ib100        # Specify the partition name

# -------- SECTION print some infos to stdout ---------------------------------
echo " "
echo "START_TIME           = $(date +'%y-%m-%d %H:%M:%S %s')"
START_TIME=$(date +%s)
echo "HOSTNAME             = $HOSTNAME"
echo "JOB_NAME             = $SLURM_JOB_NAME"
echo "JOB_ID               = $SLURM_JOB_ID"
echo "SLURM_SUBMIT_DIR     = $SLURM_SUBMIT_DIR"
echo "SLURM_JOB_NODELIST   = $SLURM_JOB_NODELIST"
echo "SLURM_NTASKS         = $SLURM_NTASKS"
echo "Creating TMP_WORK_DIR directory ..."
echo " "

#-------- SECTION environment --------------------------------------------------
module use /opt/ohpc/pub/apps/modules/all/
module load intel/2022a

date

# Safer globs: if no files match, the pattern expands to nothing (not itself)
shopt -s nullglob

run_simulation() {
    local start=$1
    local end=$2
    local step=$3
    local prefix=$4

    local checkpoint_file="checkpoint_${prefix}.txt"
    local done_file="${prefix}.done"


    local p_files=(P-*)
    (( ${#p_files[@]} )) && rm -f "${p_files[@]}"

    # If the stage is already marked as done, skip (defensive; caller also checks)
    if [[ -f $done_file ]]; then
        echo "[${prefix}] Already completed. Skipping."
        return 0
    fi

    # Resume ii from checkpoint if present, else start fresh
    local ii
    if [[ -f $checkpoint_file ]]; then
        read -r ii < "$checkpoint_file"
        echo "[${prefix}] Resuming from checkpoint: ii = $ii"
    else
        ii=$start
        echo "$ii" > "$checkpoint_file"
        echo "[${prefix}] Starting fresh: ii = $ii"
    fi

    # Loop across temperatures (supports ascending or descending)
    while (( (step > 0 && ii <= end) || (step < 0 && ii >= end) )); do
        echo "[${prefix}] La temperatura es $ii K"

        # Generate CONTROL file for this temperature
        cat > CONTROL << EOF
BNO  bulk
temperature          $ii
pressure               0
steps             80000
timestep                0.0004 ps
print every 1000
stats every 1000
cutoff         12.00
primary cutoff 10.00
rdf         10
print rdf
delr       .25
wolf sum 0.30273
rvdw 8.0
ensemble nst ber   0.2000  0.5000
job time   1520000
close time 200
traj 10000,1000,2
finish
EOF

        # Run DL_POLY
        mpirun -np "$SLURM_NTASKS" /home/mbaldassin.ifir/dlpoly/poly-wolf.x
        if [[ $? -ne 0 ]]; then
            echo "[${prefix}] ERROR: DL_POLY fallo en T=${ii}K"
            exit 1
        fi

        # Save outputs in a temperature-tagged folder
        outdir="${prefix}-${ii}"
        mkdir -p "$outdir"

        # Keep CONFIG for next step and archive relevant files
        cp CONFIG "$outdir"/
        if [[ -f REVCON ]]; then
            cp REVCON "$outdir"/
            mv -f REVCON CONFIG
        else
            echo "[${prefix}] ERROR: REVCON no generado en T=${ii}K"
            exit 1
        fi

        mv -f STATIS "$outdir"/ 2>/dev/null || true
        mv -f OUTPUT "$outdir"/ 2>/dev/null || true
        mv -f REVIVE "$outdir"/ 2>/dev/null || true
        mv -f P-* "$outdir"/ 2>/dev/null || true
        mv -f posiciones "$outdir"/ 2>/dev/null || true
        mv -f RDFDAT "$outdir"/ 2>/dev/null || true

        # Advance and update checkpoint
        ii=$((ii + step))
        echo "$ii" > "$checkpoint_file"
    done

    # Stage completed: clean checkpoint and mark done
    rm -f "$checkpoint_file"
    touch "$done_file"
    echo "[${prefix}] Completed. Checkpoint removed and ${done_file} created."
}

maybe_run_simulation() {
    local start=$1 end=$2 step=$3 prefix=$4
    local done_file="${prefix}.done"
    if [[ -f $done_file ]]; then
        echo "[${prefix}] Stage already done. Skipping."
    else
        run_simulation "$start" "$end" "$step" "$prefix"
    fi
}

# ---------- PIPELINE: will resume safely mid-pipeline on restart -------------
#maybe_run_simulation 10  90  10  "heating-1"
#maybe_run_simulation 70  10 -20  "cooling"
maybe_run_simulation 600  10 -10  "Cooling"

