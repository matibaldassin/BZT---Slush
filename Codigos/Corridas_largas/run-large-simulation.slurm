#!/bin/bash
#SBATCH --job-name=NAME_OF_THE_JOB    	# Job name
#SBATCH --nodes=1                 	# Number of nodes requested
#SBATCH --ntasks-per-node=48     	# Number of tasks (cores) per node
#SBATCH --time=240:00:00           	# Time limit hrs:min:sec
#SBATCH --output=salida_%j.out  	# Standard output file name
#SBATCH --error=salida_%j.err   	# Standard error file name
#SBATCH --partition=ib100  		# Specify the partition name

module use /opt/ohpc/pub/apps/modules/all/
module load intel/2022a

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

# ---------------- PARÁMETROS ----------------
TARGET_TEMP=70
TOTAL_STEPS=100000
BLOCK_SIZE=100000

WRITE_FREQ=10000
WRITE_START=10000

CHECKPOINT_FILE="checkpoint_steps.txt"

# ---------------- INIT ----------------
if [[ -f $CHECKPOINT_FILE ]]; then
    read -r steps_done < "$CHECKPOINT_FILE"
    echo "Reanudando desde $steps_done pasos"
else
    steps_done=0
    echo "$steps_done" > "$CHECKPOINT_FILE"
    echo "Comenzando desde 0 pasos"
fi

block_id=$((steps_done / BLOCK_SIZE))

# ---------------- LOOP ----------------
while (( steps_done < TOTAL_STEPS )); do

    remaining=$((TOTAL_STEPS - steps_done))

    if (( remaining < BLOCK_SIZE )); then
        current_block_size=$remaining
    else
        current_block_size=$BLOCK_SIZE
    fi

    block_id=$((block_id + 1))
    printf -v block_name "block-%04d" "$block_id"

    echo "----------------------------------------"
    echo "Bloque $block_id"
    echo "Pasos actuales: $steps_done"
    echo "Corriendo: $current_block_size pasos"
    echo "----------------------------------------"

    # -------- TRAJ --------
    traj_line="traj $WRITE_START,$WRITE_FREQ,2"

    echo "Usando: $traj_line"

    # -------- CONTROL --------
    cat > CONTROL << EOF
BNO  bulk
temperature          $TARGET_TEMP
pressure               0
steps             $current_block_size
timestep                0.0004 ps
print every 10000
stats every 10000
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
$traj_line
finish
EOF

    # ------ LIMPIEZA -------
    rm -f P-local P-total

    # -------- RUN --------
    mpirun -np "$SLURM_NTASKS" /home/mbaldassin.ifir/dlpoly/poly-wolf.x
    if [[ $? -ne 0 ]]; then
        echo "ERROR: DL_POLY falló en bloque $block_id"
        exit 1
    fi

    # -------- SAVE --------
    mkdir -p "$block_name"

    cp CONFIG "$block_name"/
    cp REVCON "$block_name"/ || true
    cp FIELD "$block_name"/ || true

    mv -f STATIS "$block_name"/ 2>/dev/null || true
    mv -f OUTPUT "$block_name"/ 2>/dev/null || true
    mv -f REVIVE "$block_name"/ 2>/dev/null || true
    mv -f RDFDAT "$block_name"/ 2>/dev/null || true

    # ---- P-LOCAL por bloque ----
    if [[ -f P-local ]]; then
        mv -f P-local "$block_name"/
    else
        echo "WARNING: P-local no generado en bloque $block_id"
    fi

    # ---- P-TOTAL por bloque ----
    if [[ -f P-total ]]; then
        mv -f P-total "$block_name"/
    else
        echo "WARNING: P-total no generado en bloque $block_id"
    fi
    
    # -------- CONTINUIDAD --------
    if [[ -f REVCON ]]; then
        mv -f REVCON CONFIG
    else
        echo "ERROR: REVCON no generado"
        exit 1
    fi

    # -------- UPDATE CHECKPOINT --------
    steps_done=$((steps_done + current_block_size))
    echo "$steps_done" > "$CHECKPOINT_FILE"

    echo "Checkpoint actualizado: $steps_done pasos"

done

echo "SIMULACIÓN COMPLETADA"
