#!/bin/bash

# Load variables from config
source config.sh

# Step 1: Read BAMs and create RData
jid1=$(sbatch --parsable scripts/decon/00_launcher_DECON_ReadInBam.sbatch \
    "$DIR_BAM" "$PANEL" "$FASTA" "$DECON_output")

# Set path to RData file generated in step 1
RDATA_FILE="${DECON_output}/Output_ReadInBams.RData"

# Step 2: Identify sample failures
jid2=$(sbatch --parsable --dependency=afterok:${jid1} scripts/decon/01_launcher_DECON_IdentifyFailures.sbatch \
    "$RDATA_FILE" "$DECON_output")

# Step 3: Call CNVs
jid3=$(sbatch --parsable --dependency=afterok:${jid1} scripts/decon/02_launcher_DECON_MakeCNVCalls.sbatch \
    "$RDATA_FILE" "$DECON_output")

# Print final job ID
echo "$jid3"
