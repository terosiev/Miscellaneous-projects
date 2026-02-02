#!/bin/bash
################################################################################
# Master Script: Run Complete GWAS Pipeline
#
# Executes all analysis steps in sequence for specified models
#
# Usage: 
#   ./run_pipeline.sh [model_number]
#   ./run_pipeline.sh all  # Run all models
################################################################################

# Source configuration
source scripts/00_config.sh

# Print configuration
print_config

# Parse arguments
if [ $# -eq 0 ]; then
    echo "Usage: $0 <model_number> or $0 all"
    exit 1
fi

MODEL_ARG=$1

# Function to run pipeline for a single model
run_single_model() {
    local model=$1
    
    echo ""
    echo "###############################################"
    echo "# Running GWAS Pipeline for Model ${model}"
    echo "###############################################"
    echo ""
    
    # Create directories
    create_directories $model
    
    # Step 1: Quality Control
    echo "==> Step 1: Quality Control"
    bash scripts/01_qc.sh $model
    if [ $? -ne 0 ]; then
        echo "ERROR: QC failed for model ${model}"
        return 1
    fi
    
    # Step 2: PCA
    echo ""
    echo "==> Step 2: Principal Component Analysis"
    bash scripts/02_pca.sh $model
    if [ $? -ne 0 ]; then
        echo "ERROR: PCA failed for model ${model}"
        return 1
    fi
    
    # Step 3: Association Analysis
    echo ""
    echo "==> Step 3: Association Analysis"
    bash scripts/03_association.sh $model
    if [ $? -ne 0 ]; then
        echo "ERROR: Association analysis failed for model ${model}"
        return 1
    fi
    
    # Step 4: Clumping
    echo ""
    echo "==> Step 4: LD Clumping"
    bash scripts/04_clumping.sh $model
    if [ $? -ne 0 ]; then
        echo "ERROR: Clumping failed for model ${model}"
        return 1
    fi
    
    # Step 5: VEP Annotation
    echo ""
    echo "==> Step 5: VEP Annotation"
    bash scripts/05_vep.sh $model
    if [ $? -ne 0 ]; then
        echo "WARNING: VEP annotation failed for model ${model}"
        echo "Continuing anyway..."
    fi
    
    echo ""
    echo "###############################################"
    echo "# Model ${model} Pipeline Complete!"
    echo "###############################################"
    echo ""
}

# Main execution
start_time=$(date +%s)

if [ "$MODEL_ARG" == "all" ]; then
    echo "Running pipeline for all ${N_MODELS} models..."
    for i in $(seq 1 $N_MODELS); do
        run_single_model $i
    done
else
    echo "Running pipeline for model ${MODEL_ARG}..."
    run_single_model $MODEL_ARG
fi

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo ""
echo "###############################################"
echo "# GWAS Pipeline Execution Complete!"
echo "# Total runtime: $((runtime / 60)) minutes"
echo "###############################################"
