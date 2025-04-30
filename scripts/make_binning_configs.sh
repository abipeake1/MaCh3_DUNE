#!/bin/bash
# Define paths
CONFIG_DIR="/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/configs/Samples/OA_Samples/"  # Directory where sample.yaml files are stored
BINNING_DIR="/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/configs/Samples/Binning_Sets"  # Directory containing different binning sets
OUTPUT_DIR="/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/configs/Samples/OA_Samples_projections"  # Directory to save modified YAML files
BACKUP_DIR="/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/configs/Samples/Backup"  # Directory to store backups

# Ensure necessary directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BACKUP_DIR"

# Loop through each binning set
for BINNING_FILE in "$BINNING_DIR"/*.yaml; do
    BINNING_NAME=$(basename "$BINNING_FILE" .yaml)  # Extract binning name
    
    # Loop through each sample.yaml file
    for YAML_FILE in "$CONFIG_DIR"/*.yaml; do
        SAMPLE_NAME=$(basename "$YAML_FILE" .yaml)  # Extract sample name
        OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_${BINNING_NAME}.yaml"
        
        echo "Processing: $YAML_FILE with binning set $BINNING_NAME"
        
        # Create a backup
        cp "$YAML_FILE" "$BACKUP_DIR/$(basename "$YAML_FILE").bak"
        
        # Copy the original file to output directory
        cp "$YAML_FILE" "$OUTPUT_FILE"

        # Format the binning content with proper indentation
        BINNING_CONTENT=$(awk '{print "  " $0}' "$BINNING_FILE")

        # Use awk to replace the Binning section in place
        awk -v binning="$BINNING_CONTENT" '
            BEGIN {in_binning=0}
            /^Binning:/ {print "Binning:"; print binning; print "  ForceGeneric: True"; in_binning=1; next}
            in_binning && /^[^[:space:]]/ {in_binning=0}
            !in_binning {print}
        ' "$YAML_FILE" > "$OUTPUT_FILE"

        echo "Created '$OUTPUT_FILE' with binning set '$BINNING_NAME'"
    done
done

echo "All YAML files have been processed with different binning sets."

