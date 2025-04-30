#!/bin/bash

# Set variables
DIRECTORY="$1"  # Directory with YAML files
TOTAL_POT="$2"  # Total POT value
CUSTOM_POT_FILE="$3"  # File to get the custom POT value (e.g., "0m_all.yaml")
CUSTOM_POT_PERCENTAGE="$4"  # Percentage to assign to the custom POT file (e.g., 50 for 50%)

# Ensure the directory exists
if [[ ! -d "$DIRECTORY" ]]; then
    echo "Directory not found: $DIRECTORY"
    exit 1
fi

# Get list of YAML files in the directory, excluding the custom POT file
YAML_FILES=$(find "$DIRECTORY" -name "*.yaml" -not -name "$CUSTOM_POT_FILE")

# Calculate the custom POT value based on the percentage using awk
CUSTOM_POT_VALUE=$(echo "$TOTAL_POT $CUSTOM_POT_PERCENTAGE" | awk '{printf "%.10e", $1 * ($2 / 100)}')

# Calculate the remaining POT value using awk
REMAINING_POT=$(echo "$TOTAL_POT $CUSTOM_POT_VALUE" | awk '{printf "%.10e", $1 - $2}')

# Calculate the POT per file for the other files
NUM_OTHER_FILES=$(echo "$YAML_FILES" | wc -l)
if [[ $NUM_OTHER_FILES -eq 0 ]]; then
    echo "No other YAML files found to distribute the POT."
    exit 1
fi

POT_PER_FILE=$(echo "$REMAINING_POT $NUM_OTHER_FILES" | awk '{printf "%.10e", $1 / $2}')

# Function to update the POT value in a YAML file
update_pot() {
    FILE="$1"
    NEW_POT="$2"
    # Ensure to round the POT value to a reasonable number of significant digits
    FORMATTED_POT=$(printf "%.2e" "$NEW_POT")
    sed -i "s/POT: .*/POT: $FORMATTED_POT/" "$FILE"
}

# Update the custom POT file with the calculated value
update_pot "$DIRECTORY/$CUSTOM_POT_FILE" "$CUSTOM_POT_VALUE"

# Update the remaining YAML files with the evenly distributed POT
for FILE in $YAML_FILES; do
    update_pot "$FILE" "$POT_PER_FILE"
done

echo "POT values updated successfully."
