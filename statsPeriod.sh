#!/bin/bash

# Configuration
pattern_file="/afs/cern.ch/user/g/gkainth/patterns.txt"
out_dir="/afs/cern.ch/user/g/gkainth/outputsPhast"

# Load patterns (trim empty lines)
mapfile -t patterns < <(grep -v '^[[:space:]]*$' "$pattern_file")
num_patterns=${#patterns[@]}

if (( num_patterns == 0 )); then
  echo "No patterns found in $pattern_file"
  exit 1
fi

first_pattern="${patterns[0]}"

# Collect all .out files
mapfile -t all_files < <(find "$out_dir" -maxdepth 1 -type f -name "*.out" | sort)
if [ ${#all_files[@]} -eq 0 ]; then
  echo "No .out files found in $out_dir"
  exit 1
fi

total_files=${#all_files[@]}
echo "Processing all $total_files files together..."

# Initialize totals arrays (use integers)
declare -a total1
declare -a total2
for ((i=0; i<num_patterns; i++)); do
  total1[i]=0
  total2[i]=0
done

# Progress bar function
print_progress() {
  local current=$1 total=$2
  local width=40
  local filled=$(( current * width / total ))
  local empty=$(( width - filled ))
  printf "\rProgress: ["
  for ((i=0;i<filled;i++)); do printf "#"; done
  for ((i=0;i<empty;i++)); do printf "."; done
  printf "] (%d/%d)" "$current" "$total"
}

# Process all files together
for ((idx=0; idx<total_files; idx++)); do
  file="${all_files[$idx]}"

  # Find the line number of the first occurrence of the first pattern
  match_line=$(grep -nF "$first_pattern" "$file" | head -n 1 | cut -d':' -f1)

  if [[ -n "$match_line" ]]; then
    # Extract the block of lines equal to the number of patterns
    mapfile -t nums < <(sed -n "${match_line},$((match_line + num_patterns - 1))p" "$file" | \
      awk -F':' '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2); print $1, $2}')

    for ((i=0; i<num_patterns; i++)); do
      read -r n1 n2 <<< "${nums[$i]}"
      # Safely add numbers as integers
      total1[$i]=$(( total1[i] + n1 ))
      total2[$i]=$(( total2[i] + n2 ))
    done
  fi

  print_progress $((idx+1)) $total_files
done

echo -e "\nFinished processing all files."

# Write combined output file
out_file="$out_dir/tree-all.txt"
{
  for ((i=0; i<num_patterns; i++)); do
    echo "${total1[i]} : ${total2[i]} : ${patterns[i]}"
  done
} > "$out_file"
echo "Summary written to $out_file"
