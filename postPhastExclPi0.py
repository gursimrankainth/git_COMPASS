import argparse
import glob
import os
import ROOT
from collections import defaultdict

ROOT.ROOT.EnableImplicitMT()

# ********************************
# Parse arguments
parser = argparse.ArgumentParser(description="Filter DVCS tree with bad spill, multiplicity, and duplicate event cuts.")
parser.add_argument("--tree", type=str, required=False, default=None,
                    help="Tree name (e.g. USR970 for real data, USR970_MC for MC)")
parser.add_argument("--period", type=str, required=False, default=None,
                    help="Period string (e.g. 09)")
args = parser.parse_args()

tree_name = args.tree if args.tree else input("Enter the tree name (USR970 for real data, USR970_MC for MC): ").strip()
period = args.period if args.period else input("Enter the period (e.g. 09): ").strip()

input_file = f"merged_P{period}_{tree_name}.root"
output_file = f"filtered_pi0_P{period}_{tree_name}.root"
intermediate_tree_name = f"{tree_name}_with_duplicates"
output_tree_name = f"{tree_name}_filtered"

# Determine if real data (assuming 'MC' in tree name means MC)
is_real_data = "MC" not in tree_name

# ********************************
# Load bad spill list if real data
bad_spills = set()
if is_real_data:
    bad_spill_dir = "/afs/cern.ch/user/g/gkainth/phastPackages/bad_spill"
    pattern = os.path.join(bad_spill_dir, f"P{period}*bad_spill.lst")
    matching_files = glob.glob(pattern)
    if len(matching_files) == 0:
        raise FileNotFoundError(f"No bad spill file found for period P{period} in {bad_spill_dir}")
    elif len(matching_files) > 1:
        print(f"Warning: Multiple bad spill files found for period P{period}. Using the first one: {matching_files[0]}")

    bad_spill_file = matching_files[0]
    print(f"Using bad spill file: {bad_spill_file} ... ")

    with open(bad_spill_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            run, spill = int(parts[0]), int(parts[1])
            flux, LAST, LT, MT, OT, RICH, ECAL, empty = map(int, parts[2:])
            # Add to bad_spills if bad due to flags other than LAST or RICH
            if flux == 1 or LT == 1 or MT == 1 or OT == 1 or ECAL == 1 or empty == 1:
                bad_spills.add((run, spill))

print(f"Loaded {len(bad_spills)} bad spills (excluding LAST and RICH only) ...")

# ********************************
# Open input ROOT file and get tree
in_file = ROOT.TFile.Open(input_file)
if not in_file or in_file.IsZombie():
    raise RuntimeError(f"Could not open input file '{input_file}'")

in_tree = in_file.Get(tree_name)
if not in_tree:
    raise RuntimeError(f"Tree '{tree_name}' not found in file '{input_file}'")

print(f"Processing {in_tree.GetEntries()} entries in tree '{tree_name}'...")

# ********************************
# Phase 1: Find all events with multiplicity > 1
events_with_high_mult = set()
for i in range(in_tree.GetEntries()):
    in_tree.GetEntry(i)
    if in_tree.multiplicity > 1:
        events_with_high_mult.add(in_tree.Evt)

print(f"Found {len(events_with_high_mult)} events with multiplicity > 1 ...")

# ********************************
# Phase 2: Apply cuts and save passing events
out_file = ROOT.TFile.Open(output_file, "RECREATE")
out_tree = in_tree.CloneTree(0)

n_total = in_tree.GetEntries()
n_kept = 0

for event in in_tree:
    if is_real_data:
        if (event.Run, event.Spill) in bad_spills:
            continue
        if event.TiS_flag == False:
            continue

    # Passed Hodoscope and trigger checks 
    if event.trig_flag == False: 
        continue 
    if event.hodoPass == False: 
        continue 

    # Exclusivity variables
    if abs(event.delta_pt) > 0.3:
        continue
    if ROOT.TMath.Abs(ROOT.TVector2.Phi_mpi_pi(event.delta_phi)) > 0.4:
        continue
    if abs(event.delta_Z) > 16:
        continue
    if abs(event.M2x) > 0.3:
        continue

    # Multiplicity
    if event.Evt in events_with_high_mult:
        continue

    # Kinematic fit cuts
    if event.fit_conv == False:
        continue
    if not (1 < event.Q2_fit < 10):
        continue
    if not (0.05 < event.y_fit < 0.95):
        continue
    if not (-0.64 < event.t_fit < -0.08):
        continue
    if not (10 < event.nu_fit < 144):
        continue
    chi2_red = event.chi2_fit / event.ndf_fit
    if chi2_red > 10:
        continue

    # Build pi0 candidates after applying all DVCS cuts
    if event.low_calo == -999 or event.low_calo == 2:
        continue

    if (event.low_calo == 0 and event.clusterLow_TL.T() > 0.5) or (event.low_calo == 1 and event.clusterLow_TL.T() > 0.63):
        pi0_TL = event.gamma_TL + event.gammaLow_TL
        pi0_M = pi0_TL.M()
    else:
        continue 

    # Remove bad pi0 candidates and inclusive events (more than one pi0 is identified)
    if not (0.11 < pi0_M < 0.155):
        continue 

    # Passed all cuts — write to intermediate tree
    out_tree.Fill()
    n_kept += 1

print(f"Events kept after pi0 identification: {n_kept} / {n_total} ...")

# Save intermediate tree (may contain duplicate π⁰ events)
out_tree.SetName(intermediate_tree_name)
out_tree.Write()

# ********************************
# Phase 3: Filter for exclusive π⁰ events (remove events with multiple π⁰ candidates)
if n_kept > 0:
    print("Filtering for exclusive π⁰ events (only one π⁰ per event)...")

    file = ROOT.TFile.Open(output_file, "UPDATE")
    intermediate_tree = file.Get(intermediate_tree_name)

    evt_count = defaultdict(int)

    # First pass: Count how many times each Evt appears
    for entry in intermediate_tree:
        evt_count[entry.Evt] += 1

    # Second pass: Fill only exclusive π⁰ events
    exclusive_tree = intermediate_tree.CloneTree(0)
    for entry in intermediate_tree:
        if evt_count[entry.Evt] == 1:
            exclusive_tree.Fill()

    print(f"Events kept after enforcing π⁰ exclusivity: {exclusive_tree.GetEntries()} / {n_kept}")

    # Save under expected final output name
    exclusive_tree.SetName(output_tree_name)
    exclusive_tree.Write()
    file.Close()
else:
    print("No events passed initial cuts; skipping exclusivity filtering.")

# Close 
in_file.Close()
