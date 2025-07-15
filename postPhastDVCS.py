import argparse
import glob
import os
import ROOT

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
output_file = f"filtered_P{period}_{tree_name}.root"
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

written_keys = set()  # Track written 4-momentum combinations for deduplication

for event in in_tree:
    if is_real_data:
        if (event.Run, event.Spill) in bad_spills:
            continue
        if event.TiS_flag == False:
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

    # -*THIS IS NOT A CUT*- Deduplication: only if low_calo == 2 
    if event.low_calo == 2:
        key = (
            round(event.inMu_TL.Px(), 6),
            round(event.outMu_TL.Px(), 6),
            round(event.gamma_TL.Px(), 6),
            round(event.p_camera_TL.Px(), 6)
        )
        if key in written_keys:
            continue  # Already wrote this event — skip
        written_keys.add(key)

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

    # Passed all cuts — write event
    out_tree.Fill()
    n_kept += 1

print(f"Events kept after filtering: {n_kept} / {n_total} ...")
print(f"Filtered data saved to {output_file} with tree name '{output_tree_name}'")

# Write and close
out_tree.SetName(output_tree_name)
out_tree.Write()
out_file.Close()
in_file.Close()

