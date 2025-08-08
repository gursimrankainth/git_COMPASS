import argparse
import glob
import os
import ROOT
import numpy as np
from collections import defaultdict

ROOT.ROOT.EnableImplicitMT()

# ********************************
# Parse arguments
parser = argparse.ArgumentParser(description="Filter tree to select only exclusive photon events (DVCS/BH).")
parser.add_argument(
    "--data", type=str, required=True, choices=["real", "MC"],
    help="Specify whether the input data is 'real' or 'MC'"
)
parser.add_argument("--period", type=str, required=False, default=None,
                    help="Period string (e.g. 09)")
args = parser.parse_args()

tree_name = "USR970"
period = args.period if args.period else input("Enter the period (e.g. 09): ").strip()

input_file = f"merged_P{period}.root"
output_file = f"filtered_P{period}.root"
output_tree_name = "USR970_filtered"

# Determine if real data (set using user input)
is_real_data = args.data == "real"

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
            if flux == 1 or LT == 1 or MT == 1 or OT == 1 or ECAL == 1 or empty == 1:
            #if flux == 1 or MT == 1 or LT == 1 or RICH == 1 or ECAL == 1:
                bad_spills.add((run, spill))

print(f"Loaded {len(bad_spills)} bad spills ...")

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
# Identify events with multiplicty greater than one
multiplicity = defaultdict(set)
for i in range(in_tree.GetEntries()):
    in_tree.GetEntry(i)
    event_id = in_tree.Evt
    comb = (
        in_tree.inMu_TL.Px(),
        in_tree.outMu_TL.Px(),
        in_tree.gamma_TL.Px(),
        in_tree.p_camera_TL.Px()
    )
    multiplicity[event_id].add(comb)

events_with_multiple_combs = {evt for evt, combs in multiplicity.items() if len(combs) > 1}
print(f"Found {len(events_with_multiple_combs)} events with multiple unique particle combinations ...")

# ********************************
# Apply cuts and save the passing events 
out_file = ROOT.TFile.Open(output_file, "RECREATE")
out_tree = in_tree.CloneTree(0)

n_total = in_tree.GetEntries()
n_kept = 0

for event in in_tree:
    if is_real_data:
        if (event.Run, event.Spill) in bad_spills:
            continue
        if not event.TiS_flag:
            continue 

    # Passed hodoscope and trigger checks 
    if not (event.trig_MT or event.trig_OT or event.trig_LT):
        continue 
    if not event.hodoPass:
        continue 

    # Exclusivity variables 
    if np.abs(event.delta_pt) > 0.3:
        continue 
    delta_phi = event.delta_phi
    if not ((-0.4 <= delta_phi <= 0.4) or
            (-0.4 <= delta_phi + 2 * np.pi <= 0.4) or
            (-0.4 <= delta_phi - 2 * np.pi <= 0.4)):
        continue
    if np.abs(event.delta_Z) > 16:
        continue 
    if np.abs(event.M2x) > 0.3:
        continue 

    # Remove events with multiplicity > 1 
    if event.Evt in events_with_multiple_combs: 
        continue 

    # Remove visible pi0 events 
    if (event.low_calo == 0 and event.clusterLow_TL.T() > 0.5) or (event.low_calo == 1 and event.clusterLow_TL.T() > 0.63):
        pi0_TL = event.gamma_TL + event.gammaLow_TL
        pi0_M = pi0_TL.M()
        if (0.11 < pi0_M < 0.155):
            continue
       
    # Kinematic fit 
    if not event.fit_conv:
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
    if not chi2_red < 10:
        continue
    
    # Fill the tree with passing events 
    out_tree.Fill()
    n_kept += 1

print(f"Events kept after filtering: {n_kept} / {n_total}")
print(f"Filtered events saved to: {output_file} ({n_kept} entries)")

# Write and close
out_tree.SetName(output_tree_name)
out_tree.Write()
out_file.Close()
in_file.Close()
