import ROOT
import csv

#*************************************
def root_to_csv(input_root_file, tree_name, output_csv):
    # Open ROOT file and get TTree
    file = ROOT.TFile.Open(input_root_file)
    tree = file.Get(tree_name)

    # Collect all branch names
    branch_list = [b.GetName() for b in tree.GetListOfBranches()]
    scalar_branches = []
    TL_branches = []
    vector_branches = [] 

    for bname in branch_list:
        if bname.endswith("_TL"):
            TL_branches.append(bname)
        elif bname.endswith("_vec"):
            vector_branches.append(bname)
        else:
            scalar_branches.append(bname)

    print(f"Exporting {len(scalar_branches)} scalar branches.")
    print(f"Exporting {len(TL_branches)} TLorentz branches.")
    print(f"Exporting {len(vector_branches)} vector branches.")

    # Prepare output CSV
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # CSV header: scalar branches + TL components + 3-vector components
        headers = scalar_branches.copy()

        for vname in TL_branches:
            if vname in ["cluster_TL", "clusterLow_TL"]:
                headers.extend([f"{vname}_X", f"{vname}_Y", f"{vname}_Z", f"{vname}_E"])
            else:
                headers.extend([f"{vname}_Px", f"{vname}_Py", f"{vname}_Pz", f"{vname}_E"])

        for vname in vector_branches:
            if vname == "clusterFit_vec":
                headers.extend([f"{vname}_X", f"{vname}_Y", f"{vname}_E"])  # Z holds E
            else:
                headers.extend([f"{vname}_X", f"{vname}_Y", f"{vname}_Z"])

        writer.writerow(headers)

        n_entries = tree.GetEntries()
        print(f"Processing {n_entries} entries...")

        for entry in range(n_entries):
            tree.GetEntry(entry)
            row = []

            # Scalars
            for bname in scalar_branches:
                val = getattr(tree, bname)
                try:
                    if isinstance(val, float) and abs(val - int(val)) < 1e-6:
                        row.append(int(val))
                    else:
                        row.append(val)
                except Exception:
                    row.append(val)

            # TL Vectors
            for vname in TL_branches:
                vec = getattr(tree, vname)
                if vname in ["cluster_TL", "clusterLow_TL"]:
                    row.extend([vec.X(), vec.Y(), vec.Z(), vec.T()])  # Position + Energy
                else:
                    row.extend([vec.Px(), vec.Py(), vec.Pz(), vec.E()])  # Standard TLorentzVector

            # 3-Vectors
            for vname in vector_branches:
                vec = getattr(tree, vname)
                if vname == "clusterFit_vec":
                    row.extend([vec.X(), vec.Y(), vec.Z()])  # Z holds energy
                else:
                    row.extend([vec.X(), vec.Y(), vec.Z()])

            writer.writerow(row)

    print(f"CSV export complete: {output_csv}")

#*************************************
input_root_file = "/Users/gursimran/Desktop/filtered_P09_USR970_717.root"
tree_name = "USR970_filtered"
output_csv = "filtered_P09_USR970_717.csv"

root_to_csv(input_root_file, tree_name, output_csv)