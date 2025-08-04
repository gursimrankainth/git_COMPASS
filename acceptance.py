import ROOT
from ROOT import TLorentzVector
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

#TODO: Plot the distributions for Q2, nu, x and t for the MC data and real data

# **********************************
# Reconstructed MC data after full DVCS selection
file_rec = ROOT.TFile.Open("/afs/cern.ch/user/g/gkainth/filtered_P09.root")
tree_rec = file_rec.Get("USR970_filtered")

# Generated MC data 
file_gen = ROOT.TFile.Open("/afs/cern.ch/user/g/gkainth/merged_P09.root")
tree_gen = file_gen.Get("USR970")

print("Total Entries: Reconstructed:", tree_rec.GetEntries(), ", Generated:", tree_gen.GetEntries())

# **********************************
# Define the bin edges/bins for acceptance calculation 
# nu: 11 bins of width 2 GeV between 10 and 32 GeV
nu_edges = np.linspace(10, 32, 12)
nu_bins = list(zip(nu_edges[:-1], nu_edges[1:]))
# Q2: 9 bins of width 1 (GeV/c)^2 between 1 and 10
Q2_edges = np.linspace(1, 10, 10)
Q2_bins = list(zip(Q2_edges[:-1], Q2_edges[1:]))
# |t|: 4 bins -> each bin should have roughly the same no. of events 
t_edges = [0.08, 0.136, 0.219, 0.36, 0.64]
t_bins = list(zip(t_edges[:-1], t_edges[1:]))
# phi: 8 bins of width pi/4 rad between -pi and pi
phi_edges = np.linspace(-np.pi, np.pi, 9)
phi_bins = list(zip(phi_edges[:-1], phi_edges[1:]))

# **********************************
# Define histograms
hMC_Q2nu = ROOT.TH2F("hMC_Q2nu", "Q^{2}-#nu Distribution - Reconstructed MC; Q^{2} [(GeV/c)^{2}]; #nu [GeV]", 100, 0, 11, 100, 0, 35)
hMC_Q2t = ROOT.TH2F("hMC_Q2t", "Q^{2}-|t| Distribution - Reconstructed MC; Q^{2} [(GeV/c)^{2}]; |t| [(GeV/c)^{2}]", 100, 0, 11, 100, 0, 1)
hMC_phiRV = ROOT.TH1F("hMC_phiRV", "#phi_{#gamma^{*}#gamma} Distribution - Reconstructed MC; #phi_{#gamma^{*}#gamma} [rad]; Counts", 100, -3.2, 3.2)

# Loop over events and fill histograms
for event in tree_rec: 
    Q2 = event.Q2
    nu = event.nu 
    t = event.t

    inMu = event.inMu_TL
    outMu = event.outMu_TL
    proton = event.p_camera_TL
    gamma = event.gamma_TL
    phi = event.phi_gg

    hMC_Q2nu.Fill(Q2,nu)
    hMC_Q2t.Fill(Q2,abs(t))
    hMC_phiRV.Fill(phi)

def draw_bin_lines_2D(hist, x_edges, y_edges):
    lines = []  # store lines to keep them alive

    # Vertical lines (x-axis bin edges)
    for x in x_edges:
        line = ROOT.TLine(x, y_edges[0], x, y_edges[-1])
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(1)
        line.Draw("same")
        lines.append(line)

    # Horizontal lines (y-axis bin edges)
    for y in y_edges:
        line = ROOT.TLine(x_edges[0], y, x_edges[-1], y)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(1)
        line.Draw("same")
        lines.append(line)

    return lines

# Create canvases to draw and save the histograms
""" c1 = ROOT.TCanvas("c1", "Q2 vs nu", 800, 600)
hMC_Q2nu.Draw("COLZ") 
lines1 = draw_bin_lines_2D(hMC_Q2nu, Q2_edges, nu_edges)
c1.SaveAs("Q2_vs_nu.png")

c2 = ROOT.TCanvas("c2", "Q2 vs |t|", 800, 600)
#c2.SetLogy()
hMC_Q2t.Draw("COLZ")
lines2 = draw_bin_lines_2D(hMC_Q2t, Q2_edges, t_edges)
c2.SaveAs("Q2_vs_t.png")

c3 = ROOT.TCanvas("c3", "Phig*g", 800, 600)
hMC_phiRV.Draw()
c3.SaveAs("phiRV.png") """

# **********************************
# Split the data by mu+/mu- to calculate the acceptance 
# Reconstructed data
tree_rec.Draw(">>elist_muPlus_rec", "Q_beam > 0", "entrylist")
elist_muPlus_rec = ROOT.gDirectory.Get("elist_muPlus_rec")
tree_rec.Draw(">>elist_muMinus_rec", "Q_beam < 0", "entrylist")
elist_muMinus_rec = ROOT.gDirectory.Get("elist_muMinus_rec")

# Generated data 
tree_gen.Draw(">>elist_muPlus_gen", "Q_beam > 0", "entrylist")
elist_muPlus_gen = ROOT.gDirectory.Get("elist_muPlus_gen")
tree_gen.Draw(">>elist_muMinus_gen", "Q_beam < 0", "entrylist")
elist_muMinus_gen = ROOT.gDirectory.Get("elist_muMinus_gen")

print("Mu+ events: Recsontrcuted:", elist_muPlus_rec.GetN(), ", Generated:",elist_muPlus_gen.GetN())
print("Mu- events: Recsontrcuted:", elist_muMinus_rec.GetN(), ", Generated:",elist_muMinus_gen.GetN())

# **********************************
# Define arrays to store the acceptance/weighted sum values for each phi bin (integrate over t)
# *** mu+ ***
rec_weights_muPlus = np.zeros((11,9,8)) # array that holds the weighted sums for reconstructed data 
gen_weights_muPlus = np.zeros((11,9,8)) # array that holds the weighted sums for generated data 
acceptance_muPlus  = np.zeros((11,9,8)) # array that stores the acceptance values 

acceptance_err_muPlus = np.zeros((11,9,8)) # array that stores the estimated statistical errors 
rec_weight2_muPlus = np.zeros((11,9,8))  # holds variance for each bin (reconstructed) 
gen_weight2_muPlus = np.zeros((11,9,8))  # holds variance for each bin (generated)

# *** mu- ***
rec_weights_muMinus = np.zeros((11,9,8)) # array that holds the weighted sums for reconstructed data 
gen_weights_muMinus = np.zeros((11,9,8)) # array that holds the weighted sums for generated data 
acceptance_muMinus  = np.zeros((11,9,8)) # array that stores the acceptance values 

acceptance_err_muMinus = np.zeros((11,9,8)) # array that stores the estimated statistical errors
rec_weight2_muMinus = np.zeros((11,9,8))  # holds variance for each bin (reconstructed) 
gen_weight2_muMinus = np.zeros((11,9,8))  # holds variance for each bin (generated)

# **********************************
# Function for filling the weight arrays (integrating over t)
def fill_weights_int_t(tree, elist, charge, data_type):
    """
    tree      : ROOT TTree to loop over
    elist     : ROOT TEntryList for selecting entries
    charge    : string, either "muPlus" or "muMinus"
    data_type : string, either "rec" or "gen"
    """
    tree.SetEntryList(elist)
    for event in tree:
        Q2 = event.Q2
        nu = event.nu
        phi = event.phi_gg
        weight_BH = event.weight_BH

        # Determine bins
        nu_bin = int(math.floor((nu - 10.) / 2.))
        Q2_bin = int(math.floor((Q2 - 1.) / 1.))
        phi_bin = int(math.floor((phi + np.pi) / (np.pi/4)))

        if (nu_bin >= len(nu_bins)) or (Q2_bin >= len(Q2_bins)) or (phi_bin >= len(phi_bins)):
            continue

        if data_type == "rec":
            if charge == "muPlus":
                rec_weights_muPlus[nu_bin][Q2_bin][phi_bin] += weight_BH
                rec_weight2_muPlus[nu_bin][Q2_bin][phi_bin] += weight_BH ** 2 
            elif charge == "muMinus":
                rec_weights_muMinus[nu_bin][Q2_bin][phi_bin] += weight_BH
                rec_weight2_muMinus[nu_bin][Q2_bin][phi_bin] += weight_BH ** 2 
        elif data_type == "gen":
            if charge == "muPlus":
                gen_weights_muPlus[nu_bin][Q2_bin][phi_bin] += weight_BH
                gen_weight2_muPlus[nu_bin][Q2_bin][phi_bin] += weight_BH ** 2 
            elif charge == "muMinus":
                gen_weights_muMinus[nu_bin][Q2_bin][phi_bin] += weight_BH
                gen_weight2_muMinus[nu_bin][Q2_bin][phi_bin] += weight_BH ** 2 

# **********************************
# Function for finding the statistical errors (integrating over t)
def acc_error_int_t(rec_weights, gen_weights, rec_var, gen_var, charge):
    for i in range(11):
        for j in range(9):
            for k in range(8):
                weightRec = rec_weights[i][j][k]
                weightGen = gen_weights[i][j][k]
                varRec = rec_var[i][j][k]
                varGen = gen_var[i][j][k]

                if weightGen == 0:
                    errAcc = 0.0
                else:
                    errT1 = (1 / weightGen)**2 * varRec
                    errT2 = (weightRec / (weightGen ** 2))**2 * varGen
                    errAcc = np.sqrt(errT1 + errT2)

                if charge == "muPlus":
                    acceptance_err_muPlus[i][j][k] = errAcc
                elif charge == "muMinus":
                    acceptance_err_muMinus[i][j][k] = errAcc

# **********************************
# Fill the weight arrays and calculate the acceptance for mu+ and mu- 
fill_weights_int_t(tree_rec, elist_muPlus_rec, "muPlus", "rec")
fill_weights_int_t(tree_rec, elist_muMinus_rec, "muMinus", "rec") 
fill_weights_int_t(tree_gen, elist_muPlus_gen, "muPlus", "gen")
fill_weights_int_t(tree_gen, elist_muMinus_gen, "muMinus", "gen")

# *** mu+ events *** 
nonzero_mask_muPlus = gen_weights_muPlus != 0 # protect against division by zero 
acceptance_muPlus[nonzero_mask_muPlus] = rec_weights_muPlus[nonzero_mask_muPlus] / gen_weights_muPlus[nonzero_mask_muPlus]
acc_error_int_t(rec_weights_muPlus, gen_weights_muPlus, rec_weight2_muPlus, gen_weight2_muPlus, "muPlus")

# *** mu- events ***
nonzero_mask_muMinus = gen_weights_muMinus != 0 # protect against division by zero 
acceptance_muMinus[nonzero_mask_muMinus] = rec_weights_muMinus[nonzero_mask_muMinus] / gen_weights_muMinus[nonzero_mask_muMinus]
acc_error_int_t(rec_weights_muMinus, gen_weights_muMinus, rec_weight2_muMinus, gen_weight2_muMinus, "muMinus")

# **********************************
# Plot the dsitrubtion of the events by bin (integrated over t)
def plot_acceptance_integrated_t():
    phi_bin_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    n_nu_bins, n_Q2_bins, n_phi_bins = acceptance_muPlus.shape

    fig, axes = plt.subplots(nrows=n_nu_bins, ncols=n_Q2_bins, figsize=(18, 22), sharex=True, sharey=True)

    for i in range(n_nu_bins):      # nu bins (y-axis)
        for j in range(n_Q2_bins):  # Q2 bins (x-axis)
            ax = axes[i, j]

            y_muPlus = acceptance_muPlus[i, j]
            yerr_muPlus = acceptance_err_muPlus[i, j]
            y_muMinus = acceptance_muMinus[i, j]
            yerr_muMinus = acceptance_err_muMinus[i, j]

            ax.errorbar(phi_bin_centers, y_muPlus, yerr=yerr_muPlus, fmt='o', color='red', markersize=3, label='μ⁺' if i == 0 and j == 0 else "")
            ax.errorbar(phi_bin_centers, y_muMinus, yerr=yerr_muMinus, fmt='o', color='black', markersize=3, label='μ⁻' if i == 0 and j == 0 else "")
            
            ax.set_ylim(0, 0.8)
            ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
            ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"], fontsize=10)
            ax.grid(True, linestyle='--', linewidth=0.5)
            ax.axhline(0, color='gray', linewidth=0.5)

    # Legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize=14, markerscale=1.5)

    # Create a new set of axes for the phi and acceptance scale at the top right corner
    phi_axis = fig.add_axes([0.817, 0.92, 0.076, 0.034])  # [left, bottom, width, height]
    acc_axis = fig.add_axes([0.8, 0.877, 0.11, 0.064])

    # Setup phi axis
    phi_axis.set_xlim(-np.pi, np.pi)
    phi_axis.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    phi_axis.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"], fontsize=12)
    phi_axis.set_yticks([])
    phi_axis.yaxis.set_visible(False)
    phi_axis.tick_params(axis='x', direction='in', length=5, top=True, bottom=False)
    phi_axis.xaxis.set_label_position('top')
    phi_axis.set_xlabel(r"$\phi_{\gamma^*\gamma}$ [rad]", fontsize=14, labelpad=20)
    phi_axis.xaxis.tick_top()
    phi_axis.patch.set_facecolor('none')
    for name, spine in phi_axis.spines.items():
        spine.set_visible(name == 'top')
        if name == 'top':
            spine.set_linewidth(1.0)
            spine.set_color('black')

    # Setup acceptance axis
    acc_axis.set_ylim(0, 0.8)
    acc_axis.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    acc_axis.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8"], fontsize=12)
    acc_axis.set_xticks([])
    acc_axis.xaxis.set_visible(False)
    acc_axis.set_ylabel("Acceptance", fontsize=14, labelpad=15)
    acc_axis.yaxis.set_label_position('right')
    acc_axis.yaxis.tick_right()
    acc_axis.tick_params(axis='y', direction='in', length=5)
    acc_axis.patch.set_facecolor('none')
    for name, spine in acc_axis.spines.items():
        spine.set_visible(name == 'right')
        if name == 'right':
            spine.set_linewidth(1.0)
            spine.set_color('black')

    # Global axes for Q2 and nu 
    Q2_axis = fig.add_axes([0.116, 0.09, 0.778, 0.035])
    nu_axis = fig.add_axes([0.085, 0.11, 0.035, 0.835])

    # Q2 axis setup
    Q2_axis.set_xlim(Q2_edges[0], Q2_edges[-1])
    Q2_axis.set_xticks(Q2_edges)
    Q2_axis.set_xticklabels([f"{int(a)}" for a in Q2_edges], fontsize=12)
    Q2_axis.set_yticks([])
    Q2_axis.xaxis.tick_bottom()
    Q2_axis.set_xlabel(r"$Q^2$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
    Q2_axis.patch.set_facecolor('none')
    for name, spine in Q2_axis.spines.items():
        spine.set_visible(name == 'bottom')
        spine.set_linewidth(1.0)
        spine.set_color('black')

    # nu axis setup
    nu_axis.set_ylim(nu_edges[0], nu_edges[-1])
    nu_axis.set_yticks(nu_edges)
    nu_axis.set_yticklabels([f"{int(a)}" for a in nu_edges], fontsize=12)
    nu_axis.set_xticks([])
    nu_axis.yaxis.tick_left()
    nu_axis.set_ylabel(r"$\nu$ [GeV]", fontsize=14, labelpad=20)
    nu_axis.patch.set_facecolor('none')
    for name, spine in nu_axis.spines.items():
        spine.set_visible(name == 'left')
        spine.set_linewidth(1.0)
        spine.set_color('black')

    plt.tight_layout(rect=[0.10, 0.10, 0.9, 0.95])
    plt.savefig("acceptance_integrated_t.png", dpi=300)

#plot_acceptance_integrated_t()

# **********************************
# Overwrite arrays to store the acceptance/weighted sum values for each nu bin (integrate over phi)
# *** mu+ ***
rec_weights_muPlus = np.zeros((4,9,11)) 
gen_weights_muPlus = np.zeros((4,9,11)) 
acceptance_muPlus  = np.zeros((4,9,11)) 

acceptance_err_muPlus = np.zeros((4,9,11))
rec_weight2_muPlus = np.zeros((4,9,11))
gen_weight2_muPlus = np.zeros((4,9,11)) 

# *** mu- ***
rec_weights_muMinus = np.zeros((4,9,11)) 
gen_weights_muMinus = np.zeros((4,9,11)) 
acceptance_muMinus  = np.zeros((4,9,11)) 

acceptance_err_muMinus = np.zeros((4,9,11)) 
rec_weight2_muMinus = np.zeros((4,9,11))  
gen_weight2_muMinus = np.zeros((4,9,11)) 

# **********************************
# Function for filling the weight arrays (integrating over t)
def fill_weights_int_phi(tree, elist, charge, data_type):
    tree.SetEntryList(elist)
    for event in tree:
        Q2 = event.Q2
        nu = event.nu
        t = np.abs(event.t)
        weight_BH = event.weight_BH

        # Determine bins
        nu_bin = int(math.floor((nu - 10.) / 2.))
        Q2_bin = int(math.floor((Q2 - 1.) / 1.))
        t_bin = 77
        for i, (t_low, t_high) in enumerate(t_bins):
            if t_low < t < t_high:
                t_bin = i
                break

        if (t_bin == 77) or (Q2_bin >= len(Q2_bins)) or (nu_bin >= len(nu_bins)):
            continue

        if data_type == "rec":
            if charge == "muPlus":
                rec_weights_muPlus[t_bin][Q2_bin][nu_bin] += weight_BH
                rec_weight2_muPlus[t_bin][Q2_bin][nu_bin] += weight_BH ** 2
            elif charge == "muMinus":
                rec_weights_muMinus[t_bin][Q2_bin][nu_bin] += weight_BH
                rec_weight2_muMinus[t_bin][Q2_bin][nu_bin] += weight_BH ** 2
        elif data_type == "gen":
            if charge == "muPlus":
                gen_weights_muPlus[t_bin][Q2_bin][nu_bin] += weight_BH
                gen_weight2_muPlus[t_bin][Q2_bin][nu_bin] += weight_BH ** 2
            elif charge == "muMinus":
                gen_weights_muMinus[t_bin][Q2_bin][nu_bin] += weight_BH
                gen_weight2_muMinus[t_bin][Q2_bin][nu_bin] += weight_BH ** 2

# Function for finding the statistical errors (integrating over t)
def acc_error_int_phi(rec_weights, gen_weights, rec_var, gen_var, charge):
    for i in range(4):
        for j in range(9):
            for k in range(11):
                weightRec = rec_weights[i][j][k]
                weightGen = gen_weights[i][j][k]
                varRec = rec_var[i][j][k]
                varGen = gen_var[i][j][k]

                if weightGen == 0:
                    errAcc = 0.0
                else:
                    errT1 = (1 / weightGen)**2 * varRec
                    errT2 = (weightRec / (weightGen ** 2))**2 * varGen
                    errAcc = np.sqrt(errT1 + errT2)

                if charge == "muPlus":
                    acceptance_err_muPlus[i][j][k] = errAcc
                elif charge == "muMinus":
                    acceptance_err_muMinus[i][j][k] = errAcc

# **********************************
# Fill the weight arrays and calculate the acceptance for mu+ and mu- 
fill_weights_int_phi(tree_rec, elist_muPlus_rec, "muPlus", "rec")
fill_weights_int_phi(tree_rec, elist_muMinus_rec, "muMinus", "rec") 
fill_weights_int_phi(tree_gen, elist_muPlus_gen, "muPlus", "gen")
fill_weights_int_phi(tree_gen, elist_muMinus_gen, "muMinus", "gen")

# *** mu+ events *** 
nonzero_mask_muPlus = gen_weights_muPlus != 0 # protect against division by zero 
acceptance_muPlus[nonzero_mask_muPlus] = rec_weights_muPlus[nonzero_mask_muPlus] / gen_weights_muPlus[nonzero_mask_muPlus]
acc_error_int_phi(rec_weights_muPlus, gen_weights_muPlus, rec_weight2_muPlus, gen_weight2_muPlus, "muPlus")

# *** mu- events ***
nonzero_mask_muMinus = gen_weights_muMinus != 0 # protect against division by zero 
acceptance_muMinus[nonzero_mask_muMinus] = rec_weights_muMinus[nonzero_mask_muMinus] / gen_weights_muMinus[nonzero_mask_muMinus]
acc_error_int_phi(rec_weights_muMinus, gen_weights_muMinus, rec_weight2_muMinus, gen_weight2_muMinus, "muMinus")

# **********************************
# Plot the dsitrubtion of the events by bin (integrated over phi)
def plot_acceptance_integrated_phi():
    nu_bin_centers = 0.5 * (nu_edges[:-1] + nu_edges[1:])
    n_t_bins, n_Q2_bins, n_nu_bins = acceptance_muPlus.shape

    fig, axes = plt.subplots(nrows=n_Q2_bins, ncols=n_t_bins, figsize=(18, 22), sharex=True, sharey=True)

    for i_Q2 in range(n_Q2_bins):      # Q2 bins (y-axis, rows)
        for j_t in range(n_t_bins):    # t bins (x-axis, cols)
            ax = axes[i_Q2, j_t]

            # acceptance_muPlus and acceptance_muMinus shape: (t, Q2, nu)
            y_muPlus = acceptance_muPlus[i, j]
            yerr_muPlus = acceptance_err_muPlus[i, j]
            y_muMinus = acceptance_muMinus[i, j]
            yerr_muMinus = acceptance_err_muMinus[i, j]

            ax.errorbar(nu_bin_centers, y_muPlus, fmt='o', color='red', markersize=3, label='μ⁺' if i_Q2 == 0 and j_t == 0 else "")
            ax.errorbar(nu_bin_centers, y_muMinus, fmt='o', color='black', markersize=3, label='μ⁻' if i_Q2 == 0 and j_t == 0 else "")

            ax.set_ylim(0, 0.8)
            ax.grid(True, linestyle='--', linewidth=0.5)
            ax.axhline(0, color='gray', linewidth=0.5)

    # Legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize=14, markerscale=1.5)

    # Create a new set of axes for the phi and acceptance scale at the top right corner
    nu_axis = fig.add_axes([0.71, 0.92, 0.185, 0.035])  # [left, bottom, width, height]
    acc_axis = fig.add_axes([0.8, 0.86, 0.11, 0.0825])

    # Setup phi axis
    nu_axis.set_xlim(10,32)
    nu_axis.set_xticks(nu_edges)
    nu_axis.set_xticklabels([f"{edge:.0f}" for edge in nu_edges], fontsize=12)
    nu_axis.set_yticks([])
    nu_axis.yaxis.set_visible(False)
    nu_axis.tick_params(axis='x', direction='in', length=5, top=True, bottom=False)
    nu_axis.xaxis.set_label_position('top')
    nu_axis.set_xlabel(r"$\nu$ [GeV]", fontsize=14, labelpad=20)
    nu_axis.xaxis.tick_top()
    nu_axis.patch.set_facecolor('none')
    for name, spine in nu_axis.spines.items():
        spine.set_visible(name == 'top')
        if name == 'top':
            spine.set_linewidth(1.0)
            spine.set_color('black')

    # Setup acceptance axis
    acc_axis.set_ylim(0, 0.8)
    acc_axis.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    acc_axis.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8"], fontsize=12)
    acc_axis.set_xticks([])
    acc_axis.xaxis.set_visible(False)
    acc_axis.set_ylabel("Acceptance", fontsize=14, labelpad=15)
    acc_axis.yaxis.set_label_position('right')
    acc_axis.yaxis.tick_right()
    acc_axis.tick_params(axis='y', direction='in', length=5)
    acc_axis.patch.set_facecolor('none')
    for name, spine in acc_axis.spines.items():
        spine.set_visible(name == 'right')
        if name == 'right':
            spine.set_linewidth(1.0)
            spine.set_color('black')

    # Global axes for t and Q2 
    t_axis = fig.add_axes([0.125, 0.09, 0.775, 0.035])
    Q2_axis = fig.add_axes([0.083, 0.11, 0.034, 0.835])

    # t axis setup
    tick_positions = np.arange(len(t_edges))  
    t_axis.set_xlim(0, len(t_edges) - 1)  # Ensure full width of bins is shown
    t_axis.set_xticks(tick_positions)
    t_axis.set_xticklabels([f"{a}" for a in t_edges], fontsize=12)
    t_axis.set_yticks([])
    t_axis.xaxis.tick_bottom()
    t_axis.set_xlabel(r"$|t|$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
    t_axis.patch.set_facecolor('none')
    for name, spine in t_axis.spines.items():
        spine.set_visible(name == 'bottom')
        spine.set_linewidth(1.0)
        spine.set_color('black')

    # Q2 axis setup
    Q2_axis.set_ylim(Q2_edges[0], Q2_edges[-1])
    Q2_axis.set_yticks(Q2_edges)
    Q2_axis.set_yticklabels([f"{int(a)}" for a in Q2_edges], fontsize=12)
    Q2_axis.set_xticks([])
    Q2_axis.yaxis.tick_left()
    Q2_axis.set_ylabel(r"$Q^2$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
    Q2_axis.patch.set_facecolor('none')
    for name, spine in Q2_axis.spines.items():
        spine.set_visible(name == 'left')
        spine.set_linewidth(1.0)
        spine.set_color('black')

    plt.tight_layout(rect=[0.10, 0.10, 0.9, 0.95])
    plt.savefig("acceptance_integrated_phi.png", dpi=300)

#plot_acceptance_integrated_phi()

# **********************************
# Reset tree access
tree_rec.SetEntryList(0)
tree_gen.SetEntryList(0)