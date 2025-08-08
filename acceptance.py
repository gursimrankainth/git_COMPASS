import ROOT
from ROOT import TLorentzVector
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

#TODO: Plot the distributions for Q2, nu, x and t for the MC data and real data
#TODO: Filter the reconstructed and generated data so each event is only processed once 

# *** !WARNING! ***
# Do not move the plotting lines to different locations! The arrays are overwritten before each
# set of plots is made to use less memory.  

# **********************************
# Reconstructed MC data after full DVCS selection
file_rec = ROOT.TFile.Open("/afs/cern.ch/user/g/gkainth/merged_P09_test.root")
tree_rec = file_rec.Get("USR970")

# Generated MC data 
file_gen = ROOT.TFile.Open("/afs/cern.ch/user/g/gkainth/merged_P09_gen_test.root")
tree_gen = file_gen.Get("USR970_gen")

print("Total Entries: Reconstructed:", tree_rec.GetEntries(), ", Generated:", tree_gen.GetEntries())

# **********************************
# Define the bin edges/bins for wide acceptance calculation 
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
# Function for drawing the bin lines 
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

# Function for making default plots 
def default_plots():
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

  # Create canvases to draw and save the histograms
  c1 = ROOT.TCanvas("c1", "Q2 vs nu", 800, 600)
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
  c3.SaveAs("phiRV.png")

#default_plots()

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
# Define a helper function to overwrite arrays 
def reset_arrays(shape=(4,4,8)):
  global rec_weights_muPlus, gen_weights_muPlus, acceptance_muPlus, acceptance_err_muPlus
  global rec_weight2_muPlus, gen_weight2_muPlus
  global rec_weights_muMinus, gen_weights_muMinus, acceptance_muMinus, acceptance_err_muMinus
  global rec_weight2_muMinus, gen_weight2_muMinus
  
  rec_weights_muPlus  = np.zeros(shape)
  gen_weights_muPlus  = np.zeros(shape)
  acceptance_muPlus   = np.zeros(shape)
  acceptance_err_muPlus = np.zeros(shape)
  rec_weight2_muPlus  = np.zeros(shape)
  gen_weight2_muPlus  = np.zeros(shape)

  rec_weights_muMinus = np.zeros(shape)
  gen_weights_muMinus = np.zeros(shape)
  acceptance_muMinus  = np.zeros(shape)
  acceptance_err_muMinus = np.zeros(shape)
  rec_weight2_muMinus = np.zeros(shape)
  gen_weight2_muMinus = np.zeros(shape)

# **********************************
# Define arrays to store the acceptance/weighted sum values for each phi bin (integrate over t)
# *** mu+ ***
rec_weights_muPlus = np.zeros((9,11,8)) # array that holds the weighted sums for reconstructed data 
gen_weights_muPlus = np.zeros((9,11,8)) # array that holds the weighted sums for generated data 
acceptance_muPlus  = np.zeros((9,11,8)) # array that stores the acceptance values 
acceptance_err_muPlus = np.zeros((9,11,8)) # array that stores the estimated statistical errors 
rec_weight2_muPlus = np.zeros((9,11,8))  # holds weight squared for each bin (reconstructed) 
gen_weight2_muPlus = np.zeros((9,11,8))  # holds weight squared for each bin (generated)

# *** mu- ***
rec_weights_muMinus = np.zeros((9,11,8)) # array that holds the weighted sums for reconstructed data 
gen_weights_muMinus = np.zeros((9,11,8)) # array that holds the weighted sums for generated data 
acceptance_muMinus  = np.zeros((9,11,8)) # array that stores the acceptance values 
acceptance_err_muMinus = np.zeros((9,11,8)) # array that stores the estimated statistical errors
rec_weight2_muMinus = np.zeros((9,11,8))  # holds weight squared for each bin (reconstructed) 
gen_weight2_muMinus = np.zeros((9,11,8))  # holds weight squared for each bin (generated)

# **********************************
# Function for filling the weight arrays
def fill_weights(tree, elist, charge, data_type, Q2_bins=Q2_bins, nu_bins=nu_bins, 
                t_bins=t_bins, phi_bins=phi_bins, integrate="t"):

  for i in range(elist.GetN()):
      entry_index = elist.GetEntry(i)
      tree.GetEntry(entry_index)
      event = tree # alias -> the line aboves pulls the exact event from the tree 

      Q2 = event.Q2
      nu = event.nu
      t = np.abs(event.t)
      phi = event.phi_gg
      weight_BH = event.weight_BH

      # Bin Q2
      Q2_bin = int(math.floor((Q2 - 1.) / 1.))
      if Q2_bin < 0 or Q2_bin >= len(Q2_bins):
          continue

      # Bin nu
      nu_bin = int(math.floor((nu - 10.) / 2.))
      if nu_bin < 0 or nu_bin >= len(nu_bins):
          continue

      if integrate == "t":
          # Bin phi
          phi_bin = int(math.floor((phi + np.pi) / (np.pi/4)))
          if phi_bin < 0 or phi_bin >= len(phi_bins):
              continue
          # Axis order: [Q2][nu][phi]
          i1, i2, i3 = Q2_bin, nu_bin, phi_bin

      elif integrate == "phi":
          # Bin |t|
          t_bin = -1
          for i, (t_low, t_high) in enumerate(t_bins):
              if t_low < t < t_high:
                  t_bin = i
                  break
          if t_bin < 0 or t_bin >= len(t_bins):
              continue
          # Axis order: [t][Q2][nu]
          i1, i2, i3 = t_bin, Q2_bin, nu_bin

      else:
          raise ValueError("integrate must be 't' or 'phi'")

      # Fill the appropriate arrays
      if data_type == "rec":
          if charge == "muPlus":
              rec_weights_muPlus[i1][i2][i3] += weight_BH
              rec_weight2_muPlus[i1][i2][i3] += weight_BH ** 2
          elif charge == "muMinus":
              rec_weights_muMinus[i1][i2][i3] += weight_BH
              rec_weight2_muMinus[i1][i2][i3] += weight_BH ** 2
      elif data_type == "gen":
          if charge == "muPlus":
              gen_weights_muPlus[i1][i2][i3] += weight_BH
              gen_weight2_muPlus[i1][i2][i3] += weight_BH ** 2
          elif charge == "muMinus":
              gen_weights_muMinus[i1][i2][i3] += weight_BH
              gen_weight2_muMinus[i1][i2][i3] += weight_BH ** 2

# **********************************
# Function for finding the statistical errors
def acc_error(rec_weights, gen_weights, rec_weights2, gen_weights2, charge, Q2_bins=Q2_bins, 
              nu_bins=nu_bins, t_bins=t_bins, phi_bins=phi_bins, integrate="t"):
  if integrate == "t":
      # Binning shape: (9, 11, 8) for Q2, nu, phi
      shape = (len(Q2_bins), len(nu_bins), len(phi_bins))
  elif integrate == "phi":
      # Binning shape: (4, 9, 11) for t, Q2, nu
      shape = (len(t_bins), len(Q2_bins), len(nu_bins))
  else:
      raise ValueError("integrate must be 't' or 'phi'")

  for idx in np.ndindex(shape):
      weightRec = rec_weights[idx]
      weightGen = gen_weights[idx]
      weight2Rec = rec_weights2[idx]
      weight2Gen = gen_weights2[idx]

      if weightGen == 0:
          errAcc = 0.0
      else:
          errT1 = ((1 / weightGen))**2 * weight2Rec
          errT2 = ((weightRec / weightGen))**2 * weight2Gen
          errAcc = np.sqrt(errT1 + errT2)

      if charge == "muPlus":
          acceptance_err_muPlus[idx] = errAcc
      elif charge == "muMinus":
          acceptance_err_muMinus[idx] = errAcc

# **********************************
# Fill the weight arrays and calculate the acceptance for mu+ and mu- 
fill_weights(tree_rec, elist_muPlus_rec, "muPlus", "rec", integrate="t")
fill_weights(tree_gen, elist_muPlus_gen, "muPlus", "gen", integrate="t")
fill_weights(tree_rec, elist_muMinus_rec, "muMinus", "rec", integrate="t") 
fill_weights(tree_gen, elist_muMinus_gen, "muMinus", "gen", integrate="t")

# *** mu+ events *** 
nonzero_mask_muPlus = gen_weights_muPlus != 0 # protect against division by zero 
acceptance_muPlus[nonzero_mask_muPlus] = rec_weights_muPlus[nonzero_mask_muPlus] / gen_weights_muPlus[nonzero_mask_muPlus]
acc_error(rec_weights_muPlus, gen_weights_muPlus, rec_weight2_muPlus, gen_weight2_muPlus, "muPlus", integrate="t")

# *** mu- events ***
nonzero_mask_muMinus = gen_weights_muMinus != 0 # protect against division by zero 
acceptance_muMinus[nonzero_mask_muMinus] = rec_weights_muMinus[nonzero_mask_muMinus] / gen_weights_muMinus[nonzero_mask_muMinus]
acc_error(rec_weights_muMinus, gen_weights_muMinus, rec_weight2_muMinus, gen_weight2_muMinus, "muMinus", integrate="t")

# **********************************
# Plot the dsitrubtion of the events by bin (integrated over t)
def plot_acceptance_integrated_t():
  phi_bin_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
  n_Q2_bins, n_nu_bins, n_phi_bins = acceptance_muPlus.shape

  fig, axes = plt.subplots(nrows=n_nu_bins, ncols=n_Q2_bins, figsize=(18, 22), sharex=True, sharey=True)

  for i in range(n_nu_bins):      # nu bins (y-axis)
      for j in range(n_Q2_bins):  # Q2 bins (x-axis)
          ax = axes[i, j]

          y_muPlus = acceptance_muPlus[j, i]
          yerr_muPlus = acceptance_err_muPlus[j, i]
          y_muMinus = acceptance_muMinus[j, i]
          yerr_muMinus = acceptance_err_muMinus[j, i]

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
reset_arrays(shape=(4,9,11))

# Fill the weight arrays and calculate the acceptance for mu+ and mu- 
fill_weights(tree_rec, elist_muPlus_rec, "muPlus", "rec", integrate="phi")
fill_weights(tree_gen, elist_muPlus_gen, "muPlus", "gen", integrate="phi")
fill_weights(tree_rec, elist_muMinus_rec, "muMinus", "rec", integrate="phi") 
fill_weights(tree_gen, elist_muMinus_gen, "muMinus", "gen", integrate="phi")

# *** mu+ events *** 
nonzero_mask_muPlus = gen_weights_muPlus != 0 # protect against division by zero 
acceptance_muPlus[nonzero_mask_muPlus] = rec_weights_muPlus[nonzero_mask_muPlus] / gen_weights_muPlus[nonzero_mask_muPlus]
acc_error(rec_weights_muPlus, gen_weights_muPlus, rec_weight2_muPlus, gen_weight2_muPlus, "muPlus", integrate="phi")

# *** mu- events ***
nonzero_mask_muMinus = gen_weights_muMinus != 0 # protect against division by zero 
acceptance_muMinus[nonzero_mask_muMinus] = rec_weights_muMinus[nonzero_mask_muMinus] / gen_weights_muMinus[nonzero_mask_muMinus]
acc_error(rec_weights_muMinus, gen_weights_muMinus, rec_weight2_muMinus, gen_weight2_muMinus, "muMinus", integrate="phi")

# **********************************
# Plot the dsitrubtion of the events by bin (integrated over phi)
def plot_acceptance_integrated_phi():
  nu_bin_centers = 0.5 * (nu_edges[:-1] + nu_edges[1:])
  n_t_bins, n_Q2_bins, n_nu_bins = acceptance_muPlus.shape

  fig, axes = plt.subplots(nrows=n_Q2_bins, ncols=n_t_bins, figsize=(18, 22), sharex=True, sharey=True)

  for i in range(n_Q2_bins): # rows (Q2 -> y values)    
      for j in range(n_t_bins): # columns (t -> x)   
          ax = axes[i, j]

          # acceptance_muPlus and acceptance_muMinus shape: (t, Q2, nu)
          y_muPlus = acceptance_muPlus[j, i]
          yerr_muPlus = acceptance_err_muPlus[j, i]
          y_muMinus = acceptance_muMinus[j, i]
          yerr_muMinus = acceptance_err_muMinus[j, i]

          ax.errorbar(nu_bin_centers, y_muPlus, yerr=yerr_muPlus, fmt='o', color='red', markersize=3, label='μ⁺' if i == 0 and j == 0 else "")
          ax.errorbar(nu_bin_centers, y_muMinus, yerr=yerr_muMinus, fmt='o', color='black', markersize=3, label='μ⁻' if i == 0 and j == 0 else "")

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
# Tighten bins for final 4D acceptance plots
# nu: 4 bins of width 5.5 GeV between 10 and 32 GeV
nu_edges_tight = np.linspace(10, 32, 5)
nu_bins_tight = list(zip(nu_edges_tight[:-1], nu_edges_tight[1:]))
# Q2: 4 bins of width 1 (GeV/c)^2 between 1 and 5
Q2_edges_tight = np.linspace(1, 5, 5)
Q2_bins_tight = list(zip(Q2_edges_tight[:-1], Q2_edges_tight[1:]))
# Binning for |t| and phi remains the same

# **********************************
# Define lists to store the acceptance/weighted sum arrays for each |t| bin 
# *** mu+ *** 
rec_weights_muPlus_tBins = [] 
gen_weights_muPlus_tBins = [] 
acceptance_muPlus_tBins = []
acceptance_err_muPlus_tBins = [] 
rec_weight2_muPlus_tBins = []
gen_weight2_muPlus_tBins = []

# *** mu- ***
rec_weights_muMinus_tBins = []
gen_weights_muMinus_tBins = []
acceptance_muMinus_tBins  = []
acceptance_err_muMinus_tBins = []
rec_weight2_muMinus_tBins = []
gen_weight2_muMinus_tBins = []

# **********************************
# Overwrite the arrays with the new shape (will have one copy of these for each bin)
reset_arrays((4,4,8))

# Bin the data in t
def bin_in_t(tree_rec, tree_gen):
  elist_muPlus_rec_by_t = []
  elist_muMinus_rec_by_t = []
  elist_muPlus_gen_by_t = []
  elist_muMinus_gen_by_t = []

  for t_min, t_max in t_bins:
      cut_muPlus_rec = f"Q_beam > 0 && abs(t) >= {t_min} && abs(t) < {t_max}"
      cut_muMinus_rec = f"Q_beam < 0 && abs(t) >= {t_min} && abs(t) < {t_max}"
      cut_muPlus_gen = f"Q_beam > 0 && abs(t) >= {t_min} && abs(t) < {t_max}"
      cut_muMinus_gen = f"Q_beam < 0 && abs(t) >= {t_min} && abs(t) < {t_max}"

      tree_rec.Draw(f">>elist_muPlus_rec_{t_min}_{t_max}", cut_muPlus_rec, "entrylist")
      elist_muPlus_rec_by_t.append(ROOT.gDirectory.Get(f"elist_muPlus_rec_{t_min}_{t_max}"))

      tree_rec.Draw(f">>elist_muMinus_rec_{t_min}_{t_max}", cut_muMinus_rec, "entrylist")
      elist_muMinus_rec_by_t.append(ROOT.gDirectory.Get(f"elist_muMinus_rec_{t_min}_{t_max}"))

      tree_gen.Draw(f">>elist_muPlus_gen_{t_min}_{t_max}", cut_muPlus_gen, "entrylist")
      elist_muPlus_gen_by_t.append(ROOT.gDirectory.Get(f"elist_muPlus_gen_{t_min}_{t_max}"))

      tree_gen.Draw(f">>elist_muMinus_gen_{t_min}_{t_max}", cut_muMinus_gen, "entrylist")
      elist_muMinus_gen_by_t.append(ROOT.gDirectory.Get(f"elist_muMinus_gen_{t_min}_{t_max}"))

  return elist_muPlus_rec_by_t, elist_muMinus_rec_by_t, elist_muPlus_gen_by_t, elist_muMinus_gen_by_t

elist_muPlus_rec_by_t, elist_muMinus_rec_by_t, elist_muPlus_gen_by_t, elist_muMinus_gen_by_t = bin_in_t(tree_rec, tree_gen)

# Fill the weights for each t bin 
def fill_weights_4D(tree, elist, charge, data_type):
  for elist_bin in elist:
      reset_arrays((4,4,8)) # reset the arrays for each bin 
      fill_weights(tree, elist_bin, charge, data_type, Q2_bins=Q2_bins_tight, nu_bins=nu_bins_tight, integrate="t")

      # Save copies of the filled arrays to the correct lists 
      if data_type == "rec":
          if charge == "muPlus":
              rec_weights_muPlus_tBins.append(rec_weights_muPlus.copy()) 
              rec_weight2_muPlus_tBins.append(rec_weight2_muPlus.copy())
          elif charge == "muMinus":
              rec_weights_muMinus_tBins.append(rec_weights_muMinus.copy())
              rec_weight2_muMinus_tBins.append(rec_weight2_muMinus.copy())
      elif data_type == "gen":
          if charge == "muPlus":
              gen_weights_muPlus_tBins.append(gen_weights_muPlus.copy()) 
              gen_weight2_muPlus_tBins.append(gen_weight2_muPlus.copy())
          elif charge == "muMinus":
              gen_weights_muMinus_tBins.append(gen_weights_muMinus.copy())
              gen_weight2_muMinus_tBins.append(gen_weight2_muMinus.copy())

fill_weights_4D(tree_rec, elist_muPlus_rec_by_t, "muPlus", "rec")
fill_weights_4D(tree_gen, elist_muPlus_gen_by_t, "muPlus", "gen")
fill_weights_4D(tree_rec, elist_muMinus_rec_by_t, "muMinus", "rec")
fill_weights_4D(tree_gen, elist_muMinus_gen_by_t, "muMinus", "gen")

# Calculate the acceptance and error for each bin 
for bin_idx in range(len(t_bins)):
  # *** mu+ events *** 
  acceptance_muPlus = np.zeros((4,4,8))
  nonzero_mask_muPlus = gen_weights_muPlus_tBins[bin_idx] != 0
  acceptance_muPlus[nonzero_mask_muPlus] = (
      rec_weights_muPlus_tBins[bin_idx][nonzero_mask_muPlus] / 
      gen_weights_muPlus_tBins[bin_idx][nonzero_mask_muPlus]
  )
  acceptance_muPlus_tBins.append(acceptance_muPlus.copy())

  acceptance_err_muPlus.fill(0) # reset the array before filling 
  acc_error(
      rec_weights_muPlus_tBins[bin_idx], gen_weights_muPlus_tBins[bin_idx], 
      rec_weight2_muPlus_tBins[bin_idx], gen_weight2_muPlus_tBins[bin_idx], 
      "muPlus", Q2_bins=Q2_bins_tight, nu_bins=nu_bins_tight, integrate="t"
  )
  acceptance_err_muPlus_tBins.append(acceptance_err_muPlus.copy())  # copy global array

  # *** mu- events ***
  acceptance_muMinus = np.zeros((4,4,8))
  nonzero_mask_muMinus = gen_weights_muMinus_tBins[bin_idx] != 0
  acceptance_muMinus[nonzero_mask_muMinus] = (
      rec_weights_muMinus_tBins[bin_idx][nonzero_mask_muMinus] / 
      gen_weights_muMinus_tBins[bin_idx][nonzero_mask_muMinus]
  )
  acceptance_muMinus_tBins.append(acceptance_muMinus.copy())

  acceptance_err_muMinus.fill(0) # reset the array before filling it 
  acc_error(
      rec_weights_muMinus_tBins[bin_idx], gen_weights_muMinus_tBins[bin_idx], 
      rec_weight2_muMinus_tBins[bin_idx], gen_weight2_muMinus_tBins[bin_idx], 
      "muMinus", Q2_bins=Q2_bins_tight, nu_bins=nu_bins_tight, integrate="t"
  )
  acceptance_err_muMinus_tBins.append(acceptance_err_muMinus.copy())  # copy global array

# **********************************
# Plot the dsitrubtion of the events by |t| bin (integrated over t)
def plot_acceptance_by_tBin(muPlus_tBins=acceptance_muPlus_tBins, err_muPlus_tBins=acceptance_err_muPlus_tBins, 
                            muMinus_tBins=acceptance_muMinus_tBins, err_muMinus_tBins=acceptance_err_muMinus_tBins,
                            idx=0):
  phi_bin_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
  n_Q2_bins, n_nu_bins, n_phi_bins = muPlus_tBins[idx].shape

  fig, axes = plt.subplots(nrows=n_nu_bins, ncols=n_Q2_bins, figsize=(22, 22), sharex=True, sharey=True)

  for i in range(n_nu_bins):      # nu bins (y-axis)
    for j in range(n_Q2_bins):  # Q2 bins (x-axis)
      ax = axes[i, j]

      y_muPlus = muPlus_tBins[idx][j, i, :] 
      yerr_muPlus = err_muPlus_tBins[idx][j, i, :]
      y_muMinus = muMinus_tBins[idx][j, i, :] 
      yerr_muMinus = err_muMinus_tBins[idx][j, i, :]

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
  phi_axis = fig.add_axes([0.71, 0.91, 0.18, 0.05])  # [left, bottom, width, height]
  acc_axis = fig.add_axes([0.86, 0.745, 0.05, 0.195])

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
  Q2_axis.set_xlim(Q2_edges_tight[0], Q2_edges_tight[-1])
  Q2_axis.set_xticks(Q2_edges_tight)
  Q2_axis.set_xticklabels([f"{int(a)}" for a in Q2_edges_tight], fontsize=12)
  Q2_axis.set_yticks([])
  Q2_axis.xaxis.tick_bottom()
  Q2_axis.set_xlabel(r"$Q^2$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
  Q2_axis.patch.set_facecolor('none')
  for name, spine in Q2_axis.spines.items():
    spine.set_visible(name == 'bottom')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  # nu axis setup
  nu_axis.set_ylim(nu_edges_tight[0], nu_edges_tight[-1])
  nu_axis.set_yticks(nu_edges_tight)
  nu_axis.set_yticklabels([f"{int(a)}" for a in nu_edges_tight], fontsize=12)
  nu_axis.set_xticks([])
  nu_axis.yaxis.tick_left()
  nu_axis.set_ylabel(r"$\nu$ [GeV]", fontsize=14, labelpad=20)
  nu_axis.patch.set_facecolor('none')
  for name, spine in nu_axis.spines.items():
    spine.set_visible(name == 'left')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  plt.tight_layout(rect=[0.10, 0.10, 0.9, 0.95])
  plt.savefig(f"acceptance_integrated_bin{idx+1}.png", dpi=300)

for i in range(len(t_bins)):
  plot_acceptance_by_tBin(idx=i)