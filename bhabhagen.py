import ROOT
import sys
import numpy as np
from tqdm import tqdm
import uproot

def diff_bhabha(t):
  #dsigma/dt
  return ((1 + np.cos(t/2)**4)/(np.sin(t/2)**4) - 2*(np.cos(t/2)**4)/(np.sin(t/2)**2) + (1 + np.cos(t)**2)/2)*np.sin(t)

if __name__ == "__main__":
  nev = int(float(sys.argv[1]))
  #egamma_h = ROOT.TFile("egamma_signal_h.root").Get("h")
  e_comp = np.zeros((nev, 3))
  e_mom = np.zeros((nev,))
  p_comp = np.zeros((nev, 3))
  p_mom = np.zeros((nev,))
  w = np.zeros((nev,))
  e_theta_cm = np.zeros((nev,))
  p_theta_cm = np.zeros((nev,))
  for i in tqdm(range(nev)):
    tm_vec = ROOT.TLorentzVector(0, 0, 43.7, 0.511e-3+np.sqrt(43.7*43.7 + 0.511e-3*0.511e-3))

    tm_ps = ROOT.TGenPhaseSpace()
    tm_ps.SetDecay(tm_vec, 2, np.array([0.511e-3, 0.511e-3]))

    _w = tm_ps.Generate()
    e = tm_ps.GetDecay(0)
    p = tm_ps.GetDecay(1)

    e_comp[i, :] = e.Vect()
    e_mom[i] = e.P()
    p_comp[i, :] = p.Vect()
    p_mom[i] = p.P()

    e_copy = ROOT.TLorentzVector(e)
    e_copy.Boost(0, 0, -tm_vec.Beta())
    _e_theta_cm = e_copy.Theta()
    p_copy = ROOT.TLorentzVector(p)
    p_copy.Boost(0, 0, -tm_vec.Beta())
    _p_theta_cm = p_copy.Theta()

    _w *= diff_bhabha(_e_theta_cm)
    e_theta_cm[i] = _e_theta_cm
    p_theta_cm[i] = _p_theta_cm
    w[i] = _w

  f = uproot.recreate("bhabhagen.root")
  f["events"] = {"e_mom": e_mom, "e_comp": e_comp, "p_mom": p_mom, "p_comp": p_comp, "w": w, "e_theta_cm": e_theta_cm, "p_theta_cm": p_theta_cm}
  f.close()

  '''
  presel_cuts = (theta_ep > 177/180.*np.pi)

  final_cuts = np.logical_and(dr/1000 > 5, (vl/1000)**2 + (dr/1000-4)**2 > 6**2)
  other_cuts = np.logical_and(dr/1000 < 44, vl/1000 < 40)
  indexes = np.logical_and(presel_cuts, final_cuts)
  #indexes = np.logical_and(indexes, other_cuts)
  sumw_aftercuts = (w[indexes]).sum()
  eff = sumw_aftercuts/sumw_nocuts
  print(eff)
  '''

