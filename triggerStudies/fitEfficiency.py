import ROOT as r

f 		 = r.TFile.Open("trig_eff_2016.root")
eff 	 = f.Get("trig_eff")

f_cauchy = r.TF1("f_cauchy", "TMath::ATan((x - [0]) / [1]) / TMath::Pi() + 0.5", 450., 1000.)
f_cauchy.SetParNames("median", "gamma");
f_cauchy.SetParameters(500., 50.);

f_gaus = r.TF1("f_gaus", "0.5 * (1.0 + TMath::Erf((x - [0]) / [1] / TMath::Sqrt2()))", 450., 1000.);
f_gaus.SetParNames("mean", "sigma");
f_gaus.SetParameters(500., 50.);

fitResultPtr = eff.Fit(f_cauchy);
#fitResultPtr = eff.Fit(f_gaus);
eff.Draw("AP")

name = input("End? ")