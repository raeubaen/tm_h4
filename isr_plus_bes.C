TGraph* isr;

double func(Double_t *x, Double_t *par) { return isr->Eval(x[0])*TMath::Gaus(x[0], 212, par[0])*1/(sqrt(2*3.14*par[0]*par[0]));}

void isr_plus_bes(){
  auto *ff = new TFile("isr.root");
  isr = (TGraph*)ff->Get("isr");


  TF1 *f = new TF1("func", func, 212-1.4e-3, 217, 1);
  f->SetParameters(1);
  cout << "XS with 1 MeV gaus bes + ISR; " << f->Integral(212-1.4e-3, 217) << " pb" << endl;
  f->SetParameters(0.3);
  cout << "XS with 0.3 MeV gaus bes + ISR; " << f->Integral(212-1.4e-3, 217) << " pb" << endl;
}
