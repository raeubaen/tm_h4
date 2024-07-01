double beta(double ecm){
  return 2/(3.14*137)*(2*log(ecm/0.511) - 1);
}

double f_coeff(double ecm){
  double b = beta(ecm);
  return exp(b/4 + (0.5+3.14*3.14/3)/(3.14*137) - 0.577*b)/TMath::Gamma(1+b)*b;
}

double piece2_int(double x, double b){
  return (x-1)*pow(1-x, b)*(b*x+b+x+3)/((b+1)*(b+2));
}

double f_integral(double x1, double x2, double ecm){
  double b = beta(ecm);
  double piece1 = pow(1-x1, b)/b - pow(1-x2, b)/b;
  double piece2 = piece2_int(x2, b) - piece2_int(x1, b);
  return f_coeff(ecm)*((1+b/2)*piece1 - 0.5*piece2);
}

double sigma_isr(double ecm){
  double x_min = min(1., (212-1.4e-3)/ecm), x_max = min(1., 212/ecm);
  return 66600*f_integral(x_min, x_max, ecm);
}

void isr(){
  int n=100000;
  double *ecm_list = new double[n]();
  double *sigma_ecm = new double[n]();
  for(int i=0; i<n; i++){
    ecm_list[i] = 212 - 0.0015 + (5./n)*i;
    sigma_ecm[i] = sigma_isr(ecm_list[i]);
    cout << "dsigma/decm " << ecm_list[i] << "," << sigma_ecm[i] << endl;
  }
  auto *c = new TCanvas("c");
  TGraph *g = new TGraph(n, ecm_list, sigma_ecm);
  g->SetName("isr");
  g->SetTitle("isr");
  g->Draw("AP");
  g->SaveAs("isr.root");
}
