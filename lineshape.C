
double n_k(double *x, double *par){
	float n0 = par[0];
	float alfa_plus = par[1];
	float alfa_minus = par[2];
	float beta_plus = par[3];
	float beta_minus = par[4];
	float kf = par[5];
	float k = x[0];
	if (k <= kf) return n0 - alfa_minus*pow(k/kf, beta_minus);
	else return alfa_plus * pow(k/kf, beta_plus);
}

void lineshape(){
	auto *f = new TF1("f", n_k, 0, 5147, 6);
        f->SetParameters(0.85, 0.144, 0.134, -6.21, 3.24, 0.592*5147); // momentum in eV
	auto *h = new TH1F("nobes", "h", 20000, 200, 220);
        float me = 5.11e-4;
        f->Draw();
        new TCanvas();
        for (int i=0; i<(int)1e6; i++){
                if (i % 1000 == 0) cout << i << endl;
		float k = f->GetRandom()*1e-9;
                float costh = gRandom->Uniform(-1, 1);
                float ecms = sqrt(2*me*me + 2*sqrt(me*me+43.7*43.7)*sqrt(me*me+k*k) - 2*k*43.7*costh)*1000;
		h->Fill(ecms);
        }
        h->Draw();
        new TCanvas();
        auto *hh = new TH1F("bes", "h", 20000, 200, 220);
        for (int i=0; i<(int)1e6; i++){
                if (i % 1000 == 0) cout << i << endl;
                float k = f->GetRandom()*1e-9;
                float costh = gRandom->Uniform(-1, 1);
                float eb = 43.7*gRandom->Uniform(1-1.2e-2, 1+1.2e-2);
                float ecms = sqrt(2*me*me + 2*sqrt(me*me+eb*eb)*sqrt(me*me+k*k) - 2*k*eb*costh)*1000;
                hh->Fill(ecms);
        }
        hh->Draw();
}
