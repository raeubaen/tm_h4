void tgt_trk_draw_h4(){
  auto *g = new TGraph(2, (double[2]){-1, 65}, (double[2]){-120, 120});
  g->Draw();
  for (int i=0; i<10; i++){
    for (int j=0; j<4; j++){
      TLine *l = new TLine(1.2*i+0.2*j, -2, 1.2*i+0.2*j, 2);
      l->SetLineWidth(3);
      l->Draw();
    }
  }
  for (int i=0; i<10; i++){
    for (int j=0; j<2; j++){
      TLine *l = new TLine(0.8+1.2*i+0.2*j, -4.5, 0.8+1.2*i+0.2*j, 4.5);
      l->SetLineWidth(3);
      l->SetLineColor(kBlue);
      l->Draw();
    }
  }
}
