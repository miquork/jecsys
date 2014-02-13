

const double x_pt[] =
    {8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
     507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, //1999};
     2000, 2238, 2500, 2787, 3103, 3450};
const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;
const double x_eta[] =  {-5.4,-5.0,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5.0,5.4};
const int ndiv_eta = sizeof(x_eta)/sizeof(x_eta[0])-1;     


TH2D* readHist(const char* file, const char* title) {
  double eta1, eta2,pt, val1, val2;
  int nbins;

  TH2D* hist = new TH2D(title,title,ndiv_pt, x_pt,ndiv_eta,x_eta);
  hist->SetXTitle("p_{T} [GeV]");
  hist->SetYTitle("#eta");
  hist->SetStats(0);
  hist->GetXaxis()->SetMoreLogLabels();
  hist->GetXaxis()->SetNoExponent();
  hist->SetMinimum(0);
  std::cout << "opening file " << file << '\n';
  std::ifstream f(file);
  if(! f.is_open()) return 0;
  char buff[256];
  f.getline(buff,256);
  for(int i = 0 ; i < ndiv_eta ; ++i) {
    f >> eta1 >> eta2 >> nbins;
    nbins /=3;
    assert(nbins == ndiv_pt);
    for(int j = 0 ; j < nbins ; ++j) {
      f >> pt >> val1 >> val2;
      //std::cout << pt << ", " << eta1 << ", " << eta2 << ", " << val1 << '\n';
      hist->Fill(pt,0.5*(eta1+eta2),val1);
    }    
  }
  
  f.close();

  return hist;

}


void drawUncertainty() {
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleOffset(1.10,"y");
  gStyle->SetOptTitle(1);

  std::vector<std::string> list;
  
  list.push_back("DATA_Uncertainty_AK5PFchs.txt");
  list.push_back("DATA_Uncertainty_AK5PF.txt");
  list.push_back("DATA_Uncertainty_AK7PFchs.txt");
  list.push_back("DATA_Uncertainty_AK7PF.txt");
  list.push_back("DATA_Uncertainty_AK5Calo.txt");
  list.push_back("DATA_Uncertainty_AK7Calo.txt");
  //list.push_back("MC_Uncertainty_AK5PFchs.txt");
  //list.push_back("MC_Uncertainty_AK5PF.txt");
  //list.push_back("MC_Uncertainty_AK7PFchs.txt");
  //list.push_back("MC_Uncertainty_AK7PF.txt");
  //list.push_back("MC_Uncertainty_AK5Calo.txt");
  //list.push_back("MC_Uncertainty_AK7Calo.txt");


  std::string newJEC = "Legacy11_V1";
  std::string oldJEC = "Jec11_V13";
  for(int i = 0 ; i < list.size() ; ++i) {
    TCanvas* c = new TCanvas();
    std::string file("txt/");
    file += newJEC + "_" + list[i];
    unsigned begin = file.find_last_of('/');
    unsigned end = file.find_last_of('.');
    std::string name = file.substr(begin+1,end-begin-1);
    TH2D* h2 = readHist(file.c_str(),name.c_str());
    h2->DrawClone("COLZ");
    c->SetLogx();
    std::string pname = name+".eps";
    c->Print(pname.c_str());

    std::string oldfile("txt/");
    oldfile += oldJEC + "_" + list[i];
    begin = file.find_last_of('/');
    end = file.find_last_of('.');
    std::string oldname = file.substr(begin+1,end-begin-1);

    TH2D* h2old = 0;//readHist(oldfile.c_str(),oldname.c_str());
    if(h2old) {
      TCanvas* c = new TCanvas();
      h2old->DrawClone("COLZ");
      c->SetLogx();
      pname = oldname+".eps";
      c->Print(pname.c_str());
      pname = name+"diff";
      TH2D* hdiff = (TH2D*)h2old->Clone(pname.c_str());
      pname = name+" difference (old -new)";
      hdiff->SetTitle(pname.c_str());
      hdiff->Add(h2, -1.0);
      //h2old->Scale(-1.0);
      c = new TCanvas();
      //hdiff->SetMaximum(0.2);
      //hdiff->SetMinimum(-0.0);
      hdiff->DrawClone("COLZ");
      c->SetLogx();
      //c->SetLogz();
      pname = name+"diff"+".eps";
      c->Print(pname.c_str());
      c = new TCanvas();
      hdiff->Divide(h2old);
      pname = name+" difference (old -new)/old";
      hdiff->SetTitle(pname.c_str());
      //hdiff->SetMaximum(0.2);
      //hdiff->SetMinimum(-0.1);
      hdiff->DrawClone("COLZ");
      c->SetLogx();
      pname = name+"diffrel"+".eps";
      c->Print(pname.c_str());
      delete h2old;
      delete hdiff;
    }

    delete h2;
  }
  /*
  TH2D* hold = readHist("txt/Summer13_V1_DATA_Uncertainty_AK5Calo.txt","AK5Calo old");
  TH2D* hnew = readHist("txt/Summer13_V4_DATA_Uncertainty_AK5Calo.txt","differences (old-new)");
  
  TCanvas* c = new TCanvas();
  hold->DrawClone("COLZ");
  c->SetLogx();
  c->Print("oldAK5Calo.eps");

  TH2D* hdiff = (TH2D*)hold->Clone("hdiffcalo");
  hdiff->SetTitle("AK5Calo differences (old-new)");
  hdiff->Add(hnew, -1.0);
  //hold->Scale(-1.0);
  c = new TCanvas();
  hdiff->SetMaximum(0.2);
  hdiff->SetMinimum(-0.0);
  hdiff->DrawClone("COLZ");
  c->SetLogx();
  //c->SetLogz();
  c->Print("diffchsv3calo.eps");
  c = new TCanvas();
  hdiff->Divide(hold);
  hdiff->SetTitle("AK5Calo differences (old-new)/old");
  hdiff->SetMaximum(0.2);
  hdiff->SetMinimum(-0.1);
  hdiff->DrawClone("COLZ");
  c->SetLogx();
  c->Print("diffchsv3calorel.eps");
  */
}
