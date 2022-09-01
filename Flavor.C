// Purpose: class to facilitate retrieving flavor responses and fractions
#include "Flavor.h"

double Flavor::getResp(double pt, double absetamin, double absetamax,
		       double f[5], double w, string s) {

  // eta bins used in L5Flavor / L2Relative_flavor
  // cat <file> | | awk '{print $1", "}'
  const double etabins[] =
    {-5.191, -3.489, -3.139, -2.853, -2.5, -2.322, -1.93, -1.653, -1.305,-0.783,
     0, 0.783, 1.305, 1.653, 1.93, 2.322, 2.5, 2.853, 3.139, 3.489, 5.191};
  const int neta = sizeof(etabins)/sizeof(etabins[0])-1;

  double sumw(0), sumwjes(0);
  for (int i = 0; i != neta; ++i) {
    double eta = 0.5*(etabins[i]+etabins[i+1]);
    if (fabs(eta)>=absetamin && fabs(eta)<=absetamax &&
	(!doRobin || eta>0)) {

      double jes = getResp(pt, eta, f, w, s);
      double etaw = 1;
      sumw += etaw;
      sumwjes += etaw*jes;
    }
  } // for i
  assert(sumw>0);

  return (sumwjes / sumw);
} // Flavor::getResp 

double Flavor::getResp(double pt, double absetamin, double absetamax,
		       string smix, double w, string s) {
  // eta bins used in L5Flavor / L2Relative_flavor
  // cat <file> | | awk '{print $1", "}'
  const double etabins[] =
    {-5.191, -3.489, -3.139, -2.853, -2.5, -2.322, -1.93, -1.653, -1.305,-0.783,
     0, 0.783, 1.305, 1.653, 1.93, 2.322, 2.5, 2.853, 3.139, 3.489, 5.191};
  const int neta = sizeof(etabins)/sizeof(etabins[0])-1;

  double sumw(0), sumwjes(0);
  for (int i = 0; i != neta; ++i) {
    double eta = 0.5*(etabins[i]+etabins[i+1]);
    if (fabs(eta)>=absetamin && fabs(eta)<=absetamax &&
	(!doRobin || eta>0)) {

      double jes = getResp(pt, eta, smix, w, s);
      double etaw = 1;
      sumw += etaw;
      sumwjes += etaw*jes;
    }
  } // for i
  assert(sumw>0);

  return (sumwjes / sumw);
} // Flavor::getResp

double Flavor::getResp(double pt, double eta,
		       string smix, double w, string s) {
  double f[5];
  getFracs(pt, eta, smix,f);
  return getResp(pt, eta, f, w, s);
} // Flavor::getResp smix

// w=0: Pythia8 only, w=-1: Data only, w=1: Herwig7 only
double Flavor::getResp(double pt, double eta, double f[5],
		       //double fud, double fs, double fc, double fb,
		       double w, string s) {

  assert(fabs(eta)<1.3);

  // Load JEC on first call
  if (_mjec.empty()) loadJEC("Pythia8");
  if (_mhw7.empty()) loadJEC("Herwig7");
  if (_mres.empty()) loadRES();

  // Store weights into array for easy handling
  const int nf(5);
  string sref = "qcd";
  string sf[nf] = {"ud","s","c","b", "g"};
  double wf[nf] = {f[0],  f[1], f[2], f[3], f[4]};

  // Unitarize weights just to be sure they add up
  double sumwf(0);
  for (int i = 0; i != nf; ++i) sumwf += wf[i];
  for (int i = 0; i != nf; ++i) wf[i] /= sumwf;

  // Calculate reference JEC
  FactorizedJetCorrector *jref = _mjec["all"]; assert(jref);
  jref->setJetPt(pt);
  jref->setJetEta(eta);
  double jec0 = jref->getCorrection();
  
  // Calculate alternative Herwig7 JEC for QCD
  FactorizedJetCorrector *jhw7 = _mhw7["all"]; assert(jhw7);
  jhw7->setJetPt(pt);
  jhw7->setJetEta(eta);
  double jechw7 = jhw7->getCorrection();

  // Calculate weighted JEC
  double jes(0);
  for (int i = 0; i != nf; ++i) {

    FactorizedJetCorrector *fjc = _mjec[sf[i]]; assert(fjc);
    FactorizedJetCorrector *fhw = _mhw7[sf[i]]; assert(fhw);
    TF1 *f1 = _mres[sf[i]+s]; assert(f1);
    if (wf[i]!=0) {
      fjc->setJetPt(pt);
      fjc->setJetEta(eta);
      fhw->setJetPt(pt);
      fhw->setJetEta(eta);
      // Pythia8 JEC only
      if (w==0) {
	if (doRobin) jes += wf[i] * (fjc->getCorrection());
	else         jes += wf[i] * (jec0/fjc->getCorrection());

      }
      // Data JEC (w==-1)
      if (w<0)
	jes += wf[i] * ((1+w)*(jec0/fjc->getCorrection()) - w*f1->Eval(pt));
      // Herwig JEC (w==1)
      if (w>0)
	jes += wf[i] * ((1-w)*(jec0/fjc->getCorrection()) +
			w*(jechw7/fhw->getCorrection()));
    } // wf[i]!=0
  } // for i

  return jes;
} // Flavor::getResp f[5]

void Flavor::getFracs(double pt, double eta, string smix, double (&f)[5]) {

  assert(smix=="MultijetLeading13" || smix=="MultijetRecoil25" ||
	 smix=="EMJet13" || smix=="PhotonJet13" || smix=="ZJet13" ||
	 smix=="ud" || smix=="s" || smix=="c" || smix=="b" || smix=="g");
  assert(fabs(eta)<1.3);
  
  f[0]=f[1]=f[2]=f[3]=f[4]=0;

  if (smix=="ud") { f[0]=1; f[1]=f[2]=f[3]=f[4]=0; }
  if (smix=="s")  { f[1]=1; f[0]=f[2]=f[3]=f[4]=0; }
  if (smix=="c")  { f[2]=1; f[0]=f[1]=f[3]=f[4]=0; }
  if (smix=="b")  { f[3]=1; f[0]=f[1]=f[2]=f[4]=0; }
  if (smix=="g")  { f[4]=1; f[0]=f[1]=f[2]=f[3]=0; }

  if (smix=="MultijetRecoil25") {

    const int nf = 5;
    // Parameters from minitools/drawMultijetFlavor.C
    double p[nf][3] = {
      //double pmjrecoil[nf][3] = {
      // [p0]+[p1]*pow(x,[p2])
      {0.097845, 0.0017581, 0.65565}, // (ud 102.2/75)
      {0.057069, -0.003282, 0.3333}, // (s 51.0/72)
      {0.052562, 0.0052541, 0.3333}, // (c 31.2/72)
      {0.029487, 0.0019683, 0.3333}, // (b 51.7/72)
      {0.75931, -0.0034231, 0.58796}}; // (g 56.7/75)
    
    double sum(0);
    for (int i = 0; i != nf; ++i) {
      f[i] = max(0.,min(1.,p[i][0] + p[i][1]*pow(pt,p[i][2])));
      sum += f[i];
    }
    for (int i = 0; i != nf; ++i)
      f[i] /= sum;
  } // MultijetRecoil25

  if (smix=="MultijetLeading13") {

    const int nf = 5;
    // Parameters from minitools/drawMultijetFlavor.C
    double p[nf][3] = {
      //double pmjlead[nf][3] = {
      // [p0]+[p1]*pow(x,[p2])
      {0.020545, 0.0068051, 0.64066}, // (ud 93.1/75)
      {0.064103, -0.003282, 0.3333}, // (s 106.6/72)
      {0.071159, 0.0027745, 0.3333}, // (c 89.7/72)
      {0.026364, 0.0019683, 0.3333}, // (b 185.0/71)
      {0.89715, -0.023898, 0.47537}}; // (g 66.1/75)

    double sum(0);
    for (int i = 0; i != nf; ++i) {
      f[i] = max(0.,min(1.,p[i][0] + p[i][1]*pow(pt,p[i][2])));
      sum += f[i];
    }
    for (int i = 0; i != nf; ++i)
      f[i] /= sum;
  } // MultijetLeading13

  // Decent fits to EMJet13, could be better for ud & g
  if (smix=="EMJet13") {
    
    const int nf = 5;
    // Parameters from gamjet/drawPurityEstimates.C
    double p[nf][4] = {
      // [p0]+[p1]*pow(x,[p2])+[p3]/x
      {-0.033376, 0.012214, 0.53486, 3.7079}, // (ud 74.2/15)
      {0.076392, -0.0059514, 0.29532, -0.23598}, // (s 18.0/15)
      {-0.061882, 0.13575, 0.028896, -0.98938}, // (c 24.8/15)
      {0.049648, -9.007e-06, 0.91925, -1.0576}, // (b 25.1/15)
      {0.98949, -0.043615, 0.38531, -4.9996}}; // (g 41.4/15)

    double sum(0);
    for (int i = 0; i != nf; ++i) {
      f[i] = max(0.,min(1., p[i][0] + p[i][1]*pow(pt,p[i][2]) + p[i][3]/pt ));
      sum += f[i];
    }
    for (int i = 0; i != nf; ++i)
      f[i] /= sum;
  } // EMJet13

  // These are rough fits, bad fit and chi2 for c & g
  if (smix=="PhotonJet13") {

    const int nf = 5;
    // Parameters from gamjet/drawPurityEstimates.C (gamma+jet)
    double p[nf][5] = {
      // [p0]+[p1]*pow(x,[p2])+[p3]*pow(x,[p4])
      {0.75602, -0.027673, -0.41329, -1.1398, -0.42332}, // (ud 59.5/14)
      {0.25212, 6.1313e-06, 1.178, -0.087696, 0.1573}, // (s 61.9/14)
      {2.8802, -0.0083099, 0.002725, -2.5256, 0.013656}, // (c 1622.5/14)
      {0.77408, 0.0012412, 0.49022, -0.68889, 0.018864}, // (b 66.4/14)
      {0.77766, -7.343e-05, 0.90022, -0.9518, -0.073982}}; // (g 487.4/14)

    double sum(0);
    for (int i = 0; i != nf; ++i) {
      f[i] = max(0.,min(1., p[i][0] + p[i][1]*pow(pt,p[i][2]) +
			p[i][3]*pow(pt,p[i][4]) ));
      sum += f[i];
    }
    for (int i = 0; i != nf; ++i)
      f[i] /= sum;
  } // PhotonJet13

  // ud, g could be a bit better, but otherwise decent fits with sum = 1 +/- 1%
  if (smix=="ZJet13") {

    // Parameters from minitools/drawZJetFlavor.C
    const int nf = 5;
    double p[nf][6] = {
      // [p0]+[p1]*x/(x+[p2])+[p3]*x/(x+[p4])+[p5]/x
      { -1.8475, 0.19962, 6.7418,  2.2337,  9.046,  9.6497}, // (ud 74.4/14)
      {-0.63151,  2.0528, 26.039, -1.3965, 54.643,   4.557}, // (s 14.8/14)
      {-0.12862, 0.99525, 32.874, -0.8012, 50.882,  1.0251}, // (c 18.7/13)
      {-0.11083, 0.52386, 41.228, -0.4095, 105.48, 0.75193}, // (b 24.7/13)
      {  3.9934, -7.9394, 21.159,  4.2641, 43.495, -20.404}}; // (g 39.2/13)
    
    double sum(0);
    double x = max(20.,min(350.,pt));
    for (int i = 0; i != nf; ++i) {
      f[i] = max(0.,min(1., p[i][0]+p[i][1]*x/(x+p[i][2])
			+ p[i][3]*x/(x+p[i][4]) + p[i][5]/x));
      sum += f[i];
    }
    for (int i = 0; i != nf; ++i)
      f[i] /= sum;
  } // ZJet13
} // getFracs

void Flavor::loadJEC(string smc) {
  
  assert(smc=="Pythia8" || smc=="Herwig7");

  const char *cm = smc.c_str();
  cout << Form("Flavor::loadJEC(\"%s\")...",cm) << endl << flush; 

  const int nf(5+1);
  string sf[nf] = {"ud","s","c","b","g","all"};

  const char *cs = "%s/Autumn18_V3_MC_%s_%s_L2Relative_AK4PFchs.txt";
  const char *cd = "rootfiles/flavor";
  for (int i = 0; i != nf; ++i) {
    
    const char *cf = sf[i].c_str();
    string s = Form(cs,cd,cm,cf);
    cout << "  " << s << endl << flush;

    vector<JetCorrectorParameters> v;
    JetCorrectorParameters *p = new JetCorrectorParameters(s);
    v.push_back(*p);
    FactorizedJetCorrector *j = new FactorizedJetCorrector(v);

    if (smc=="Pythia8") _mjec[cf] = j;
    if (smc=="Herwig7") _mhw7[cf] = j;
  } // for i

} // Flavor::loadJEC

void Flavor::loadRES() {
  
  const int nf = 5;
  string vf[] = {"ud","s","c","b","g"};

  if (false) { // old drawZflavor.C results instead of Zflavor.C
  // minitools/drawZflavor.C HDM results
  /*
  double ph[nf][3] = {
    //{-5.7333e-15, 0, 0},  // i (0.0/9)
    {0.66088, -1, 0},       // q (14.5/8)
    {0.66088, -1, 0},       // q (14.5/8) => copy for s
    {0.13851, 0, 0},        // c (24.6/9)
    {0.1541, 0, 0},         // b (18.6/9)
    {-1.5505, -0.6583, 0}}; // g (10.9/8)
  */
  /*
  double ph[nf][3] = { // variant with q,g slope fixed to -1
    //{-5.7333e-15, 0, 0}, // i (0.0/9) HDM
    {0.66088, -1, 0}, // q (14.5/9) HDM
    {0.66088, -1, 0}, // q (14.5/9) HDM => copy for s
    {0.13851, 0, 0}, // c (24.6/9) HDM
    {0.1541, 0, 0}, // b (18.6/9) HDM
    {-1.2968, -1, 0}}; // g (14.7/9) HDM
  */
  double ph[nf][3] = { // variant with q,g slope fixed to -1, scaled by -50%
    //{-5.7333e-15, 0, 0}, // i (0.0/9) HDM
    {0.5*0.66088, -1, 0}, // q (14.5/9) HDM
    {0.5*0.66088, -1, 0}, // q (14.5/9) HDM => copy for s
    {0.13851, 0, 0}, // c (24.6/9) HDM
    {0.1541, 0, 0}, // b (18.6/9) HDM
    {0.5*-1.2968, -1, 0}}; // g (14.7/9) HDM

  // minitools/drawZflavor.C MPF results
  double pm[nf][3] = {
    //{-5.7333e-15, 0, 0},     // i (0.0/9) MPF
    {0.10385, -1, 0},          // q (10.8/8) MPF
    {0.10385, -1, 0},          // q (10.8/8) MPF => copy for s
    {0.2903, 0, 0},            // c (11.2/9) MPF
    {0.38388, 0, 0},           // b (9.2/9) MPF
    {-0.48713, -0.085302, 0}}; // g (6.9/8) MPF

  double (&p)[nf][3] = ph;

  for (int i = 0; i != nf; ++i) {

    string sf = vf[i];
    const char *cf = sf.c_str();

    //TF1 *f1 = new TF1(Form("f1%s",cf),"[0]+[1]*pow(x,[2])",15,6500);
    //f1->SetParameters(1,0,0);
    //if (sf=="g")
    //f1->SetParameter(0,0.985);

    // minitools/drawZflavor.C HDM function, scaled from % difference
    // Cap pT dependence at 400 GeV, since data only extended to ~200 GeV?
    TF1 *f1 = new TF1(Form("f1%s",cf),
		      //"1+0.01*([2]+[0]*pow(0.01*min(400.,x),[1]))",
		      "1+0.01*([2]+[0]*pow(0.01*x,[1]))",
		      15,6500);
    f1->SetParameters(p[i][0],p[i][1],p[i][2]);
    _mres[sf] = f1;
  } // for i
  } // false
  else { // new functions from data

  // Flavor fit parameters from minitools/Zflavor.C (Eta13)
  TF1 *f1q = new TF1("f1q","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);
  TF1 *f1g = new TF1("f1g","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);
  TF1 *f1c = new TF1("f1c","1+0.01*[p0]",45,300);
  TF1 *f1b = new TF1("f1b","1+0.01*[p0]",45,300);
  TF1 *f1z = new TF1("f1z","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);

  TF1 *f1q3 = new TF1("f1q3","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);
  TF1 *f1g3 = new TF1("f1g3","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);
  TF1 *f1c3 = new TF1("f1c3","1+0.01*[p0]",45,300);
  TF1 *f1b3 = new TF1("f1b3","1+0.01*[p0]",45,300);
  TF1 *f1z3 = new TF1("f1z3","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);

  TF1 *f1q5 = new TF1("f1q5","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);
  TF1 *f1g5 = new TF1("f1g5","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);
  TF1 *f1c5 = new TF1("f1c5","1+0.01*[p0]",45,300);
  TF1 *f1b5 = new TF1("f1b5","1+0.01*[p0]",45,300);
  TF1 *f1z5 = new TF1("f1z5","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",45,300);

  if (true) {

    // Eta13
    /*
    f1q->SetParameters(0.7966, 0.9311, -1); // chi2/NDF=9.6/7
    f1g->SetParameters(-1.764, -1.194, -1); // chi2/NDF=12.5/7
    f1c->SetParameter(0,-0.6218); // chi2/NDF=25.2/8
    f1b->SetParameter(0,0.3651); // chi2/NDF=15.3/8
    f1z->SetParameters(0, 0.5, -1); // chi2/NDF=0.0/0
    */
    // New Run2Test values
    f1q->SetParameters(0.2746, 0.1031, -1); // chi2/NDF=20.3/7
    f1g->SetParameters(-1.07, -0.5904, -1); // chi2/NDF=14.3/7
    f1c->SetParameter(0,-0.2341); // chi2/NDF=16.1/8
    f1b->SetParameter(0,0.5346); // chi2/NDF=24.1/8
    f1z->SetParameters(0, 0.5, -1); // chi2/NDF=0.0/0

    _mres["ud"] = f1q;
    _mres["s"] = f1q;
    _mres["c"] = f1c;
    _mres["b"] = f1b;
    _mres["g"] = f1g;
  }
  if (true) {

    // Eta13
    /*
    f1q3->SetParameters(0.7966, 0.9311, -1); // chi2/NDF=9.6/7
    f1g3->SetParameters(-1.764, -1.194, -1); // chi2/NDF=12.5/7
    f1c3->SetParameter(0,-0.6218); // chi2/NDF=25.2/8
    f1b3->SetParameter(0,0.3651); // chi2/NDF=15.3/8
    f1z3->SetParameters(0, 0.5, -1); // chi2/NDF=0.0/0
    */
    // New Run2Test values
    f1q3->SetParameters(0.2746, 0.1031, -1); // chi2/NDF=20.3/7
    f1g3->SetParameters(-1.07, -0.5904, -1); // chi2/NDF=14.3/7
    f1c3->SetParameter(0,-0.2341); // chi2/NDF=16.1/8
    f1b3->SetParameter(0,0.5346); // chi2/NDF=24.1/8
    f1z3->SetParameters(0, 0.5, -1); // chi2/NDF=0.0/0

    _mres["ud13"] = f1q3;
    _mres["s13"] = f1q3;
    _mres["c13"] = f1c3;
    _mres["b13"] = f1b3;
    _mres["g13"] = f1g3;
  }
  if (true) {

    // Eta25
    /*
    f1q5->SetParameters(0.9382, 1.03, -1); // chi2/NDF=9.6/7
    f1g5->SetParameters(-2.264, -1.981, -1); // chi2/NDF=5.3/7
    f1c5->SetParameter(0, -0.1573); // chi2/NDF=19.8/8
    f1b5->SetParameter(0, 0.1104); // chi2/NDF=16.3/8
    f1z5->SetParameters(0, 0.5, -1); // chi2/NDF=0.0/0
    */
    // New Run2Test values
    f1q5->SetParameters(0.269, 0.0725, -1); // chi2/NDF=16.6/7
    f1g5->SetParameters(-1.37, -1.056, -1); // chi2/NDF=12.1/7
    f1c5->SetParameter(0,0.0836); // chi2/NDF=10.4/8
    f1b5->SetParameter(0,0.314); // chi2/NDF=23.7/8
    f1z5->SetParameters(0, 0.5, -1); // chi2/NDF=0.0/0

    _mres["ud25"] = f1q5;
    _mres["s25"] = f1q5;
    _mres["c25"] = f1c5;
    _mres["b25"] = f1b5;
    _mres["g25"] = f1g5;
  }
  
  } // else
} // Flavor::loadRES


