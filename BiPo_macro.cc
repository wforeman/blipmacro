  //////////////////////////////////////////////////////////////////
  // 
  //  Analysis ROOT Macro
  //
  ////////////////////////////////////////////////////////////////// 
  #include "core/tools.h"
  #include "core/vars.h"
  #include "core/bipo.h"
  #include <time.h>
  
  // --- Choose run configuration ---
  int   fConfig   = 0;
  float fMinHr    = 0; 
  float fMaxHr    = 40e9; 
  
  float chargeScaling = 1.0; //0.905; 
  TRandom2* fRand;

  // --- Input files ---
  infile_t fInputFiles[4] = {
  /*0*/  { "BlipAna_20221017_BiPo_Overlay_QY75.root", "blipanaTrkMask/anatree",true,0,0},
         //{ "BlipAna_20221017_BiPo_MC.root", "blipana/anatree",true,0,0},
         //{ "BlipAna_20221017_Data_Overlay.root","blipanaTrkMask/anatree",false,0,0},
  /*1*/  { "BlipAna_20221017_RadonData_Phase1.root","blipanaTrkMask/anatree",false,1627415210,1627592728},
  /*2*/  { "BlipAna_20221017_RadonData_Phase2.root","blipanaTrkMask/anatree",false,1627594380,1627761265},
  /*3*/  { "BlipAna_20221017_Run3.root","blipanaTrkMask/anatree",false,1528526500,1532448800}
  };
  
  float fRecomb       = 0.584; 

  const float binPeriodHrs[4] = { 2,  2,   2,   48 };
  const float binPeriodMax[4] = { 44, 44, 44,  1104 };

  int   fCaloMode         = 0;

  // --- General selection options ---
  bool  fFidVolCut        = 1; //1;      // Fiducialize beta
  bool  fPickyBetaMode    = 1;      // Require beta blips match on all 3 planes
  int   fBetaWires_max    = 4;
  float fBetaCharge_min   = 12e3;   // Min charge of beta candidate blip [e-]
  float fBetaCharge_max   = 90e3;   // Max charge of beta candidate blip [e-]
  float fBetaEnergy_min   = 0.5;    // 
  float fBetaEnergy_max   = 3.5;    
  int   fWireRange        = 0;      // +/- range to look for alpha candidate;
  int   fAlphaWires_max   = 2;
  float fAlphaCharge_min  = 0e3;    // Min charge of alpha candidate cluster [e-]
  float fAlphaCharge_max  = 6e3;    // Max charge of alpha candidate cluster [e-]
  float fAlphaEnergy_min  = 0;
  float fAlphaEnergy_max  = 0.15; //0.24; //0.24 nominal;   //
  bool  fLinearizeCorr    = 1;      // Linearize detector effects control facto
  float fdT_binSize       = 20.;    // Bin width for all dT spectra plots [us]
  float fdT_min           = 60.;    // Min dT for looking for candidate [us]
  float fdT_max           = 340.;   // Max dT for looking for candidate [us]
  
  // --- MC efficiency for equiv activity calc ---
  double fEfficiencyMC      = 0.11;
  double fEfficiencyMC_err  = 0.03;

  // --- Special MC options ---
  bool  fIgnoreTrueAlphas = false;
  bool  fIgnoreTrueGammas = false;
  bool  fIgnoreNonMC      = false; //true;

  // --- Detector properties ---
  int   nWiresColl        = 3455;
  float fSamplePeriod     = 0.5; // microseconds
  float fZlim[2]          = {50,985}; // Z range (0 to 1037 cm)
  float fYlim[2]          = {-80,80}; // Y range (-120 to 120 cm)
  
  // --- Special switches ---
  int   fRandomWireShift  = 0; 
  int   fBackgroundMode   = 1;      // 0= wire-shift, 1= dT-flip
  
  // --- Noisy wires to skip (collection plane) ---
  std::vector<int> fNoisyWires{
  0, 1, 2, 3, 6, 7, 83, 384, 1151, 1247, 1535, 1540, 1919, 1920, 2303, 2311, 2334, 2335, 
  2400, 2415, 2687, 2733, 2753, 2783, 2879, 3071, 3072, 3215, 3263, 3274, 3286, 3299, 
  3318, 3327, 3385, 3391, 3408, 3409, 3410, 3411, 3412, 3413, 3414, 3415, 3416, 3417, 
  3418, 3419, 3420, 3421, 3422, 3423, 3424, 3425, 3426, 3427, 3428, 3429, 3430, 3431, 
  3432, 3433, 3434, 3435, 3436, 3437, 3438, 3439, 3440, 3441, 3442, 3443, 3444, 3445, 
  3446, 3447, 3448, 3449, 3451, 3452, 3453, 3454, 3455
  };

  // --- Vector of bad/missing wires (collection plane)
  std::vector<int>fBadWires{
  18, 19, 37, 38, 39, 71, 82, 84, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 
  187, 188, 189, 190, 191, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 
  349, 350, 351, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 
  415, 432, 433, 434, 435, 436, 437, 438, 439, 440, 442, 443, 444, 446, 527, 528, 754, 756, 
  816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 961, 962, 
  963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 1376, 1377, 1378, 1379, 1380, 1381, 1382, 
  1383, 1384, 1385, 1386, 1387, 1388, 1389, 1390, 1391, 1592, 2227, 2231, 2336, 2337, 2338, 
  2339, 2340, 2341, 2342, 2343, 2344, 2345, 2346, 2347, 2348, 2349, 2350, 2351, 2352, 2353, 
  2354, 2355, 2356, 2357, 2358, 2359, 2360, 2361, 2362, 2363, 2364, 2365, 2366, 2367, 2368, 
  2369, 2370, 2371, 2372, 2373, 2374, 2375, 2376, 2377, 2378, 2379, 2380, 2381, 2382, 2383, 
  2384, 2385, 2386, 2387, 2388, 2389, 2390, 2391, 2392, 2393, 2394, 2395, 2396, 2397, 2398, 
  2399, 2401, 2402, 2403, 2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412, 2413, 2414, 
  2416, 2417, 2418, 2419, 2420, 2421, 2422, 2423, 2424, 2425, 2426, 2427, 2428, 2429, 2430, 
  2431, 2432, 2433, 2434, 2435, 2436, 2437, 2438, 2439, 2440, 2441, 2442, 2443, 2444, 2445, 
  2446, 2447, 2448, 2449, 2450, 2451, 2452, 2453, 2454, 2455, 2456, 2457, 2458, 2459, 2460, 
  2461, 2462, 2463, 2625, 2688, 2689, 2690, 2691, 2692, 2693, 2694, 2695, 2696, 2697, 2698, 
  2699, 2700, 2701, 2702, 2703, 2736, 2737, 2738, 2739, 2740, 2741, 2742, 2743, 2744, 2745, 
  2746, 2747, 2748, 2749, 2750, 2751, 2912, 2913, 2914, 2915, 2916, 2917, 2918, 2919, 2920, 
  2921, 2922, 2923, 2924, 2925, 2926, 2927, 3392, 3393, 3394, 3395, 3396, 3397, 3398, 3399, 
  3400, 3401, 3402, 3403, 3404, 3405, 3406, 3407
  };
  


  //#######################################################################
  // Derived parameters
  //######################################################################
  // Fiducial vol correction factor
  float dz = fZlim[1]-fZlim[0];
  float dy = fYlim[1]-fYlim[0];
  float _fiducialFrac = (fFidVolCut) ? std::min(1., (dz*dy)/(1037.*233.) ) : 1.0;
  // Counters / maps / etc
  bool  _isMC           = false;
  int   _numEvents      = 0;
  int   _numBiPo        = 0;
  int   _numBiPo_mcmatch   = 0;
  int   _numBiPo_true_perfectReco = 0;
  std::vector<bool> _clustAvailable;
  std::map<int,std::vector<int>> _map_wire_clusters;
  std::vector<bool> wireIsNoisy (nWiresColl,false);
  std::vector<bool> wireIsBad   (nWiresColl,false);

  // Live time
  int   _minTick         = 0;
  int   _maxTick         = 6400 - (int)fdT_max*2;
  float _liveTimePerEvt  = _maxTick*fSamplePeriod*1e-6; //sec
  float _totalLiveTime  = 0;
  // Blip fidelity requirements
  int _betaMinPlanes = 2;
  int _betaMaxDiff = 9999;
      
  std::vector<bool> clust_isAlpha; //(nclusts,false);
  std::vector<bool> clust_isBeta; //(nclusts,false);
  std::vector<bool> clust_isGamma; //(nclusts,false);

  TF1* f_backward_corr = new TF1("backward_corr","[0]+[1]*x",fdT_min,fdT_max);

  //##########################################################################
  // Functions and ROOT objects
  //##########################################################################
  void                        makePlots();
  void                        makeHistograms();
  void                        setRootStyle();
  std::vector<BiPoCandidate>  FindCandidates(int, int, int, bool, int&, int&, int&);
  FitResult                   fitdT(TH1D*,bool,bool);
  //int                         FindG4Index(int);
  //double                      calcActivity(double,double);
  
  // ROOT objects
  TTree*      fTree;
  TFile*      fOutFile;

  // Histograms
  TDirectory* tdir_util;
  TDirectory* tdir_plots;
 
  TH1D* h_cuts;

  TH1D* h_nclusts_inwindow;
  TH1D* h_nclusts_20us;
  TH1D* h_ncands_inwindow;
  
  TH1D* h_cand_dT;
  TH1D* h_cand_dT_shift;
  TH1D* h_cand_dT_flip;
  TH1D* h_cand_dT_sub;
  
  TH1D* h_control_dT;
  TH1D* h_control_dT_flip;
  TH1D* h_control_dT_ratio;
  TH1D* h_control_ratio;

  TH2D* h_wt_clusts;
  TH2D* h_wt_blips;
  TH2D* h_wt_blips_filt;
  TH2D* h_wt_bipos;
  TH2D* h_zy_bipos;
  TH2D* h_zy_bipos_bg;
  TH2D* h_zy_bipos_sub;

  TH1D* h_time_vs_N;        
  TH2D* h_2D_time_vs_dT;
  TH2D* h_2D_time_vs_dT_bg;

  TH1D* h_time_vs_rate; 
  TH1D* h_time_vs_rate_bg;
  TH1D* h_time_vs_activity; 
 
  TH1D* h_alpha_charge;
  TH1D* h_alpha_charge_bg;
  TH1D* h_alpha_charge_sub;
  TH1D* h_beta_charge;
  TH1D* h_beta_charge_bg;
  TH1D* h_beta_charge_sub;
  
  TH1D* h_alpha_energy;
  TH1D* h_alpha_energy_bg;
  TH1D* h_alpha_energy_sub;
  TH1D* h_beta_energy;
  TH1D* h_beta_energy_bg;
  TH1D* h_beta_energy_sub;
 
  TH2D* h_alpha_energyVsdT;

  TH1D* h_true_alpha_depne;
  TH1D* h_true_alpha_charge;
  TH1D* h_matched_alpha_charge;
  
  TH1D* h_beta_trueEnergy;
  TH1D* h_beta_trueEnergySum;
  TH1D* h_beta_trueEnergySum_reco;
  TH1D* h_beta_trueEnergySum_recoCuts;
  TH1D* h_beta_nwires;
  TH1D* h_alpha_nwires;
  
  //##########################################################################
  // Initialize histograms
  //##########################################################################
  void makeHistograms()
  {
    fOutFile->cd();
    
    //fOutTree = new TTree("bipoTree","bipoTree");
    //fOutTree ->SetBranchAddress("event",);
    //fOutTree ->SetBranchAddress("event_hr");

    //int ncuts = 5;
    //h_cuts = new TH1D("cuts","Blip cut analysis",
    
    // E = (Q/R) * 23.6E-6 MeV
    float alphaQmax   = fAlphaCharge_max; 
    int   alphaQbins  = alphaQmax/200.;
    float betaQmax    = fBetaCharge_max;
    int   betaQbins   = betaQmax/2000;
    float alphaEmax   = fAlphaEnergy_max;
    int   alphaEbins    = fAlphaEnergy_max/0.01;
    float betaEmax    = fBetaEnergy_max;
    int   betaEbins   = fBetaEnergy_max/0.10;

    h_beta_nwires = new TH1D("beta_nwires","True beta cluster;Number of collection plane wires;Entries",10,0,10);
    h_alpha_nwires = new TH1D("alpha_nwires","True alpha cluster;Number of collection plane wires;Entries",10,0,10);

    h_beta_trueEnergy             = new TH1D("beta_trueEnergy",             "All true decays;Electron true energy [MeV];Entries per bin",70,0,betaEmax);
    h_beta_trueEnergySum          = new TH1D("beta_trueEnergySum",          "All true decays;Decay vertex true energy [MeV];Entries per bin",70,0,betaEmax);
    h_beta_trueEnergySum_reco     = new TH1D("beta_trueEnergySum_reco",     "Reco'd on coll plane;Decay vertex true energy [MeV];Entries per bin",70,0,betaEmax);
    h_beta_trueEnergySum_recoCuts = new TH1D("beta_trueEnergySum_recoCuts", "blip cuts;Decay vertex true energy [MeV];Entries per bin",70,0,betaEmax);
    
    //h_beta_trueEnergy = new TH1D("beta_trueEnergy","True #beta energy;Decay vertex true energy [MeV];Entries per bin",175,0,betaEmax);
    //h_beta_trueEnergy_reco = new TH1D("beta_trueEnergy_reco","True #beta energy;Decay vertex true energy [MeV];Entries per bin",175,0,betaEmax);

    //h_alpha_nwires       = new TH1D("alpha_nwires","Candidate alpha clusters (coll plane);Wire extent",20,0,20);
    h_nclusts_inwindow  = new TH1D("nclusts_inwindow","Mean clusters per wire in time window following Bi-candidate",50,0,5);
    h_nclusts_20us      = new TH1D("nclusts_20us","Mean clusters per wire <20#mus following Bi-candidate",50,0,5);
    h_ncands_inwindow   = new TH1D("ncands_inwindow","Number of Po candidates in time window following Bi-candidate",10,0,10);
    
    float Zmin = -100;  float Zmax = 1100;  int Zbins = 120;
    float Ymin = -150;  float Ymax = 150;   int Ybins = 30;
    float Tmin = -1000;     float Tmax = 6000;  int Tbins = 700; 
    float Wmin = -100;  float Wmax = 3500;  int Wbins = 1800; 
    h_zy_bipos      = new TH2D("zy_bipos","BiPo candidates;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_bipos      ->SetOption("colz");
    h_zy_bipos_bg      = new TH2D("zy_bipos_bg","Background BiPo candidates;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_bipos_bg      ->SetOption("colz");
    h_zy_bipos_sub      = new TH2D("zy_bipos_sub","Backgrounds-subtracted BiPo candidates;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_bipos_sub      ->SetOption("colz");
    
    h_wt_clusts     = new TH2D("wt_clusts","2D clusts;Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_blips      = new TH2D("wt_blips","3D blips;Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_blips_filt = new TH2D("wt_blips_filt","3D blips (quality cuts);Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_bipos      = new TH2D("wt_bipos","BiPo candidates;Collection Plane Wire;Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_clusts     ->SetOption("colz");
    h_wt_blips      ->SetOption("colz");
    h_wt_blips_filt ->SetOption("colz");
    h_wt_bipos      ->SetOption("colz");
   
    int dTbins = fdT_max / fdT_binSize;
    h_cand_dT       = new TH1D("cand_dT","Selected BiPo Candidates;#DeltaT [#mus];Candidates per second / 20 #mus", dTbins,0.,fdT_max);
    h_cand_dT_shift = (TH1D*)h_cand_dT->Clone("cand_dT_shift"); h_cand_dT_shift ->SetTitle("Shifted-wire region candidates");
    h_cand_dT_flip  = (TH1D*)h_cand_dT->Clone("cand_dT_flip");  h_cand_dT_flip  ->SetTitle("Opposite dT candidates");
    h_cand_dT_sub   = (TH1D*)h_cand_dT->Clone("cand_dT_sub");   h_cand_dT_sub   ->SetTitle("Background-subtracted spectrum");
   
    h_control_dT        = new TH1D("control_dT","OFFSET REGION;Time difference [#mus];Number of candidates", dTbins,0.,fdT_max);
    h_control_dT_flip   = new TH1D("control_dT_flip","OFFSET REGION;Time difference [#mus];Number of candidates", dTbins,0.,fdT_max);
    h_control_dT_ratio  = new TH1D("control_dT_ratio","OFFSET REGION;Time difference [#mus];Number of candidates", dTbins,0.,fdT_max);

    h_alpha_charge      = new TH1D("alpha_charge","Candidate alphas;Collected charge [e^{-}];Entries per second", alphaQbins, 0, alphaQmax);
    h_alpha_charge_bg   = (TH1D*)h_alpha_charge->Clone("alpha_charge_bg");
    h_alpha_charge_sub  = (TH1D*)h_alpha_charge->Clone("alpha_charge_sub");
    h_alpha_charge_sub  ->SetTitle("Candidate alphas after background subtraction");
    h_beta_charge       = new TH1D("beta_charge","Candidate betas;Collected charge [e^{-}];Events", betaQbins, 0, betaQmax);
    h_beta_charge_bg    = (TH1D*)h_beta_charge->Clone("beta_charge_bg");
    h_beta_charge_sub   = (TH1D*)h_beta_charge->Clone("beta_charge_sub");
    h_beta_charge_sub   ->SetTitle("Candidate betas after background subtraction");
   
    h_alpha_energy      = new TH1D("alpha_energy","Candidate alphas;Electron-equivalent energy [MeVee];Entries per second", alphaEbins, 0, alphaEmax);
    h_alpha_energy_bg   = (TH1D*)h_alpha_energy->Clone("alpha_energy_bg");
    h_alpha_energy_sub  = (TH1D*)h_alpha_energy->Clone("alpha_energy_sub");
    h_alpha_energy_sub  ->SetTitle("Candidate alphas after background subtraction");
    h_beta_energy       = new TH1D("beta_energy","Candidate betas;Energy [MeV];Entries per second", betaEbins, 0, betaEmax);
    h_beta_energy_bg    = (TH1D*)h_beta_energy->Clone("beta_energy_bg");
    h_beta_energy_sub   = (TH1D*)h_beta_energy->Clone("beta_energy_sub");
    h_beta_energy_sub   ->SetTitle("Candidate betas after background subtraction");

    h_alpha_energyVsdT = new TH2D("alpha_energyVsdT","Alpha candidates;#DeltaT [#mus];Energy [MeVee]",dTbins,0.,fdT_max,alphaEbins,0,alphaEmax);
    h_alpha_energyVsdT->SetOption("colz");

    int   timeBins  = binPeriodMax[fConfig]/binPeriodHrs[fConfig];
    float timeMax   = binPeriodMax[fConfig];
    h_time_vs_rate      = new TH1D("time_vs_rate",";Time [hr];Rate per 3.2 ms readout",timeBins,0,timeMax);
    h_time_vs_activity  = (TH1D*)h_time_vs_rate->Clone("time_vs_activity");
    h_time_vs_activity  ->GetYaxis()->SetTitle("Equivalent activity [mBq/kg]");
    h_time_vs_rate_bg   = (TH1D*)h_time_vs_rate->Clone("time_vs_rate_BG");
    h_time_vs_rate_bg   ->SetTitle("Background component");

    if( _isMC ) {
      h_true_alpha_depne      = (TH1D*)h_alpha_charge->Clone("true_alpha_depne");
      h_true_alpha_depne      ->SetTitle("True ionization electrons from alpha");
      h_true_alpha_charge     = (TH1D*)h_alpha_charge->Clone("true_alpha_charge");
      h_true_alpha_charge     ->SetTitle("True alpha charge at anode");
      h_matched_alpha_charge  = (TH1D*)h_alpha_charge->Clone("matched_alpha_charge");
      h_matched_alpha_charge  ->SetTitle("Cluster charge matched to alpha");
    }

    // =====================================================
    // Diagnotic and utility histograms
    tdir_util->cd();
    h_2D_time_vs_dT   = new TH2D("2D_time_vs_dT",";Time [hr];#DeltaT [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_2D_time_vs_dT_bg= (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_bg");
    h_time_vs_N       = new TH1D("time_vs_N",";Time [hr];Number of entries into dT plot",timeBins,0,timeMax);
  }
  
  
  //#################################################################################
  // Primary macro
  //#################################################################################
  void BiPo_macro()
  {

    // ******************************* 
    // Initial configurations
    // *******************************
    infile_t inFile = fInputFiles[fConfig];
    _isMC           = inFile.isMC;
    
    printf("Reading input file: %s : %s\n",inFile.fileName.c_str(),inFile.treeName.c_str());
    if( !_isMC ) chargeScaling = 1.0;

    if( chargeScaling < 1. ) 
      printf("WARNING: scalinig charge by %f\n",chargeScaling);
    // open the file and set up the TTree
    std::string   _fileName = "files/" + inFile.fileName;
    TFile* file = new TFile(_fileName.c_str(),"READ");
    fTree = (TTree*)file->Get(inFile.treeName.c_str());
  
    // set branches
    fTree->SetBranchAddress("lifetime",&lifetime);
    fTree->SetBranchAddress("timestamp",&timestamp);
    fTree->SetBranchAddress("nclusts",&nclusts);                     
    fTree->SetBranchAddress("clust_nwires",&clust_nwires);
    fTree->SetBranchAddress("clust_plane",&clust_plane);              
    fTree->SetBranchAddress("clust_startwire",&clust_startwire);                
    fTree->SetBranchAddress("clust_endwire",&clust_endwire);                
    fTree->SetBranchAddress("clust_nhits",&clust_nhits);              
    fTree->SetBranchAddress("clust_charge",&clust_charge);            
    fTree->SetBranchAddress("clust_time",&clust_time);                
    fTree->SetBranchAddress("clust_blipid",&clust_blipid);            
    fTree->SetBranchAddress("nblips",&nblips);                       
    fTree->SetBranchAddress("blip_nplanes",&blip_nplanes);
    fTree->SetBranchAddress("blip_energy",&blip_energy);
    fTree->SetBranchAddress("blip_y",&blip_y);                        
    fTree->SetBranchAddress("blip_z",&blip_z);                        
    fTree->SetBranchAddress("blip_charge",&blip_charge);              
    fTree->SetBranchAddress("blip_pl0_clustid",&blip_clustid[0]);
    fTree->SetBranchAddress("blip_pl1_clustid",&blip_clustid[1]);
    fTree->SetBranchAddress("blip_pl2_clustid",&blip_clustid[2]);
    if( _isMC ) {
      fTree->SetBranchAddress("nedeps",&nedeps);
      fTree->SetBranchAddress("clust_edepid",&clust_edepid);            
      fTree->SetBranchAddress("nparticles",&nparticles);
      fTree->SetBranchAddress("part_isPrimary",&part_isPrimary);
      fTree->SetBranchAddress("part_mother",&part_mother);
      //fTree->SetBranchAddress("part_trackID",&part_trackID);
      fTree->SetBranchAddress("part_pdg",&part_pdg);
      fTree->SetBranchAddress("part_KE",&part_KE);
      //fTree->SetBranchAddress("part_startT",&part_startT);
      fTree->SetBranchAddress("edep_g4id",&edep_g4id);
      fTree->SetBranchAddress("edep_energy",&edep_energy);
      fTree->SetBranchAddress("edep_pdg",&edep_pdg);
      fTree->SetBranchAddress("edep_electrons",&edep_electrons);
      fTree->SetBranchAddress("edep_charge",&edep_charge);
      fTree->SetBranchAddress("edep_tdrift",&edep_tdrift);
      //fTree->SetBranchAddress("edep_x",&edep_x);
      //fTree->SetBranchAddress("edep_y",&edep_y);
      //fTree->SetBranchAddress("edep_z",&edep_z);
    }

    // make output file to store plots
    std::string _outFileName = "output/plots_bipo_" + inFile.fileName;
    fOutFile = new TFile(_outFileName.c_str(), "recreate");
    tdir_plots  = fOutFile->mkdir("plots");
    tdir_util   = fOutFile->mkdir("util");

    // initialize all histograms
    setRootStyle();
    makeHistograms();
    
    // for picky blip mode, set proper restrictions
    if( fPickyBetaMode ) {
      _betaMinPlanes    = 3;      // Min number of matched planes (must be 2 or 3)
      _betaMaxDiff      = 1;      // Difference in wire intersection points [cm]
    }

    if( fCaloMode ) {
      fBetaCharge_min = 0;
      fBetaEnergy_min = 0;
    }
  
    if( fBackgroundMode == 1 ) _minTick = (int)fdT_max*2;
    _liveTimePerEvt  = (_maxTick-_minTick)*fSamplePeriod*1e-6; //sec

    // Print-interval counter
    int print_counter = 0;

    for(auto iwire : fNoisyWires )  wireIsNoisy[iwire] = true;
    for(auto iwire : fBadWires )    wireIsBad[iwire] = true;




    // ****************************************************
    // Loop over the events
    // ****************************************************
    std::time_t loopStart = time(0);
    size_t totalEntries = fTree->GetEntries();
    for(size_t iEvent=0; iEvent < totalEntries; iEvent++){

      print_counter++;
      if( print_counter > 1000 ) {
        printf("========== EVENT %lu / %lu, %6.2f %%, BiPo count: %i =====================\n",
          iEvent,
          totalEntries,
          100*iEvent/float(totalEntries),
          _numBiPo);
        print_counter = 1;
      }
      
      // ..... quick-test options ...........
      //int maxEvt    = 1000; if(  iEvent >= maxEvt ) break;
      //int sparsify  = 1000; if(  (iEvent % sparsify) != 0 ) continue; 
      //..................................

      // Retrieve event info
      fTree->GetEntry(iEvent);
     

      // Record the event time relative to start of dataset period
      double eventHr = ( timestamp - inFile.t0 ) / 3600.;
      if( !_isMC && (eventHr < fMinHr || eventHr > fMaxHr) ) continue;
      h_time_vs_N->Fill(eventHr);
      _numEvents++;
     

      // ===============================================================
      // Clear masks to track cluster availability to avoid double-counts
      // ===============================================================
      _clustAvailable.assign(nclusts, true);
  
      // ======================================
      // Check truth info
      // ======================================
      int alphaPDG = 1000020040;

      std::vector<bool> g4_isBeta(nparticles,false);
      
      // keep track of the number of primary electrons per decay
      // (we want to exclude events with Auger electrons, which 
      // are unlikely to be selected in our sample)
      int nPrimaryElectrons = 0;
      float sumBetaEnergy = 0;
      for(int i=0; i<nparticles; i++){
        int pdg = part_pdg[i];
        int isPrimary = part_isPrimary[i];
        int mother = part_mother[i];
        float KE = part_KE[i];
        


        if( !part_isPrimary[i] || part_mother[i] > 0 ) continue;
        //printf("  %i   isPrimary? %i   PDG: %i    Mother: %i    KE: %f\n",i,(int)isPrimary,pdg,mother,KE);
        
        if( part_pdg[i]==11 ) {
          g4_isBeta[i] = true;
          nPrimaryElectrons++;
          sumBetaEnergy += KE;
          //h_beta_trueEnergy->Fill(KE);
        }

        if( part_pdg[i] == alphaPDG ) {
          //if( nPrimaryElectrons == 1 ) 
          h_beta_trueEnergySum->Fill(sumBetaEnergy);
          nPrimaryElectrons = 0;
          sumBetaEnergy = 0;
          //std::cout<<"Found alpha\n";
        }

      
      }
      
        for(int i=0; i<nedeps; i++){
          int g4index = edep_g4id[i]; //FindG4Index(edep_g4id[i]);
          int   q_dep   = edep_electrons[i];
          int   q_drift = edep_charge[i];
          if( !part_isPrimary[g4index] ) continue;
          if( edep_pdg[i] == alphaPDG ) {
            if( q_drift > 0 ) h_true_alpha_charge ->Fill(q_drift/0.826);
            if( q_dep > 0 )   h_true_alpha_depne  ->Fill(q_dep);
            float driftTime = edep_tdrift[i];
            float readoutTime = driftTime + part_startT[i];
            int   readoutTick = readoutTime/fSamplePeriod;
            if( readoutTick > _minTick && readoutTick < _maxTick ) 
              _numBiPo_true_perfectReco++;
          }
        }

      
        // look for collection plane clusts matched to an alpha, beta, or gamma
        clust_isAlpha.assign(nclusts, false); 
        clust_isBeta.assign(nclusts, false); 
        clust_isGamma.assign(nclusts, false); 
        for(int i=0; i < nclusts; i++){
          //if( clust_plane[i] != 2 ) continue;
          int eid = clust_edepid[i];
          if( eid < 0 ) {
            if( fIgnoreNonMC ) _clustAvailable[i] = false;
            continue;
          }
          int g4index = edep_g4id[eid]; //FindG4Index(edep_g4id[eid]);
          if( part_isPrimary[g4index] && part_mother[g4index] == 0 ) {
            if(part_pdg[g4index] == alphaPDG )  clust_isAlpha[i] = true;
            if(part_pdg[g4index] == 11)         clust_isBeta[i]  = true;
          } else {
            if(part_pdg[g4index] == 11)         clust_isGamma[i] = true;
          }
          if( clust_isAlpha[i] ) h_matched_alpha_charge->Fill( clust_charge[i] );
          // Option to ignore alpha blips for background assessment
          if( fIgnoreTrueAlphas && clust_isAlpha[i] ) _clustAvailable[i] = false;
          if( fIgnoreTrueGammas && clust_isGamma[i] ) _clustAvailable[i] = false;

          if( clust_plane[i] == 2 ) {
          
            // ---beta plots ---
            if( clust_isBeta[i] ) {
              h_beta_nwires->Fill(clust_nwires[i]);
              h_beta_trueEnergySum_reco->Fill(edep_energy[eid]);
            }
            

            // --- alpha plots ---
            if( clust_isAlpha[i] ) h_alpha_nwires->Fill(clust_nwires[i]);

          }
        }



      // ====================================================
      // Map of clust IDs per wire on collection plane
      // ====================================================
      _map_wire_clusters.clear();
      for(int i=0; i < nclusts; i++){
        if( clust_plane[i] != 2 ) continue;
        clust_charge[i] *= chargeScaling;
        h_wt_clusts->Fill(clust_startwire[i],clust_time[i]);
        for(int j=clust_startwire[i]; j<=clust_endwire[i]; j++)
          _map_wire_clusters[j].push_back(i);
      }

      
      /*
      // ==============================================================
      // Create list of blip IDs sorted by charge
      // ==============================================================
      std::vector<int> sortedBlips;
      std::vector<bool> flag(nblips, false);
      for(int i=0; i<nblips; i++){
        int leadID = -9;
        for(int j=0; j<nblips; j++){
          if( flag[j] ) continue;
          if( leadID < 0 ) leadID = j;
          if( blip_charge[j] > blip_charge[leadID] ) leadID = j;
        }
        sortedBlips.push_back(leadID);
        flag[leadID] = true;
      }
      */
      
      // ==============================================================
      // Loop over 3D blips...
      // ==============================================================
      //for(auto& iBlip : sortedBlips ) {
      for(int iBlip=0; iBlip<nblips; iBlip++){

        // find associated cluster on collection
        int   ic    = blip_clustid[2][iBlip];
        int   ied   = clust_edepid[ic];
  
        // plot wire-time coordinate
        h_wt_blips->Fill( clust_startwire[ic], clust_time[ic] );
  
        // skip if this cluster was already included in a BiPo candidate 
        if( !_clustAvailable[ic] ) continue;
      
        
        //float corrFact    = 1.0;
        float beta_charge = clust_charge[ic]; //*corrFact;
        float beta_energy = blip_energy[iBlip]; //( beta_charge / fRecomb ) * (23.6e-6);
        int   beta_wires  = clust_nwires[ic];

        // blip plane-matching requirements
        if( blip_nplanes[iBlip] < _betaMinPlanes ) continue;
        //if( blip_sigmayz[iBlip] > _betaMaxDiff ) continue;
        
        // evaluate if in fiducial volume
        if( fFidVolCut ) {
          if( blip_z[iBlip] < fZlim[0] || blip_z[iBlip] > fZlim[1] ) continue;
          if( blip_y[iBlip] < fYlim[0] || blip_y[iBlip] > fYlim[1] ) continue;
        }
        
        // apply charge/size cuts on beta
        if(   beta_wires  > fBetaWires_max  ||
              //beta_charge < fBetaCharge_min || 
              //beta_charge > fBetaCharge_max ) 
              beta_energy < fBetaEnergy_min || 
              beta_energy > fBetaEnergy_max ) 
          continue;


  			// skip if any wire in cluster is from a noisy wire
        // skip if this cluster is at the very edge of the wireplane
        int w_start = std::max(0,clust_startwire[ic]-1);
        int w_end   = std::min(clust_endwire[ic]+1,(int)nWiresColl-1);
        if( w_start < 10 || w_end > int(nWiresColl-10) ) continue;
        bool flag = false;
        for(int iwire = w_start; iwire <= w_end; iwire++){
          if( wireIsNoisy[iwire] || wireIsBad[iwire] ) {
            flag = true;
            break;
          }
        }
        if( flag ) continue;
        
        
        
        // require semi-isolation (no other clusters within +/-20us)
        for(int iwire = w_start; iwire <= w_end; iwire++){
          for(auto& jc : _map_wire_clusters[iwire] ) {
            if( ic == jc ) continue;
            float dt = fabs(clust_time[jc]-clust_time[ic]);
            if( dt < 20 ) {flag = true; break;}
          }
        } 
        if( flag ) continue;


        // fill some blip/cluster location histograms
        h_wt_blips_filt->Fill( clust_startwire[ic], clust_time[ic] );
       
        if( clust_isBeta[ic] ) {
          h_beta_trueEnergySum_recoCuts->Fill( edep_energy[ied] );
        }

        // skip if we are near end of wire (account for 400us/800tick trigger offset)
        float peakT = clust_time[ic]+800;
        if( peakT < 0 )         continue;
        if( peakT > _maxTick )  continue;
        if( peakT < _minTick )  continue; 
        
        
        //std::cout<<"beta charge "<<beta_charge<<"\n";
        


        // ------------------------------------------------------------
        // Begin search for alpha candidates
        // ------------------------------------------------------------
        int nclusts_inwindow        = 0;
        int nclusts_20us            = 0;

        //std::cout<<"Standard cand search\n";
        // Search for standard candidates first
        int nwires = 0;
        std::vector<BiPoCandidate> v_cands = FindCandidates(ic, 0, fWireRange, false, nclusts_inwindow,nclusts_20us,nwires);
        h_nclusts_inwindow->Fill(float(nclusts_inwindow)/nwires);
        h_nclusts_20us->Fill(float(nclusts_20us)/nwires);
        h_ncands_inwindow->Fill((int)v_cands.size());
       
        //std::cout<<"FOUND "<<v_cands.size()<<" CANDIDATES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        int ngood = 0;
        for(auto& thisCand : v_cands) {
          _clustAvailable[ic] = false;
          _clustAvailable[thisCand.id2] = false;
          if( thisCand.dT < fdT_min ) continue;
          h_beta_charge         ->Fill(thisCand.q1);
          h_beta_energy         ->Fill(thisCand.e1);
          h_alpha_charge        ->Fill(thisCand.q2);
          h_alpha_energy        ->Fill(thisCand.e2);
          h_alpha_energyVsdT    ->Fill(thisCand.dT,thisCand.e2);
          //if( thisCand.q1 < fBetaCharge_min ) continue;
          if( thisCand.e1 < fBetaEnergy_min ) continue;
          ngood++;
          h_zy_bipos->Fill( blip_z[iBlip], blip_y[iBlip]);
          h_wt_bipos->Fill( clust_startwire[ic], clust_time[ic]);
          h_cand_dT             ->Fill(thisCand.dT);
          h_2D_time_vs_dT       ->Fill(eventHr,thisCand.dT);
          //h_alpha_nwires        ->Fill(clust_nwires[thisCand.id2]);
        }
        if( ngood ) {
            _numBiPo++;
            if( clust_isBeta[ic] ) _numBiPo_mcmatch++;
        }
        
        /*
        // Evaluate standard candidates
        if( v_cands.size() ) {
           
          // flag this beta candidate as unavailable
          _clustAvailable[ic] = false;
          
          // skim just the alphas passing the dT cut
          std::vector<BiPoCandidate> v_cands_dtcut;
          for(auto& thisCand : v_cands) {
            _clustAvailable[thisCand.id2] = false;
            if( thisCand.dT < fdT_min ) continue;
            v_cands_dtcut.push_back( thisCand );
          }
          
          
          if( v_cands_dtcut.size() ) {
            _numBiPo++;
            if( clust_isBeta[ic] ) _numBiPo_mcmatch++;
            // plot locations
            h_zy_bipos->Fill( blip_z[iBlip], blip_y[iBlip]);
            h_wt_bipos->Fill( clust_wire[ic], clust_time[ic]);
            //float weight = 1./v_cands_dtcut.size();
            for(auto& thisCand : v_cands_dtcut) {
              //h_beta_ratio          ->Fill(thisCand.beta_ratio);
              h_beta_charge         ->Fill(thisCand.q1);
              h_beta_energy         ->Fill(thisCand.e1);
              h_alpha_charge        ->Fill(thisCand.q2);
              h_alpha_energy        ->Fill(thisCand.e2);
              h_cand_dT             ->Fill(thisCand.dT);
              h_2D_time_vs_dT       ->Fill(eventHr,thisCand.dT);
              h_alpha_nwires        ->Fill(clust_nwires[thisCand.id2]);
            }
          }
        }//end evaluation of standard cands
        */
       
        //std::cout<<"Control regions\n";
        // Search for candidates in control region shifted +5 wires away
        int a, b, n;
        std::vector<BiPoCandidate> v_control      = FindCandidates(ic, 5, 2, false, a, b, n );
        std::vector<BiPoCandidate> v_control_flip = FindCandidates(ic, 5, 2, true, a, b, n );
        for( auto& vc : v_control )       {
          h_control_dT ->Fill( vc.dT );
          //if( vc.q1 > fBetaCharge_min ) h_control_dT      ->Fill( vc.dT );
        }
        for( auto& vc : v_control_flip )  {
          h_control_dT_flip ->Fill( vc.dT );
          //if( vc.q1 > fBetaCharge_min ) h_control_dT_flip ->Fill( vc.dT );
        }
        
        // Search for background candidates (wire shift)
        //int nwires_shift = 0;
        //int shift = 2*fWireRange+1;
        //int nclusts_inwindow_shift  = 0;
        //int nclusts_20us_shift      = 0;
        //std::vector<BiPoCandidate> v_cands_shift = FindCandidates(ic,shift,fWireRange,false,nclusts_inwindow_shift,nclusts_20us_shift,nwires_shift);

        //std::cout<<"background region\n";
        // Search for background candidates (same wire, but dT-flip)
        int nclusts_inwindow_flip;
        int nclusts_20us_flip;
        int nwires_flip = 0;
        std::vector<BiPoCandidate> v_cands_flip = FindCandidates(ic,0, fWireRange, true, nclusts_inwindow_flip,nclusts_20us_flip,nwires_flip);


       
        /*
        // --------------------------------------------
        // Evaluate wire shift candidates
        // ---------------------------------------------
        if(     v_cands_shift.size() ) {

          std::vector<BiPoCandidate> v_cands_dtcut;
          for(auto& thisCand : v_cands_shift) {
            if( thisCand.dT > fdT_min ) v_cands_dtcut.push_back( thisCand );
          }
          
          if( v_cands_dtcut.size() ) {
            //float weight_n = nwires/float(nwires_shift);
            //float weight = weight_n/v_cands_dtcut.size();
            for(auto& thisCand : v_cands_dtcut) {
              h_cand_dT_shift->Fill(thisCand.dT);
              if( fBackgroundMode == 0 ) {
                h_2D_time_vs_dT_bg       ->Fill(eventHr,thisCand.dT);
                h_beta_charge_bg         ->Fill(thisCand.q1);
                h_beta_energy_bg         ->Fill(thisCand.e1);
                h_alpha_charge_bg        ->Fill(thisCand.q2);
                h_alpha_energy_bg        ->Fill(thisCand.e2);
              }
            }
          }
        }
        */
        
        // --------------------------------------------
        // Evaluate dT-flip candidates
        // ---------------------------------------------
        for(auto& thisCand : v_cands_flip) {
          _clustAvailable[ic] = false;
          _clustAvailable[thisCand.id2] = false;
          if( thisCand.dT < fdT_min ) continue;
          if( fBackgroundMode == 1 ) {
            h_beta_charge_bg         ->Fill(thisCand.q1);
            h_beta_energy_bg         ->Fill(thisCand.e1);
            h_alpha_charge_bg        ->Fill(thisCand.q2);
            h_alpha_energy_bg        ->Fill(thisCand.e2);
            h_zy_bipos_bg->Fill( blip_z[iBlip], blip_y[iBlip]);
          }
          //if( thisCand.q1 < fBetaCharge_min ) continue;
          h_cand_dT_flip             ->Fill(thisCand.dT);
          if( fBackgroundMode == 1 )
            h_2D_time_vs_dT_bg    ->Fill(eventHr,thisCand.dT);
        }

        /*
        if(     v_cands_flip.size() ) {
          std::vector<BiPoCandidate> v_cands_dtcut;
          for(auto& thisCand : v_cands_flip) {
            _clustAvailable[thisCand.id2] = false;
            if( thisCand.dT > fdT_min ) v_cands_dtcut.push_back( thisCand );
          }
          
          if( v_cands_dtcut.size() ) {
            //float weight = 1./v_cands_dtcut.size();
            for(auto& thisCand : v_cands_dtcut) {
              h_cand_dT_flip->Fill(thisCand.dT);
              if( fBackgroundMode == 1 ) {
                h_2D_time_vs_dT_bg       ->Fill(eventHr,thisCand.dT);
                h_beta_charge_bg         ->Fill(thisCand.q1);
                h_alpha_charge_bg        ->Fill(thisCand.q2);
                h_beta_energy_bg         ->Fill(thisCand.e1);
                h_alpha_energy_bg        ->Fill(thisCand.e2);
                //h_beta_ratio_bg          ->Fill(thisCand.beta_ratio);
              }
            }
          }
        }
        */
      

      }//end loop over 3D blips

    }//endloop over events
    double loopDuration = ( time(NULL) - loopStart );


    // Control region ratio: ratio of forward to backward
    

    // ***************************************************
    // Scale everything so it's 'per second in full AV'
    // ***************************************************
    _totalLiveTime = float(_numEvents) * _liveTimePerEvt;
    float scaleFact = (1./_totalLiveTime)*(1./_fiducialFrac);
    h_cand_dT           ->Scale( scaleFact); //, "width");
    h_cand_dT_shift     ->Scale( scaleFact); // "width");
    h_cand_dT_flip      ->Scale( scaleFact); //, "width");
    //h_2D_time_vs_dT     ->Scale( scaleFact );
    //h_2D_time_vs_dT_bg  ->Scale( scaleFact );
    
    std::cout<<"Total live time: "<<_totalLiveTime<<"\n";
    std::cout<<"Fiducial fraction: "<<_fiducialFrac<<"\n";
    std::cout<<"Scale factor: "<<scaleFact<<"\n";
    
    h_alpha_charge      ->Scale( scaleFact );
    h_alpha_charge_bg   ->Scale( scaleFact );
    h_alpha_energy      ->Scale( scaleFact );
    h_alpha_energy_bg   ->Scale( scaleFact );
    h_beta_charge       ->Scale( scaleFact );
    h_beta_charge_bg    ->Scale( scaleFact );
    h_beta_energy       ->Scale( scaleFact );
    h_beta_energy_bg    ->Scale( scaleFact );
    
    // Keep sum squares
    h_control_dT        ->Sumw2();
    h_control_dT_flip   ->Sumw2();
    h_control_dT_ratio  ->Divide(h_control_dT,h_control_dT_flip);
    h_control_dT_ratio  ->Fit(f_backward_corr,"R");
    if( fLinearizeCorr )  h_cand_dT_flip->Multiply( f_backward_corr );
    else                  h_cand_dT_flip->Multiply( h_control_dT_ratio );
   
    // ***************************************************
    // Histogram subtraction time!
    // ***************************************************
    h_cand_dT_sub         ->Add(h_cand_dT,        1);
    h_alpha_charge_sub    ->Add(h_alpha_charge,   1);
    h_beta_charge_sub     ->Add(h_beta_charge,    1);
    h_alpha_energy_sub    ->Add(h_alpha_energy,   1);
    h_beta_energy_sub     ->Add(h_beta_energy,    1);
    h_zy_bipos_sub         ->Add(h_zy_bipos,      1);
    //h_beta_ratio_sub      ->Add(h_beta_ratio,     1);

    switch(fBackgroundMode){
      case 0: h_cand_dT_sub->Add(h_cand_dT_shift,-1.); break;
      case 1: h_cand_dT_sub->Add(h_cand_dT_flip,-1.); break;
    } 

    h_alpha_charge_sub->Add(h_alpha_charge_bg,-1.);
    h_beta_charge_sub->Add(h_beta_charge_bg,-1.);
    h_alpha_energy_sub->Add(h_alpha_energy_bg,-1.);
    h_beta_energy_sub->Add(h_beta_energy_bg,-1.);
    h_zy_bipos_sub->Add(h_zy_bipos_bg,-1);
    //h_beta_ratio_sub->Add(h_beta_ratio_bg,-1.);

    // ***************************************************
    // Write all histos currently in stack
    // ***************************************************
    fOutFile->Write(); 
    makePlots();
    
    printf("\n*******************************************\n");
    printf("File                : %s\n",      inFile.fileName.c_str()); 
    printf("Total events        : %i\n",      _numEvents); 
    printf("Total live time     : %f sec\n",  _totalLiveTime);
    printf("Live time per evt   : %f us\n",   _liveTimePerEvt*1e6);
    printf("Fiducial fraction   : %f\n",      _fiducialFrac); 
    printf("dT min/max          : %.2f-%.2f us\n", fdT_min,fdT_max);
    printf("Ave cands / evt     : %f\n",h_cand_dT->GetEntries()/(float)_numEvents );
    printf("  - truth-matched   : %f\n",_numBiPo_mcmatch/(float)_numEvents );
    printf("BG subtrctn method  : %i\n",fBackgroundMode);
    printf("Processing time     : %f sec (%f sec/evt)\n", loopDuration, loopDuration/float(_numEvents));
    printf("Excluded %i noisy wires \n",     (int)fNoisyWires.size());	
    printf("*******************************************\n\n");
    fOutFile->Close();
   
    std::cout<<"Num bipos if perfect reco: "<<_numBiPo_true_perfectReco<<"\n";
    std::cout<<_numBiPo_true_perfectReco/_totalLiveTime<<"\n";
  }
  
  
  
  //#################################################################################
  // Make plots here
  //#################################################################################
  void makePlots()
  {
    
    std::cout<<"MAKING PLOTS "<<_isMC<<"\n";
    // Histograms needed:
    //   - h_cand_dT_sub
    //   - h_2D_time_vs_dT
    //   - h_2D_time_vs_dT_bg
    //   - h_time_vs_N
    //   - h_alpha_charge_sub
    //   - h_beta_charge_sub
  

    std::string name;
    float       range;
    float       min;
    float       max;
   
    if( !_isMC ) {

      // ============================================
      // Do slice-by-slice dT fit
      // ============================================
      TH1D* h_time_vs_p0 = (TH1D*)h_time_vs_rate->Clone("time_vs_p0");
      TH1D* h_time_vs_p1 = (TH1D*)h_time_vs_rate->Clone("time_vs_p1");

      TH1D* h_slice;
      TH1D* h_bg;
      int nbins = h_time_vs_N->GetXaxis()->GetNbins();
      for(int i=1; i<=nbins; i++){
        h_slice     = Make1DSlice( h_2D_time_vs_dT, i, i, Form("time_vs_dT_%i",i) );
        h_bg        = Make1DSlice( h_2D_time_vs_dT_bg, i, i, Form("time_vs_dT_flip_%i",i) );
        double liveTime   = h_time_vs_N->GetBinContent(i)*_liveTimePerEvt;
        float scaleFact = (1./liveTime)*(1./_fiducialFrac);
        h_slice ->Scale(scaleFact);//, "width");
        h_bg    ->Scale(scaleFact);//, "width");
        if( fLinearizeCorr )  h_bg->Multiply( f_backward_corr );
        else                  h_bg->Multiply( h_control_dT_ratio );
        h_slice ->Add( h_bg, -1. );
        
        //tdir_util->cd();
        FitResult fr = fitdT(h_slice,true,false);
        if( fr.rate_signal != -9 ) {
          h_time_vs_activity  ->SetBinContent(  i,  fr.activity);
          h_time_vs_activity  ->SetBinError(    i,  fr.activity_err);
          h_time_vs_rate      ->SetBinContent(  i,  fr.rate_signal);
          h_time_vs_rate      ->SetBinError(    i,  fr.rate_signal_err);
          h_time_vs_rate_bg   ->SetBinContent(  i,  fr.N_bg);
          h_time_vs_rate_bg   ->SetBinError(    i,  fr.N_bg_err);
          h_time_vs_p0        ->SetBinContent(  i,  fr.p0);
          h_time_vs_p0        ->SetBinError(    i,  fr.p0_err);
          h_time_vs_p1        ->SetBinContent(  i,  fr.p1);
          h_time_vs_p1        ->SetBinError(    i,  fr.p1_err);
        }
      }
      
      fOutFile->cd();
      h_time_vs_activity  ->Write(0, TObject::kOverwrite );
      h_time_vs_rate      ->Write(0, TObject::kOverwrite );
      h_time_vs_rate_bg   ->Write(0, TObject::kOverwrite );
      
      std::cout<<"Rate at 2 hours: "<<h_time_vs_rate->Interpolate(2)<<"\n";
      std::cout<<"Rate at 4 hours: "<<h_time_vs_rate->Interpolate(4)<<"\n";
      std::cout<<"Rate at 6 hours: "<<h_time_vs_rate->Interpolate(6)<<"\n";


      // ============================================
      // Plot fit parameters vs time
      // ============================================
      TH1D* h_p0  = h_time_vs_p0; //h_time_vs_rate_bg; //h_time_vs_p0;
      TH1D* h_p1  = h_time_vs_p1; //h_time_vs_rate; //h_time_vs_p1;
      FormatTH1D(h_p1, kBlue, 1, 2, 20, 1);
      FormatTH1D(h_p0, kRed, 1, 2, 20, 1);
      name = "c_time_vs_params";
      min = std::min(GetHistMin(h_p0),GetHistMin(h_p1));
      max = std::max(GetHistMax(h_p0),GetHistMax(h_p1));
      range = (max-min);
      TCanvas* c = new TCanvas(name.c_str(),name.c_str(),500,380);
      gStyle->SetOptStat(0);
      gPad->SetGridy(1);
      h_p0->GetYaxis()->SetRangeUser(min-0.3*range, max+0.5*range);
      h_p0->GetYaxis()->SetTitleOffset(1.1); 
      //h_p0->GetYaxis()->SetTitle("#DeltaT fit parameter value");
      h_p0->GetYaxis()->SetTitle("Fit component integral");
      h_p0->DrawCopy();
      h_p1->DrawCopy("same");
      
      //tdir_plots->cd();
      c->Write();


      // ============================================
      // Plot BiPo rate vs time
      // ============================================
      name = "c_time_vs_rate";
      TCanvas* c2 = new TCanvas(name.c_str(),name.c_str(),550,400);
      gStyle->SetOptStat(0);
      gPad->SetGridy(1);
      TH1D* h1  = h_time_vs_rate;
      FormatTH1D(h1, kBlue, 1, 1, 20, 0.7);
      FormatAxes(h1, 0.05, 0.045, 1.1, 1.1);
      h1->GetYaxis()->SetTitleOffset(1.1); 
      h1->GetYaxis()->SetTitle("Candidates per 3.2 ms readout");
      h1->DrawCopy();
      c2->Write();
      
      // ============================================
      // Plot activity vs time
      // ============================================
      name = "c_time_vs_activity";
      TCanvas* c3 = new TCanvas(name.c_str(),name.c_str(),550,400);
      gStyle->SetOptStat(0);
      gPad->SetGridy(1);
      TH1D* h2 = h_time_vs_activity;
      FormatTH1D(h2, kBlack, 1, 1, 20, 0.7);
      FormatAxes(h2, 0.05, 0.045, 1.1, 1.1);
      h2->GetYaxis()->SetTitle("Equivalent activity [mBq/kg]");
      h2->DrawCopy();
      c3->Write();
    
    }


    // ============================================
    // Do final fit on dT spectrum
    // ============================================
    TH1D* h3  = h_cand_dT_sub;
    h3->GetXaxis()->SetTitle("#DeltaT [#mus]");
    h3->GetYaxis()->SetTitle("Candidates per second / 20 #mus");
    FormatTH1D(h3, kBlue+2, 1, 1, 21, 0.5);
    FormatAxes(h3, 0.05, 0.045, 1.1, 1.1);
    //fitdT( h3, false, true );
    tdir_plots->cd();
    fitdT( h3, true, false );

  } 
  
  
  
  //################################################################################
  // Function that performs the dT fit
  //#################################################################################
  FitResult fitdT(TH1D* h, bool writeCanvas = false, bool constrainNorm = false ){
    
    std::cout
    <<"\n\nFitting dT spectrum "<<h->GetTitle()<<", "<<h->GetEntries()<<"\n";
  
    FitResult out;
  
    std::string label = h->GetName();
    TCanvas* c = new TCanvas(Form("c_fit_%s",label.c_str()),Form("c_fit_%s",label.c_str()),500,420);
    //gPad->SetMargin(mar_l, mar_r, mar_b, mar_t ); 
    float histMax = GetHistMax(h);
    float histMin = GetHistMin(h);
    float range   = (histMax-histMin);
    
    TGraphErrors* gr = MakeGraph(h);
    
    /*
    // Define fit function, initialize parameters
    TF1* fit = new TF1("FullFit","[0] + [1]*exp(-x/[2]) + [3]*exp(-x/[4])",fdT_min,fdT_max);
    fit->SetParameter(0, histMax/5 );
    fit->SetParLimits(0, 0, histMax );
    fit->SetParameter(3, histMax/5 );
    fit->SetParLimits(3, 0, histMax );
    fit->SetParameter(4, 15 );
    fit->SetParLimits(4, 5, 40);
    fit->SetParameter(1, histMax );
    fit->SetParLimits(1, 0, histMax*2 );
    fit->FixParameter(2, 164.5 );
    */

    // Define single exp fit + flat BG
    TF1* fit = new TF1("FullFit","[0] + [1]*exp(-x/[2])",fdT_min,fdT_max);
    fit->SetParameter(0, (histMax+histMin)/2 );
    fit->SetParameter(1, (histMax-histMin) );
    if( constrainNorm ) {
      fit->SetParLimits(1, 0., 100.);
    }
    fit->FixParameter(2, 164.3 );

    // Draw plot and fit
    c->cd();
    gr->Fit(fit,"RE");
  
    /*
    bool isAtLimit = gMinuit->fLimset;
    std::cout<<"Is MINUIT reaching limit? "<<gMinuit->fLimset<<"\n";

    if( isAtLimit && constrainNorm ) {
      float p1      = fit->GetParameter(1);
      float p1_err  = fit->GetParError(1);
      std::cout<<"Par1: "<<p1<<" +/- "<<p1_err<<"\n";
      std::cout<<"Histogram range: "<<range<<"\n";
      if( p1_err > range ) {
        std::cout<<"!!! error is beyond range of histogram, so probably bad !!!\n";
        p1_err = 0;
      }
      TF1* fit2 = (TF1*)fit->Clone("FullFit2");
      fit2->FixParameter(0,fit->GetParameter(0)-fit->GetParError(0));
      hc->Fit(fit2,"REN");
      double newp1 = fit2->GetParameter(1);
      double diff = fabs(p1 - newp1);
      double err = sqrt( pow(p1_err,2) + pow(diff,2) );
      std::cout<<"new norm: "<<newp1<<" +/- "<<fit2->GetParError(1) <<"\n"; 
      std::cout<<"Re-assigning fit error: "<<err<<"\n";
      fit->SetParError(1,err);
     
      //std::cout<<"new par1: "<<fit->GetParameter(1)<<" +/- "<<fit->GetParError(1)<<"\n";
      //TF1* fit2 = (TF1*)fit->Clone("FullFit2");
      //fit2->ReleaseParameter(1);
      //hc->Fit(fit2,;
      //double norm1 = fit->GetParameter(1);
      //double norm2 = fit2->GetParameter(1);
      //double diff = fabs(norm2-norm1);
      //double err = sqrt( pow( diff, 2) + pow( fit2->GetParError(1),2 ) );
      //std::cout<<"new norm: "<<norm2<<" +/- "<<fit2->GetParError(1) <<"\n"; 
      //fit->SetParError(1,err);
    }
    */


    //hc->DrawCopy();
    gr->Draw("AP");
    gStyle->SetOptFit(1112);
    //gPad->Modified(); gPad->Update();
    if( writeCanvas ) {
      //hc->Write(0, TObject::kOverwrite );
      tdir_util->cd();
      gr->Write(label.c_str());
      tdir_plots->cd();
      c->Write();
      fOutFile->cd();
    }
  
    // Factor out the signal and BG components into
    // separate TF1 functions
    TF1* f_bipo = new TF1("bipo","[0]*exp(-x/[1])");
    f_bipo->FixParameter(0, fit->GetParameter(1) );
    f_bipo->FixParameter(1, fit->GetParameter(2) );
    TF1* f_flat = new TF1("flat","[0]");
    f_flat->FixParameter(0, fit->GetParameter(0));
  
    // Full extrapolated BiPo component; integral of A*exp(-x/B) from 0-infinity is A*B
    double p0       = fit->GetParameter(0);
    double p0_err   = fit->GetParError(0);
    double p1       = fit->GetParameter(1);
    double p1_err   = fit->GetParError(1);
    double tau      = fit->GetParameter(2);
    double chi2     = fit->GetChisquare();
    int    ndf      = fit->GetNDF();
    double n_bipo   = (p1*tau)/fdT_binSize;
    double n_bipo_err = (p1_err*tau)/fdT_binSize;
    
    // Components within the dT selection window
    double N_total    = h ->Integral(1,h->GetXaxis()->GetNbins());
    double N_bipo     = f_bipo->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_bipo_err = f_bipo->IntegralError(fdT_min,fdT_max)/fdT_binSize;
    double N_bg       = f_flat->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_bg_err   = f_flat->IntegralError(fdT_min,fdT_max)/fdT_binSize;
    
    // As extra systematic, fit the flat comp to zero
    std::cout<<"Re-fitting with flat component set to 0:\n";
    fit->FixParameter(0,0);
    gr->Fit(fit,"QN");
    double n2 = (fit->GetParameter(1)*fit->GetParameter(2))/fdT_binSize;
    double n_bipo_syst = fabs(n_bipo-n2);
    std::cout<<"Nom rate: "<<n_bipo<<"\n";
    std::cout<<"Syst diff: "<<n_bipo_syst<<"\n";
   
    // efficiency factor + error
    float effMC        = (fEfficiencyMC>0) ? fEfficiencyMC : 1.;
    float effMC_syst    = (fEfficiencyMC>0) ? std::max(0.,fEfficiencyMC_err) : 0.;
  
    // Error terms for specific activity calculation
    double err1   = (n_bipo!=0) ? fabs(n_bipo_err/n_bipo) : 0.;
    double err2a  = (n_bipo!=0) ? fabs(n_bipo_syst/n_bipo) : 0.;
    double err2b  = (effMC>0) ? effMC_syst/effMC : 0.;

    double activity       = 1e3 * n_bipo / effMC / 85000.; //calcActivity(n_bipo,effMC);
    double activity_err   = fabs(activity) * err1;
    double activity_syst  = fabs(activity) * sqrt( pow(err2a,2) + pow(err2b,2) );
    double activity_errTot = sqrt( pow(activity_err,2) + pow(activity_syst,2) );
    
    // normalize "rate" to be per 3.2ms readout
    n_bipo      *= 0.0032;
    n_bipo_err  *= 0.0032;
    n_bipo_syst *= 0.0032;
    N_bg        *= 0.0032;
    N_bg_err    *= 0.0032;
    
    printf("================ dT fit =================\n");
    printf("p0          : %f +/- %f\n", p0, p0_err);
    printf("p1          : %f +/- %f\n", p1, p1_err);
    printf("Chi2/ndf    : %f / %i = %f\n",chi2,ndf,chi2/ndf);
    printf("Total rate  : %f per sec\n",N_total);
    printf(" - BiPo     : %f per sec\n",N_bipo);
    printf(" - Flat     : %f per sec\n",N_bg);
    printf("Rate        : %f +/- %f (stat) +/- %f (syst) BiPos per 3.2 ms\n",n_bipo,n_bipo_err,n_bipo_syst);
    printf("            = %f +/- %f BiPos per 3.2 ms\n",n_bipo,sqrt(pow(n_bipo_err,2)+pow(n_bipo_syst,2)));
    printf("Activity    : %f +/- %f (stat) +/- %f (syst) mBq/kg\n",activity,activity_err,activity_syst);
    printf("            = %f +/- %f mBq/kg\n",activity,activity_errTot);
    printf("(assumed eff = %f)\n",effMC);

    out.rate_signal       = n_bipo;
    out.rate_signal_err   = sqrt(pow(n_bipo_err,2) + pow(n_bipo_syst,2) );
    out.N_signal          = N_bipo;
    out.N_signal_err      = N_bipo_err;
    out.N_bg              = N_bg;
    out.N_bg_err          = N_bg_err;
    out.p0                = p0;
    out.p0_err            = p0_err;
    out.p1                = p1;
    out.p1_err            = p1_err;
    out.activity          = activity;
    out.activity_err      = activity_errTot;
    out.activity_err_stat = activity_err;
    out.activity_err_syst = activity_syst;
    return out;
  
  }
  

  //################################################################################
  // Function to search for alpha candidates relative to a beta candidate
  //#################################################################################
  std::vector<BiPoCandidate> FindCandidates(int ic, int wire_shift, int range, bool flipdT, int& npileup, int& npileup_20us, int& nwires ) {
   
    //std::cout<<"looking for cands; beta clust "<<ic<<", beta range "<<clust_startwire[ic]<<", "<<clust_endwire[ic]<<"... shift "<<wire_shift<<"    range "<<range<<"\n";
    std::vector<BiPoCandidate> v;
    std::vector<int> cand_IDs;
    npileup = 0;
    npileup_20us = 0;

    std::set<int> wires;
    int w1  = clust_startwire[ic] - wire_shift;
    int w2  = clust_endwire[ic] + wire_shift;
    
    //std::cout<<"w1 w2 "<<w1<<"  "<<w2<<"\n";
     

    for(int i=w1-range; i<=w1+range; i++) wires.insert(i);
    for(int i=w2-range; i<=w2+range; i++) wires.insert(i);

    nwires = 0;

    for(auto iWire : wires ) {
      //std::cout<<"checking wire "<<iWire<<"\n";

      if( wireIsBad[iWire] || wireIsNoisy[iWire] ) continue;
      nwires++;
      

      //std::cout<<"There are "<<_map_wire_clusters[iWire].size()<<" clusts on this wire\n";
      for(auto& jc : _map_wire_clusters[iWire] ) {
        
        if( ic == jc ) continue;
        if( !_clustAvailable[jc] ) continue;
        if( std::count(cand_IDs.begin(),cand_IDs.end(),jc) ) continue;

        float dT        = (clust_time[jc]-clust_time[ic])*fSamplePeriod;
        if (flipdT) dT  *= -1.;
        
        //std::cout<<"  dT "<<dT<<"\n";
       
        

        // Count clusters in forward window, and fill some pre-cut histos
        if( dT > 0 && dT < fdT_max ) {
          npileup++;
          if( dT < 20 ) npileup_20us++;
          // --- alpha charge/nhits cut ---
          //std::cout<<"clust charge "<<clust_charge[ic]<<"\n";
          float beta_E = (23.6e-6)*(clust_charge[ic]/0.584);
          float alpha_E = (23.6e-6)*(clust_charge[jc]/0.584);
          //std::cout<<"dT = "<<dT<<"    beta "<<beta_E<<"    alpha "<<alpha_E<<"\n";
          if( clust_nwires[jc] <= fAlphaWires_max
            //&&  clust_charge[jc] > fAlphaCharge_min 
            //&&  clust_charge[jc] < fAlphaCharge_max ) {
            &&  alpha_E > fAlphaEnergy_min 
            &&  alpha_E < fAlphaEnergy_max ) {
              //std::cout<<"Found candidate\n";
              //BiPoCandidate c = { clust_blipid[ic], ic, jc, dT, clust_charge[ic], clust_charge[jc], clust_rms[ic]/clust_amp[ic]};
              BiPoCandidate c = { clust_blipid[ic], ic, jc, dT, clust_charge[ic], clust_charge[jc], beta_E, alpha_E};
              v.push_back(c);
              cand_IDs.push_back(jc);
            }
          
         
        }//end min dT req 
      }//<-- end loop over clusters on this wire
    }//<-- endloop over wires
  
    return v;
  
  }


  //###################################################################
  //int FindG4Index(int g4id) {
  //  for(size_t i=0; i<nparticles; i++) {
  //    if( part_trackID[i] == g4id ) return i;
  //  }
  //  return -9;
  //}


