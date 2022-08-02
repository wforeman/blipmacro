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
  int   fConfig   = 1;
  float fMinHr    = 0;
  float fMaxHr    = 999;
  
  // --- Input files ---
  infile_t fInputFiles[4] = {
  /*0*/  { "BlipAna_20220802_BiPo_MC_Overlay_DefaultQY.root", "blipanaTrkMask/anatree",true,0,0},
  ///*0*/  { "BlipAna_20220802_BiPo_MC_Overlay.root", "blipanaTrkMask/anatree",true,0,0},
  //*0*/  { "BlipAna_20220728_BiPo_MC.root", "blipana/anatree",true,0,0},
  /*1*/  { "BlipAna_20220728_RadonData_Phase1.root","blipanaTrkMask/anatree",false,1627415210,1627592728},
  /*2*/  { "BlipAna_20220731_RadonData_Phase2.root","blipanaTrkMask/anatree",false,1627594380,1627761265},
  /*3*/  { "BlipAna_20220731_Data_Run3.root","blipanaTrkMask/anatree",false,1528520000,1531760000}
                                                                          // prev: 1627594369 1627594412
                                                                          // 7/31: 1627592820
                                                                          //       Mike Z.'s elog entry 126161: 
                                                                          //       bypass started 7/29/21 16:33
                                                                          //       --> 1627594380
  };

  double fEfficiencyMC = 0.104;
  double fEfficiencyMC_err = 0.028;

  const float binPeriodHrs[4] = { 2,   2,  2,  48 };
  const float binPeriodMax[4] = { 44, 44, 44,  864 };
  // Run 3 dataset: 
  // T0 = 1528520000 (9 June 2018, 23:53)
  // T1 = 1531760000 (16 July 2018, 11:53)
  //const int timeBins = 22;
  //float timeMax = 44;

  // --- General macro parameters ---
  int   fBackgroundMode   = 1;      // 0= wire-shift, 1= dT-flip
  bool  fPickyBlipMode    = false;  // Require blips match on all 3 planes
  int   fWireRange        = 1;      // +/- range to look for alpha candidate
  float fBetaCharge_min   = 5e3; //3.5e3;  // Min charge of beta candidate blip [e-]
  float fBetaCharge_max   = 100e3;   // Max charge of beta candidate blip [e-]
  float fAlphaCharge_min  = 0e3;    // Min charge of alpha candidate cluster [e-]
  float fAlphaCharge_max  = 5e3;    // Max charge of alpha candidate cluster [e-]
  int   fAlphaHits_max    = 999;      // Max hits in alpha candidate cluster
  int   fAlphaWires_max   = 2;      // Max wire extent on coll plane
  bool  fAlphaReq3D       = false;  // Require alpha pulse be 3D-matched
  float fdT_binSize       = 20.;    // Bin width for all dT spectra plots [us]
  float fdT_min           = 20.;    // Min dT for looking for candidate [us]
  float fdT_max           = 500.;   // Max dT for looking for candidate [us]
  int   fMaxClustMult     = 999;    // Max number of clusters in time window
  int   fMaxCands         = 999;
  float fSphereRadius     = 25;    
  bool  fFidVolCut        = false;
  float fZmin             = 50; //50;     // Z range (0 to 1037 cm)
  float fZmax             = 985; //1035; //985;    //
  float fYmin             = -80; //-120; //-70;    // Y range (-120 to 120 cm)
  float fYmax             = 80; //120;     //

  // --- MC options ---
  bool  fIgnoreTrueAlphas = false;
  bool  fIgnoreTrueGammas = false;

  // --- Detector properties ---
  int   nWiresColl        = 3455;
  float fSamplePeriod     = 0.5; // microseconds
  
  // --- Special switches ---
  int   fRandomWireShift  = 0; 
  
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
  float dz = fZmax-fZmin;
  float dy = fYmax-fYmin;
  //float _fidCorFactor = (fFidVolCut) ? std::max(1., (1037.*230.) / (dz*dy) ) : 1.0;
  float _fiducialFrac = (fFidVolCut) ? std::max(1., (dz*dy)/(1037.*230.) ) : 1.0;
  float _effMC        = (fEfficiencyMC>0) ? fEfficiencyMC : 1.;
  float _effMC_err    = (fEfficiencyMC>0) ? std::max(0.,fEfficiencyMC_err) : 0.;
  // Counters / maps / etc
  bool  _isMC           = false;
  int   _numEvents      = 0;
  int   _numBiPo        = 0;
  int   _numBiPo_6_312  = 0;
  int   _numBiPo_true   = 0;
  int   _numBiPo_6_312_true = 0;
  float _totalLiveTime  = 0;
  int   _nwires         = 3456;
  //std::vector<bool> _blipAvailable;
  std::vector<bool> _clustAvailable;
  std::map<int,std::vector<int>> _map_wire_clusters;
  std::vector<bool> wireIsNoisy (_nwires,false);
  std::vector<bool> wireIsBad   (_nwires,false);
  // Live time
  int   _minTick         = 0;
  int   _maxTick         = 6400 - (int)fdT_max*2;
  float _liveTimePerEvt  = _maxTick*fSamplePeriod*1e-6; //sec
  // Blip fidelity requirements
  int _betaMinPlanes = 2;
  int _betaMaxDiff = 9999;

  int _numBiPo_true_perfectReco = 0;

  //##########################################################################
  // Functions and ROOT objects
  //##########################################################################
  void                        makePlots();
  void                        makeHistograms();
  void                        setRootStyle();
  std::vector<BiPoCandidate>  FindCandidates(int, int, bool, int&, int&);
  FitResult                   fitdT(TH1D*,bool,bool);
  int                         FindG4Index(int);
  
  // ROOT objects
  TTree*      fTree;
  TFile*      fOutFile;

  // Histograms
  TDirectory* tdir_util;
  TDirectory* tdir_plots;
  TH1D* h_blip_charge;
  TH1D* h_blip_nhits;
  TH1D* h_alpha_nwires;
  TH1D* h_nclusts_wire_ave;
  TH1D* h_nclusts_inwindow;
  TH1D* h_ncands_inwindow;
  TH1D* h_cand_dT;
  TH1D* h_cand_dT_shift;
  TH1D* h_cand_dT_flip;
  TH1D* h_cand_dT_sub;
  //TH1D* h_clust_dT;
  TH1D* h_hit_gof_2D;
  TH1D* h_hit_gof_3D;
  TH2D* h_zy_blips;
  TH2D* h_zy_blips_filt;
  TH2D* h_zy_blips_picky;
  TH2D* h_zy_bipos;
  TH2D* h_wt_clusts;
  TH2D* h_wt_blips;
  TH2D* h_wt_blips_filt;
  TH2D* h_wt_bipos;
  TH1D* h_blip_sphere_mult;



  /*
  TH1D* h_time_vs_lifetime;
  TH1D* h_time_vs_lifetime_N;
  */

  TH1D* h_time_vs_N;        // basic event count
  
  TH2D* h_2D_time_vs_dT;       // dT for candidate region
  //TH2D* h_2D_time_vs_dT_shift;
  //TH2D* h_2D_time_vs_dT_flip;
  TH2D* h_2D_time_vs_dT_bg;
  //TH2D* h_2D_time_vs_dT_sub;

  TH1D* h_time_vs_rate_bipo; 
  TH1D* h_time_vs_activity; 
  TH1D* h_time_vs_rate_bg;
  TH1D* h_time_vs_ratio;
 
  TH1D* h_timeFine_vs_evts;
  TH1D* h_timeFine_vs_rate;
  
  TH2D* h_time_vs_ntrks;
  
  TH1D* h_alpha_charge;
  TH1D* h_alpha_charge_bg;
  TH1D* h_alpha_charge_sub;

  TH1D* h_beta_charge;
  TH1D* h_beta_charge_bg;
  TH1D* h_beta_charge_sub;

  /* 
  TH2D* h_poszy[timeBins];
  TH2D* h_poszy_shift[timeBins];
  TH2D* h_poszy_sub[timeBins];
  */

  float mar_l = 0.13;
  float mar_r = 0.10;
  float mar_b = 0.12;
  float mar_t = 0.10;

  // MC-truth histograms
  TH1D* h_true_alpha_depne;
  TH1D* h_true_alpha_charge;
  TH1D* h_matched_alpha_charge;
  TH1D* h_true_attenuation;

  //##########################################################################
  // Set default ROOT plot style
  //##########################################################################
  void setRootStyle() {
    
    // set margin sizes
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    // set title offsets (for axis label)
    gStyle->SetTitleXOffset(1.2);
    gStyle->SetTitleYOffset(1.3);

    // use large fonts
    //Int_t font=42; // Helvetica
    Double_t tsize=0.040;
    Double_t lsize=0.040;
    gStyle->SetTitleSize(tsize,"x");
    gStyle->SetLabelSize(lsize,"x");
    gStyle->SetTitleSize(tsize,"y");
    gStyle->SetLabelSize(lsize,"y");
    gStyle->SetTitleSize(tsize,"z");
    gStyle->SetLabelSize(lsize,"z");
    gStyle->SetTextSize(tsize);
  
    // use bold lines and markers
    //gStyle->SetLineWidth(2);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(0.8);
    gStyle->SetMarkerColor(kAzure-6);

    // get rid of X error bars 
    gStyle->SetErrorX(0.001);

    // do not display any of the standard histogram decorations
    //atlasStyle->SetOptTitle(0);
    //atlasStyle->SetOptStat(1111);
    //atlasStyle->SetOptStat(0);
    //atlasStyle->SetOptFit(1111);
    //atlasStyle->SetOptFit(0);

    // put tick marks on top and RHS of plots
    //atlasStyle->SetPadTickX(1);
    //atlasStyle->SetPadTickY(1);
  }

  //##########################################################################
  // Initialize histograms
  //##########################################################################
  void makeHistograms()
  {
    // =====================================================
    // Important final histograms
    fOutFile->cd();
    
    h_blip_charge       = new TH1D("blip_charge","3D Blips;Collection Plane Charge [e];Events", 500,0,50e3);
    h_blip_nhits        = new TH1D("blip_nhits","Candidate blips;Collection Plane Hits",20,0,20);
    h_alpha_nwires       = new TH1D("alpha_nwires","Candidate alpha clusters (coll plane);Wire extent",20,0,20);
    h_nclusts_inwindow  = new TH1D("nclusts_inwindow","Number of clusters in time window following Bi-candidate",20,0,20);
    h_ncands_inwindow   = new TH1D("ncands_inwindow","Number of Po candidates in time window following Bi-candidate",10,0,10);
    h_blip_sphere_mult  = new TH1D("blip_sphere_mult",Form("Sphere radius %f cm;3D blip multiplicity",fSphereRadius),50,0,50);
    
    float Zmin = -100;  float Zmax = 1100;  int Zbins = 600;
    float Ymin = -150;  float Ymax = 150;   int Ybins = 150;
    float Tmin = -1000;     float Tmax = 6000;  int Tbins = 700; 
    float Wmin = -100;  float Wmax = 3500;  int Wbins = 1800; 
    h_zy_blips      = new TH2D("zy_blips","3D blips;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_blips_filt = new TH2D("zy_blips_filt","3D blips (quality cuts);Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_blips_picky = new TH2D("zy_blips_picky","3D blips (PICKY);Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_bipos      = new TH2D("zy_bipos","BiPo candidates;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_blips      ->SetOption("colz");
    h_zy_blips_filt ->SetOption("colz");
    h_zy_blips_picky->SetOption("colz");
    h_zy_bipos      ->SetOption("colz");
    h_wt_clusts     = new TH2D("wt_clusts","2D clusts;Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_blips      = new TH2D("wt_blips","3D blips;Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_blips_filt = new TH2D("wt_blips_filt","3D blips (quality cuts);Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_bipos      = new TH2D("wt_bipos","BiPo candidates;Collection Plane Wire;Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
    h_wt_clusts     ->SetOption("colz");
    h_wt_blips      ->SetOption("colz");
    h_wt_blips_filt ->SetOption("colz");
    h_wt_bipos      ->SetOption("colz");
   
    int dTbins = (fdT_max) / fdT_binSize;
    //h_cand_dT       = new TH1D("cand_dT","Selected BiPo Candidates;Time difference [#mus];Candidates per event", dTbins,0.,fdT_max);
    h_cand_dT       = new TH1D("cand_dT","Selected BiPo Candidates;Time difference [#mus];Candidates per second", dTbins,0.,fdT_max);
    h_cand_dT_shift = (TH1D*)h_cand_dT->Clone("cand_dT_shift"); h_cand_dT_shift ->SetTitle("Shifted-wire region candidates");
    h_cand_dT_flip  = (TH1D*)h_cand_dT->Clone("cand_dT_flip");  h_cand_dT_flip  ->SetTitle("Opposite dT candidates");
    h_cand_dT_sub   = (TH1D*)h_cand_dT->Clone("cand_dT_sub");   h_cand_dT_sub   ->SetTitle("Background-subtracted spectrum");
    
    float alphaQmax = fAlphaCharge_max;
    int   alphaQbins = alphaQmax/200.;
    h_alpha_charge      = new TH1D("alpha_charge","Candidate alphas;Collection Plane Charge [e];Events", alphaQbins, 0, alphaQmax);
    //h_alpha_charge_shift = (TH1D*)h_alpha_charge->Clone("alpha_charge_shift");
    h_alpha_charge_bg   = (TH1D*)h_alpha_charge->Clone("alpha_charge_bg");
    h_alpha_charge_sub = (TH1D*)h_alpha_charge->Clone("alpha_charge_sub");
    h_alpha_charge_sub  ->SetTitle("Candidate alphas after background subtraction");
    
    float betaQmax = 100e3;
    int   betaQbins = 50;
    h_beta_charge      = new TH1D("beta_charge","Candidate betas;Collection Plane Charge [e];Events", betaQbins, 0, betaQmax);
    //h_beta_charge_shift = (TH1D*)h_beta_charge->Clone("beta_charge_shift");
    h_beta_charge_bg = (TH1D*)h_beta_charge->Clone("beta_charge_bg");
    h_beta_charge_sub = (TH1D*)h_beta_charge->Clone("beta_charge_sub");
    h_beta_charge_sub  ->SetTitle("Candidate betas after background subtraction");
   
    int timeBins = binPeriodMax[fConfig]/binPeriodHrs[fConfig];
    float timeMax = binPeriodMax[fConfig];
    h_time_vs_rate_bipo = new TH1D("time_vs_rate_bipo",";Time [hr];Rate per 3.2 ms readout",timeBins,0,timeMax);
    h_time_vs_activity   = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_activity");
    h_time_vs_activity    ->GetYaxis()->SetTitle("Equivalent activity [mBq/kg]");
    h_time_vs_rate_bg   = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_rate_BG");
    h_time_vs_rate_bg   ->SetTitle("Background component");
   
    //for(int i=0; i<timeBins;i++){
      //h_poszy[i]        = new TH2D(Form("zy_%i",i),      "Bipo candidates;Z [cm]; Y [cm]", Zbins/10,  Zmin, Zmax, Ybins/10,  Ymin, Ymax);
      //h_poszy[i]        ->SetOption("colz");
      //h_poszy_shift[i]  = new TH2D(Form("zy_shift_%i",i),"Bipo BG;Z [cm]; Y [cm]",         Zbins/10,  Zmin, Zmax, Ybins/10,  Ymin, Ymax);
      //h_poszy_shift[i]  ->SetOption("colz");
      //h_poszy_sub[i]    = new TH2D(Form("zy_sub_%i",i),  "Post subtraction",               Zbins/10,  Zmin, Zmax, Ybins/10,  Ymin, Ymax);
      //h_poszy_sub[i]    ->SetOption("colz");
    //}

    if( _isMC ) {
      h_true_alpha_depne  = (TH1D*)h_alpha_charge->Clone("true_alpha_depne");
      h_true_alpha_depne  ->SetTitle("True ionization electrons from alpha");
      h_true_alpha_charge = (TH1D*)h_alpha_charge->Clone("true_alpha_charge");
      h_true_alpha_charge ->SetTitle("True alpha charge at anode");
      h_matched_alpha_charge = (TH1D*)h_alpha_charge->Clone("matched_alpha_charge");
      h_matched_alpha_charge ->SetTitle("Cluster charge matched to alpha");
      h_true_attenuation    = new TH1D("true_attenuation","Electron attenuation at truth-level;Calculated lifetime [#mus]",500,0,100e3);
    }

    // =====================================================
    // Diagnotic and utility histograms
    tdir_util->cd();
    h_2D_time_vs_dT      = new TH2D("2D_time_vs_dT",";Time [hr];#DeltaT [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_2D_time_vs_dT_bg= (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_bg");
    //h_2D_time_vs_dT_shift= (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_shift");
    //h_2D_time_vs_dT_flip= (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_flip");
    //h_2D_time_vs_dT_sub  = (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_sub");
    h_time_vs_N       = new TH1D("time_vs_N",";Time [hr];Number of entries into dT plot",timeBins,0,timeMax);
    h_time_vs_ratio     = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_ratio");
    h_time_vs_ratio     ->SetTitle("Signal-to-background ratio");
    h_time_vs_ratio     ->GetYaxis()->SetTitle("Signal-to-background ratio");
    h_timeFine_vs_evts     = new TH1D("timeFine_vs_evts",";Event time [hr];Number of events",timeBins,0,timeMax);
    h_timeFine_vs_rate  = new TH1D("timeFine_vs_rate",";Event time [hr];Number of candidates per readout",timeBins,0,timeMax);
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
    if( _isMC && fIgnoreTrueAlphas ) {
      printf("******** WARNING: IGNORING TRUE ALPHAS ***********************\n");
    }

    // open the file and set up the TTree
    std::string   _fileName = "files/" + inFile.fileName;
    TFile* file = new TFile(_fileName.c_str(),"READ");
    fTree = (TTree*)file->Get(inFile.treeName.c_str());
  
    // set branches
    fTree->SetBranchAddress("timestamp",&timestamp);
    fTree->SetBranchAddress("nclusts",&nclusts);                     
    fTree->SetBranchAddress("clust_nwires",&clust_nwires);
    fTree->SetBranchAddress("clust_plane",&clust_plane);              
    fTree->SetBranchAddress("clust_wire",&clust_wire);                
    fTree->SetBranchAddress("clust_startwire",&clust_startwire);                
    fTree->SetBranchAddress("clust_endwire",&clust_endwire);                
    fTree->SetBranchAddress("clust_nhits",&clust_nhits);              
    fTree->SetBranchAddress("clust_charge",&clust_charge);            
    fTree->SetBranchAddress("clust_time",&clust_time);                
    fTree->SetBranchAddress("clust_blipid",&clust_blipid);            
    fTree->SetBranchAddress("nblips",&nblips);                       
    fTree->SetBranchAddress("blip_nplanes",&blip_nplanes);
    fTree->SetBranchAddress("blip_sigmayz",&blip_sigmayz);
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
      fTree->SetBranchAddress("part_trackID",&part_trackID);
      fTree->SetBranchAddress("part_pdg",&part_pdg);
      fTree->SetBranchAddress("part_startT",&part_startT);
      //fTree->SetBranchAddress("part_depElectrons",&part_depElectrons);
      //fTree->SetBranchAddress("part_numElectrons",&part_numElectrons);
      fTree->SetBranchAddress("edep_g4id",&edep_g4id);
      fTree->SetBranchAddress("edep_pdg",&edep_pdg);
      fTree->SetBranchAddress("edep_electrons",&edep_electrons);
      fTree->SetBranchAddress("edep_charge",&edep_charge);
      fTree->SetBranchAddress("edep_tdrift",&edep_tdrift);
      fTree->SetBranchAddress("edep_x",&edep_x);
      fTree->SetBranchAddress("edep_y",&edep_y);
      fTree->SetBranchAddress("edep_z",&edep_z);
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
    if( fPickyBlipMode ) {
      _betaMinPlanes    = 3;      // Min number of matched planes (must be 2 or 3)
      _betaMaxDiff      = 2;      // Difference in wire intersection points [cm]
    } 
  
    if( fBackgroundMode == 1 ) _minTick = (int)fdT_max*2;
    _liveTimePerEvt  = (_maxTick-_minTick)*fSamplePeriod*1e-6; //sec

    // Keep list of cluster counts per wire/channel
    std::map<int,int> map_wn;
    
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
      _numEvents++;
      
      // Record the event time relative to start of dataset period
      double eventHr = ( timestamp - inFile.t0 ) / 3600.;
      if( !_isMC && (eventHr < fMinHr || eventHr > fMaxHr) ) continue;
      h_time_vs_N->Fill(eventHr);
      h_timeFine_vs_evts->Fill(eventHr);
      
      // ===============================================================
      // Clear masks to track blip/cluster availability to avoid double-counts
      // ===============================================================
      //_blipAvailable.assign(nblips, true);
      _clustAvailable.assign(nclusts, true);
      
      // ======================================
      // Check truth info
      // ======================================
      std::vector<int> alpha_edeps;
      std::vector<int> beta_edeps;
      int alphaPDG = 1000020040;
      for(int i=0; i<nedeps; i++){
        int g4index = FindG4Index(edep_g4id[i]);
        float q_dep   = edep_electrons[i];
        float q_drift = edep_charge[i];
        if( !part_isPrimary[g4index] ) continue;
        //std::cout<<"Edep "<<i<<"   PDG "<<edep_pdg[i]<<"   G4id "<<edep_g4id[i]<<"   qdep "<<q_dep<<"\n";
        if( edep_pdg[i] == alphaPDG ) {
          if( q_drift > 0 ) h_true_alpha_charge ->Fill(q_drift/0.826);
          if( q_dep > 0 )   h_true_alpha_depne  ->Fill(q_dep);
          //float driftTime = edep_x[i] / 0.109793;
          float driftTime = edep_tdrift[i];
          float readoutTime = driftTime + part_startT[i];
          int   readoutTick = readoutTime/fSamplePeriod;
//          std::cout<<"Found alpha, drift time "<<edep_tdrift[i]<<"\n";
          // does time fall within readout tick window
          if( readoutTick > _minTick && readoutTick < _maxTick ) {
            _numBiPo_true_perfectReco++;
          }
        }
      }
      
      // look for collection plane clusts matched to an alpha
      std::vector<bool> clust_isAlpha(nclusts,false);
      std::vector<bool> clust_isBeta(nclusts,false);
      std::vector<bool> clust_isGamma(nclusts,false);
      for(int i=0; i < nclusts; i++){
        if( clust_plane[i] != 2 ) continue;
        int eid = clust_edepid[i];
        if( eid < 0 ) continue;
        int g4index = FindG4Index(edep_g4id[eid]);
        if( part_isPrimary[g4index] ) {
          if(part_pdg[g4index] == alphaPDG )  clust_isAlpha[i] = true;
          if(part_pdg[g4index] == 11)         clust_isBeta[i]  = true;
        } else {
          if(part_pdg[g4index] == 11)         clust_isGamma[i] = true;
        }
        
        if( clust_isAlpha[i] ) h_matched_alpha_charge->Fill( clust_charge[i] );
        
        // Option to ignore alpha blips for background assessment
        if( fIgnoreTrueAlphas && clust_isAlpha[i] ) _clustAvailable[i] = false;
        if( fIgnoreTrueGammas && clust_isGamma[i] ) _clustAvailable[i] = false;
      }



      // ====================================================
      // Map of clust IDs per wire on collection plane
      // ====================================================
      _map_wire_clusters.clear();
      for(int i=0; i < nclusts; i++){
        if( clust_plane[i] != 2 ) continue;
        for(int j=clust_startwire[i]; j<=clust_endwire[i]; j++)
          _map_wire_clusters[clust_wire[i]].push_back(i);
        map_wn[clust_wire[i]] += 1;
        h_wt_clusts->Fill(clust_wire[i],clust_time[i]);
      }

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
      

      // ==============================================================
      // Loop over 3D blips...
      // ==============================================================
      for(auto& iBlip : sortedBlips ) {
     
        // mark this blip as used
        //_blipAvailable[iBlip] = false;
     
        // find associated cluster on collection
        int   ic    = blip_clustid[2][iBlip];
  
        // plot ZY position
        h_zy_blips->Fill( blip_z[iBlip], blip_y[iBlip] );
        h_wt_blips->Fill( clust_wire[ic], clust_time[ic] );
  
        // skip if this cluster was already included in a BiPo candidate 
        if( !_clustAvailable[ic] ) continue;
        
  			// skip if any wire in cluster is from a noisy wire
        // skip if this cluster is at the very edge of the wireplane
        int w_start = std::max(0,clust_startwire[ic]-1);
        int w_end   = std::min(clust_endwire[ic]+1,(int)_nwires-1);
        if( w_start < 10 || w_end > int(_nwires-10) ) continue;
        bool flag = false;
        for(int iwire = w_start; iwire <= w_end; iwire++){
          if( wireIsNoisy[iwire] || wireIsBad[iwire] ) {
            flag = true;
            break;
          }
        }

        if( flag ) continue;
        

        // 3D match cuts for beta
        bool isWorthy = (blip_nplanes[iBlip] >= _betaMinPlanes && blip_sigmayz[iBlip] < _betaMaxDiff);
        
        if((blip_nplanes[iBlip] == 3 && blip_sigmayz[iBlip] < 3 ) ) h_zy_blips_picky->Fill( blip_z[iBlip], blip_y[iBlip] );
        
        if( !isWorthy ) continue;
        
        // fill some blip/cluster location histograms
        h_zy_blips_filt->Fill( blip_z[iBlip], blip_y[iBlip] );
        h_wt_blips_filt->Fill( clust_wire[ic], clust_time[ic] );
        h_blip_charge->Fill(clust_charge[ic]);
        h_blip_nhits->Fill(clust_nhits[ic]);

        // apply charge/size cuts on beta
        if(   clust_charge[ic] < fBetaCharge_min 
          ||  clust_charge[ic] > fBetaCharge_max ) continue;
        
        // skip if adjacent wires are missing
        if( wireIsBad[clust_wire[ic]+1] ) continue;
        if( wireIsBad[clust_wire[ic]-1] ) continue;
       

        // evaluate if in fiducial volume
        if( fFidVolCut ) {
          if( blip_z[iBlip] < fZmin || blip_z[iBlip] > fZmax ) continue;
          if( blip_y[iBlip] < fYmin || blip_y[iBlip] > fYmax ) continue;
        }
      
        // check for nearby 3D blips
        //int nMultSphere3D = 0;
        //for(int ii=0; ii<nblips; ii++){
        //  if( ii == iBlip ) continue;
        //  if( blip_nplanes[ii] < _betaMinPlanes || blip_maxdiff[ii] > _betaMaxDiff ) continue;
        //  TVector3 p1(blip_x[iBlip],blip_y[iBlip],blip_z[iBlip]);
        //  TVector3 p2(blip_x[ii],blip_y[ii],blip_z[ii]);
        //  float ds = (p2-p1).Mag();
        //  if(ds<fSphereRadius) nMultSphere3D++; 
        //}
        //h_blip_sphere_mult->Fill(nMultSphere3D);
        //if( fMinSphereMult3D > 0 && nMultSphere3D < fMinSphereMult3D ) continue;
        //if( fMaxSphereMult3D > 0 && nMultSphere3D > fMaxSphereMult3D ) continue;

        // skip if we are near end of wire (account for 400us/800tick trigger offset)
        float peakT = clust_time[ic]+800;
        if( peakT < 0 ) continue;
        if( peakT > _maxTick ) continue;
        if( peakT < _minTick ) continue; 
       

        // ------------------------------------------------------------
        // Begin search for alpha candidates
        // ------------------------------------------------------------
        int nclusts_inwindow        = 0;
        int nclusts_inwindow_shift  = 0;
        int nclusts_inwindow_shift2  = 0;

        // -------------------------------------------
        // Search for standard candidates
        int nwires = 0;
        std::vector<BiPoCandidate> v_cands = FindCandidates(ic, 0, false, nclusts_inwindow,nwires);
        h_nclusts_inwindow->Fill(nclusts_inwindow);
        h_ncands_inwindow->Fill((int)v_cands.size());
        
        // -------------------------------------------
        // Search for background candidates (wire shift)
        int nwires_shift = 0;
        int shift = 2*fWireRange+1;
        std::vector<BiPoCandidate> v_cands_shift = FindCandidates(ic,shift,false,nclusts_inwindow_shift,nwires_shift);

        // ------------------------------------------------------------
        // Search for background candidates (same wire, but dT-flip)
        int nclusts_inwindow_flip;
        int nwires_flip = 0;
        std::vector<BiPoCandidate> v_cands_flip = FindCandidates(ic, 0, true, nclusts_inwindow_flip,nwires_flip);
      
        // --------------------------------------------
        // Evaluate standard candidates
        // ---------------------------------------------
        if( v_cands.size() && v_cands.size()<=fMaxCands && nclusts_inwindow <=  fMaxClustMult ) {
          
          // flag this beta candidate as unavailable
          _clustAvailable[ic] = false;

          std::vector<BiPoCandidate> v_cands_dtcut;
          bool flag = false;
          for(auto& thisCand : v_cands) {
            
            _clustAvailable[thisCand.id2] = false;
            //int blipid = clust_blipid[thisCand.id2];
            //if( blipid >= 0 ) _blipAvailable[blipid] = false;
            // count how many are withiin dT window
            if( thisCand.dT > fdT_min ) v_cands_dtcut.push_back( thisCand );
            // radon mitigation paper-like count
            if( !flag && thisCand.dT > 6.25 && thisCand.dT < 312.5 ) {
              flag = true;
              _numBiPo_6_312++;
              if( clust_isBeta[ic] ) _numBiPo_6_312_true++;
              h_timeFine_vs_rate  ->Fill(eventHr);
            }
          }
          
          
          if( v_cands_dtcut.size() ) {
            _numBiPo++;
            if( clust_isBeta[ic] ) _numBiPo_true++;
            // plot locations
            h_zy_bipos->Fill( blip_z[iBlip], blip_y[iBlip]);
            h_wt_bipos->Fill( clust_wire[ic], clust_time[ic]);
            float weight = 1./v_cands_dtcut.size();
            for(auto& thisCand : v_cands_dtcut) {
              h_beta_charge         ->Fill(thisCand.q1,weight);
              h_alpha_charge        ->Fill(thisCand.q2,weight);
              h_cand_dT             ->Fill(thisCand.dT,weight);
              h_2D_time_vs_dT       ->Fill(eventHr,thisCand.dT,weight);
              h_alpha_nwires        ->Fill(clust_nwires[thisCand.id2],weight);
              //h_poszy[eventTimeBin] ->Fill(blip_z[thisCand.blipID],blip_y[thisCand.blipID]);
            }
            
          }
         
          
          /*
          auto& thisCand = v_cands.at(0);
          
          // plot locations
          h_zy_bipos->Fill( blip_z[iBlip], blip_y[iBlip] );
          h_wt_bipos->Fill( clust_wire[ic], clust_time[ic] );
         
          _clustAvailable[thisCand.id1] = false;
          _clustAvailable[thisCand.id2] = false;
          int blipid = clust_blipid[thisCand.id2];
          if( blipid >= 0 ) _blipAvailable[blipid] = false;
        
          if( thisCand.dT > 6.25 && thisCand.dT < 312.5 ) {
            _numBiPo_6_312++;
            if( clust_isBeta[ic] ) _numBiPo_6_312_true++;
            h_timeFine_vs_rate  ->Fill(eventHr);
          }
          if( thisCand.dT >= fdT_min ) {
            _numBiPo++;
            if( clust_isBeta[ic] ) _numBiPo_true++;
            h_beta_charge         ->Fill(thisCand.q1);
            h_alpha_charge        ->Fill(thisCand.q2);
            h_cand_dT             ->Fill(thisCand.dT);
            h_2D_time_vs_dT       ->Fill(eventHr,thisCand.dT);
            h_alpha_nwires        ->Fill(clust_nwires[thisCand.id2]);
            //h_poszy[eventTimeBin] ->Fill(blip_z[thisCand.blipID],blip_y[thisCand.blipID]);
          }
          */

        }//end evaluation of standard cands
        
        // --------------------------------------------
        // Evaluate wire shift candidates
        // ---------------------------------------------
        if( v_cands_shift.size() && v_cands_shift.size()<=fMaxCands && nclusts_inwindow_shift <=  fMaxClustMult ) {

          std::vector<BiPoCandidate> v_cands_dtcut;
          for(auto& thisCand : v_cands_shift) {
            // count how many are withiin dT window
            if( thisCand.dT >= fdT_min ) v_cands_dtcut.push_back( thisCand );
          }
          
          if( v_cands_dtcut.size() ) {
            float weight_n = nwires/float(nwires_shift);
            float weight = weight_n/v_cands_dtcut.size();
            for(auto& thisCand : v_cands_dtcut) {
              h_cand_dT_shift->Fill(thisCand.dT,weight);
              if( fBackgroundMode == 0 ) {
                h_2D_time_vs_dT_bg       ->Fill(eventHr,thisCand.dT,weight);
                h_beta_charge_bg         ->Fill(thisCand.q1,weight);
                h_alpha_charge_bg        ->Fill(thisCand.q2,weight);
              }
            }
          }
        }
        
        /*
        // --------------------------------------------
        // Evaluate wire shift candidates
        // ---------------------------------------------
        //float weight = 1./float(n_regions_checked);
        //std::cout<<"applying weight "<<weight<<" to background\n";
        float weight = nwires/float(nwires_shift);
        if( v_cands_shift.size() && v_cands_shift.size() ==1 && nclusts_inwindow_shift <=  fMaxClustMult ) {
          for(auto& thisCand : v_cands_shift ) {
            if( thisCand.dT >= fdT_min ) {
              h_cand_dT_shift             ->Fill(thisCand.dT,weight);
              if( fBackgroundMode == 0 ) {
                h_2D_time_vs_dT_bg       ->Fill(eventHr,thisCand.dT,weight);
                h_beta_charge_bg         ->Fill(thisCand.q1,weight);
                h_alpha_charge_bg        ->Fill(thisCand.q2,weight);
              }
            }
          }//endloop over candidates
        }//end evaluation of wire-shift candidates
        */
        


        // --------------------------------------------
        // Evaluate dT-flip candidates
        // ---------------------------------------------
        if( v_cands_flip.size() && v_cands_flip.size() <= fMaxCands && nclusts_inwindow_flip <=  fMaxClustMult ) {
          std::vector<BiPoCandidate> v_cands_dtcut;
          for(auto& thisCand : v_cands_flip) {
            //_clustAvailable[thisCand.id2] = false;
            // count how many are withiin dT window
            if( thisCand.dT >= fdT_min ) v_cands_dtcut.push_back( thisCand );
          }
          
          if( v_cands_dtcut.size() ) {
            float weight = 1./v_cands_dtcut.size();
            for(auto& thisCand : v_cands_dtcut) {
              h_cand_dT_flip->Fill(thisCand.dT,weight);
              if( fBackgroundMode == 1 ) {
                h_2D_time_vs_dT_bg       ->Fill(eventHr,thisCand.dT,weight);
                h_beta_charge_bg         ->Fill(thisCand.q1,weight);
                h_alpha_charge_bg        ->Fill(thisCand.q2,weight);
              }
            }
          }
        }
        
        /*
        // --------------------------------------------
        // Evaluate dT-flip candidates
        // ---------------------------------------------
        if( v_cands_flip.size() && v_cands_flip.size() ==1 && nclusts_inwindow_flip <=  fMaxClustMult ) {
          for(auto& thisCand : v_cands_flip ) {
            if( thisCand.dT >= fdT_min ) {
              h_cand_dT_flip             ->Fill(thisCand.dT);
              if( fBackgroundMode == 1 ) {
                h_2D_time_vs_dT_bg       ->Fill(eventHr,thisCand.dT);
                h_beta_charge_bg       ->Fill(thisCand.q1);
                h_alpha_charge_bg        ->Fill(thisCand.q2);
              }
            }
          }//endloop over candidates
        }//end evaluation of dt-flip candidates
      
        */
      

      }//end loop over 3D blips

    }//endloop over events
    double loopDuration = ( time(NULL) - loopStart );

   
    // ***************************************************
    // Scale dT plots so they're per readout
    // ***************************************************
    _totalLiveTime = float(_numEvents) * _liveTimePerEvt;
    //float scaleFact = 1./(float(_numEvents));
    // Scale everything so it's 'per second in AV'
    float scaleFact = 1./float(_totalLiveTime * _fiducialFrac);
    h_cand_dT         ->Scale( scaleFact );
    h_cand_dT_shift   ->Scale( scaleFact );
    h_cand_dT_flip    ->Scale( scaleFact );
    //h_2D_time_vs_dT   ->Scale( scaleFact );
    //h_2D_time_vs_dT_shift ->Scale( scaleFact );
    //h_2D_time_vs_dT_flip  ->Scale( scaleFact );
    h_timeFine_vs_rate  ->Divide( h_timeFine_vs_evts );

    // ***************************************************
    // Histogram subtraction time!
    // ***************************************************
    h_cand_dT_sub         ->Add(h_cand_dT,        1);
    //h_2D_time_vs_dT_sub   ->Add(h_2D_time_vs_dT,  1);
    h_alpha_charge_sub    ->Add(h_alpha_charge,   1);
    h_beta_charge_sub     ->Add(h_beta_charge,    1);
    
    if( fBackgroundMode == 0 ) {
      h_cand_dT_sub       ->Add(h_cand_dT_shift,      -1.);
      //h_2D_time_vs_dT_sub ->Add(h_2D_time_vs_dT_shift,-1.);
    } else
    if( fBackgroundMode == 1 ) {
      h_cand_dT_sub       ->Add(h_cand_dT_flip,       -1.);
      //h_2D_time_vs_dT_sub ->Add(h_2D_time_vs_dT_flip, -1.);
    }

    h_alpha_charge_sub->Add(h_alpha_charge_bg,-1.);
    h_beta_charge_sub->Add(h_beta_charge_bg,-1.);

    /*
    h_cand_dT         ->Scale( scaleFact );
    h_cand_dT_shift   ->Scale( scaleFact );
    h_cand_dT_flip    ->Scale( scaleFact );
    //h_2D_time_vs_dT   ->Scale( scaleFact );
    //h_2D_time_vs_dT_shift ->Scale( scaleFact );
    //h_2D_time_vs_dT_flip  ->Scale( scaleFact );
    h_timeFine_vs_rate  ->Divide( h_timeFine_vs_evts );
    h_cand_dT_sub->Scale(scaleFact);
    h_2D_time_vs_dT_sub ->Scale(scaleFact);
    */

    //for(int i=0; i<timeBins; i++){
      //h_poszy_sub[i]->Add(h_poszy[i],1);
      //h_poszy_sub[i]->Add(h_poszy_shift[i],-1);
    //}
   
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
    printf("dT min/max          : %.2f-%.2f us\n",  fdT_min,fdT_max);
    printf("Ave cands / evt     : %f (%f in 6.25-312.us)\n",h_cand_dT->GetEntries()/(float)_numEvents, _numBiPo_6_312/(float)_numEvents );
    printf("  - truth-matched   : %f (%f in 6.25-312.us)\n",_numBiPo_true/(float)_numEvents, _numBiPo_6_312_true/(float)_numEvents );
    printf("BG subtrctn method  : %i\n",fBackgroundMode);
    printf("Processing time     : %f sec (%f sec/evt)\n", loopDuration, loopDuration/float(_numEvents));
    printf("Excluding %i noisy wires \n",     (int)fNoisyWires.size());	
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
  
    std::string name;
    float       range;
    float       min;
    float       max;
   
    if( !_isMC ) {

      // ============================================
      // Do slice-by-slice dT fit
      // ============================================
      TH1D* h_time_vs_p0 = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_p0");
      TH1D* h_time_vs_p1 = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_p1");

      TH1D* h_slice;
      TH1D* h_bg;
      int nbins = h_time_vs_N->GetXaxis()->GetNbins();
      for(int i=1; i<=nbins; i++){
        //std::cout<<"SLICE "<<i<<"   \n";
        h_slice     = Make1DSlice( h_2D_time_vs_dT, i, i, Form("time_vs_dT_%i",i) );
        h_bg    = Make1DSlice( h_2D_time_vs_dT_bg, i, i, Form("time_vs_dT_flip_%i",i) );
        //float scaleFact = 1./(float)h_time_vs_N->GetBinContent(i);
        float scaleFact = 1./float(h_time_vs_N->GetBinContent(i)*_liveTimePerEvt);
        h_slice->Scale(scaleFact);
        h_bg->Scale(scaleFact);
        h_slice->Add( h_bg, -1. );
        float t = h_time_vs_N->GetXaxis()->GetBinCenter(i);
        float dt = h_time_vs_N->GetXaxis()->GetBinWidth(i)/sqrt(12);
        FitResult fr = fitdT(h_slice,true,false);
        if( fr.rate_signal != -9 ) {
          h_time_vs_activity  ->SetBinContent(  i,  fr.activity);
          h_time_vs_activity  ->SetBinError(    i,  fr.activity_err);
          h_time_vs_rate_bipo ->SetBinContent(  i,  fr.rate_signal);
          h_time_vs_rate_bipo ->SetBinError(    i,  fr.rate_signal_err);
          h_time_vs_rate_bg   ->SetBinContent(  i,  fr.rate_bg);
          h_time_vs_rate_bg   ->SetBinError(    i,  fr.rate_bg_err);
          h_time_vs_p0        ->SetBinContent(  i,  fr.p0);
          h_time_vs_p0        ->SetBinError(    i,  fr.p0_err);
          h_time_vs_p1        ->SetBinContent(  i,  fr.p1);
          h_time_vs_p1        ->SetBinError(    i,  fr.p1_err);
        }
      }
      
      fOutFile->cd();
      h_time_vs_activity ->Write(0, TObject::kOverwrite );
      h_time_vs_rate_bipo ->Write(0, TObject::kOverwrite );
      h_time_vs_rate_bg   ->Write(0, TObject::kOverwrite );
      h_time_vs_ratio     ->Write(0, TObject::kOverwrite );
      
      std::cout<<"Rate at 2 hours: "<<h_time_vs_rate_bipo->Interpolate(2)<<"\n";
      std::cout<<"Rate at 4 hours: "<<h_time_vs_rate_bipo->Interpolate(4)<<"\n";
      std::cout<<"Rate at 6 hours: "<<h_time_vs_rate_bipo->Interpolate(6)<<"\n";


      // ============================================
      // Plot fit parameters vs time
      // ============================================
      TH1D* h_p0  = h_time_vs_rate_bg; //h_time_vs_p0;
      TH1D* h_p1  = h_time_vs_rate_bipo; //h_time_vs_p1;
      FormatTH1D(h_p1, kBlue, 1, 2, 20, 1);
      FormatTH1D(h_p0, kRed, 1, 2, 20, 1);
      name = "c_time_vs_rate";
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
      
      tdir_plots->cd();
      c->Write();


      // ============================================
      // Plot BiPo rate vs time
      // ============================================
      TH1D* h1  = h_time_vs_rate_bipo;
      FormatTH1D(h1, kBlue+2, 1, 2, 20, 1);
   
      name = "c_time_vs_rate";
      range = (GetHistMax(h1)-GetHistMin(h1));
      TCanvas* c2 = new TCanvas(name.c_str(),name.c_str(),500,380);
      gStyle->SetOptStat(0);
      gPad->SetGridy(1);
      //h1->GetYaxis()->SetRangeUser(GetHistMin(h1)-0.3*range, GetHistMax(h1)+0.5*range);
      h1->GetYaxis()->SetTitleOffset(1.1); 
      h1->GetYaxis()->SetTitle("Decay candidates per 3.2ms readout");
      h1->DrawCopy();
      c2->Write();
      
      // ============================================
      // Plot activity vs time
      // ============================================
      TH1D* h2  = h_time_vs_activity;
      FormatTH1D(h2, kBlack, 1, 2, 20, 1);
   
      name = "c_time_vs_activity";
      range = (GetHistMax(h2)-GetHistMin(h2));
      TCanvas* c3 = new TCanvas(name.c_str(),name.c_str(),500,380);
      gStyle->SetOptStat(0);
      gPad->SetGridy(1);
      //h1->GetYaxis()->SetRangeUser(GetHistMin(h1)-0.3*range, GetHistMax(h1)+0.5*range);
      h2->GetYaxis()->SetTitleOffset(1.3); 
      h2->GetYaxis()->SetTitle("Equivalent activity [mBq/kg]");
      h2->DrawCopy();
      c3->Write();
    
    }


    // ============================================
    // Do final fit on dT spectrum
    // ============================================
    h_cand_dT_sub->GetXaxis()->SetTitle("#DeltaT [#mus]");
    h_cand_dT_sub->GetYaxis()->SetTitle("Decay candidates per sec");
    fitdT( h_cand_dT_sub, true, true );

  } 
  
  
  
  //################################################################################
  // Function that performs the dT fit
  //#################################################################################
  FitResult fitdT(TH1D* h, bool writeCanvas = false, bool constrainNorm = false ){
  
    FitResult out;
  
    std::cout<<"\n\nFitting dT spectrum "<<h->GetTitle()<<", "<<h->GetEntries()<<"\n";
    std::string label = h->GetName();
    TCanvas* c = new TCanvas(Form("c_fit_%s",label.c_str()),Form("c_fit_%s",label.c_str()),500,420);
    //gPad->SetMargin(mar_l, mar_r, mar_b, mar_t ); 
    //TH1D* hc = (TH1D*)h->Clone(Form("c_fit_%s",label.c_str()));
    TH1D* hc = (TH1D*)h->Clone();
    float histMax = GetHistMax(hc);
    float histMin = GetHistMin(hc);
    float range   = (histMax-histMin);
    std::cout<<"Hist min/max: "<<histMin<<", "<<histMax<<"\n";
   
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

    // *** single exp fit + flat BG ***
    TF1* fit = new TF1("FullFit","[0] + [1]*exp(-x/[2])",fdT_min,fdT_max);
    fit->SetParameter(0, 0. );
    //fit->SetParLimits(0, -histMax, histMax );
    fit->SetParameter(1, (histMax-histMin)/2 );
    //fit->SetParLimits(1, histMin-range, histMax+range);
    if( constrainNorm ) fit->SetParLimits(1, 0., histMax+range);
    fit->FixParameter(2, 164.3 );
  
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.0001);

    // Draw plot and fit
    c->cd();
    hc->Fit(fit,"RE"); // "WL" = log likelihood w/ weighted bins
   
    bool isAtLimit = gMinuit->fLimset;
    std::cout<<"Is MINUIT reaching limit? "<<gMinuit->fLimset<<"\n";

    // For more accurate error, let normalization dip below
    // zero and re-fit. Add the difference between this and
    // the constrained fit in quadrature with the previous error.
    // This can be considered a systematic.
    if( isAtLimit && constrainNorm ) {
      std::cout<<"Par1: "<<fit->GetParameter(1)<<" +/- "<<fit->GetParError(1)<<"\n";
      TF1* fit2 = (TF1*)fit->Clone("FullFit2");
      fit2->FixParameter(0,fit->GetParameter(0)-fit->GetParError(0));
      hc->Fit(fit2,"REN");
      double norm1 = fit->GetParameter(1);
      double norm2 = fit2->GetParameter(1);
      double diff = fabs(norm2-norm1);
      double err = sqrt( pow( diff, 2) + pow( fit2->GetParError(1),2 ) );
      std::cout<<"new norm: "<<norm2<<" +/- "<<fit2->GetParError(1) <<"\n"; 
      fit->SetParError(1,err);
      std::cout<<"new par1: "<<fit->GetParameter(1)<<" +/- "<<fit->GetParError(1)<<"\n";
    }


    hc->DrawCopy();
    //fit->DrawCopy("SAME");
    //
    gStyle->SetOptFit(1112);
    gPad->Modified(); gPad->Update();
    //gPad->SetLogy();
    if( writeCanvas ) {
      tdir_util->cd();
      hc->Write(0, TObject::kOverwrite );
      tdir_plots->cd();
      c->Write();
      fOutFile->cd();
    }
  
    // Factor out the signal and BG components into
    // separate TF1 functions
    TF1* f_exp = new TF1("expGeneric","[0]*exp(-x/[1])");
    //f_exp->SetParameter(0,1); 
    //f_exp->SetParameter(1,999);
    
    TF1* f_flat = new TF1("flat","[0]");
    f_flat->FixParameter(0, fit->GetParameter(0));
  
    //TF1* f_expBG = (TF1*)f_exp->Clone("expBG");
    //f_expBG->FixParameter(0, fit->GetParameter(3) );
    //f_expBG->FixParameter(1, fit->GetParameter(4) );
    
    TF1* f_bipo = (TF1*)f_exp->Clone("bipo");
    f_bipo->FixParameter(0, fit->GetParameter(1) );
    f_bipo->FixParameter(1, fit->GetParameter(2) );
    
    // Full extrapolated BiPo component; integral of A*exp(-x/B) from 0-infinity is A*B
    double A       = fit->GetParameter(1);
    double dA      = fit->GetParError(1);
    double B       = fit->GetParameter(2);
    double dB      = fit->GetParError(2);
    double n_bipo  = (A*B)/fdT_binSize;
    double n_bipo_err = sqrt( pow(dA/A,2)+pow(dB/B,2) )*fabs(n_bipo);
   

    // Components within the dT selection window
    double N_total = hc      ->Integral(0,hc->GetXaxis()->GetNbins());
    double N_fit   = fit     ->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_flat  = f_flat  ->Integral(fdT_min,fdT_max)/fdT_binSize;
    //double N_expbg = f_expBG ->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_bipo  = f_bipo  ->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_bg    = N_flat; //N_expbg;
    double N_bg_err = fit->IntegralError(fdT_min,fdT_max)/fdT_binSize;
    double sbratio = N_bipo / N_bg;

    //double activity_mBq = 1e3 * _fidCorFactor * (n_bipo/_liveTimePerEvt ) / 85000.;
    //double activity_mBq_err = (n_bipo_err/n_bipo)*activity_mBq;
   
    // Error terms
    double err1 = (n_bipo!=0) ? n_bipo_err/n_bipo : 0.;
    double err2 = (_effMC>0) ? _effMC_err/_effMC : 0.;

    double activity_mBq = 1e3 * (1./_effMC) * n_bipo / 85000.; 
    double activity_mBq_err = activity_mBq * sqrt( pow(err1,2) + pow(err2,2) ); 
    
    // normalize "rate" to be per 3.2ms readout
    //double rate_norm_factor = (3.2e-3)/_liveTimePerEvt;
    double rate_norm_factor = 0.0032;
    n_bipo *= rate_norm_factor;
    n_bipo_err *= rate_norm_factor;
    
    printf("================ dT fit =================\n");
    printf("p0                  : %f +/- %f\n", fit->GetParameter(0),fit->GetParError(0));
    printf("p1                  : %f +/- %f\n", fit->GetParameter(1),fit->GetParError(1));
    printf("p2                  : %f +/- %f\n", fit->GetParameter(2),fit->GetParError(2));
    printf("Chi2/ndf            : %f\n",        fit->GetChisquare()/fit->GetNDF());
    printf("Total entries       : %f\n",        hc->GetEntries());
    printf("Total rate          : %f per sec\n",N_total);
    printf(" - BiPo             : %f per sec\n",N_bipo);
    printf(" - Flat             : %f per sec\n",N_flat);
    printf("Extrapolated rate   : %f +/- %f BiPos per 3.2ms readout\n",n_bipo,n_bipo_err);
    printf("Equivalent activity : %f +/- %f mBq/kg\n",activity_mBq,activity_mBq_err);
    printf("(assumed eff = %f)\n",_effMC);


    out.rate_signal       = n_bipo;
    out.rate_signal_err   = n_bipo_err;
    out.rate_bg           = N_bg;
    out.rate_bg_err       = N_bg_err;
    out.ratio             = sbratio;
    out.p0                = fit->GetParameter(0);
    out.p0_err            = fit->GetParError(0);
    out.p1                = fit->GetParameter(1);
    out.p1_err            = fit->GetParError(1);
    out.activity          = activity_mBq;
    out.activity_err      = activity_mBq_err;
    return out;
  
  }
  

  //################################################################################
  // Function to search for alpha candidates relative to a beta candidate
  //#################################################################################
  std::vector<BiPoCandidate> FindCandidates(int ic, int wire_shift, bool flipdT, int& npileup, int& nwires ) {
    
    //std::cout<<"Searching for cands: "<<(int)flipdT<<"   "<<start<<"-"<<end<<"\n";
    std::vector<BiPoCandidate> v;
    std::vector<int> cand_IDs;
    npileup = 0;
    
    int w_start = clust_startwire[ic] - wire_shift;
    int w_end   = clust_endwire[ic] + wire_shift;

    std::set<int> wires;

    for(int i=w_start-fWireRange; i<=w_start+fWireRange; i++)
      wires.insert(i);

    if( clust_startwire[ic] != clust_endwire[ic] ) {
    
      for(int i=w_end-fWireRange; i<w_end+fWireRange; i++)
        wires.insert(i);
    }

    nwires = wires.size();

    for(auto iWire : wires ) {
    
      for(auto& jc : _map_wire_clusters[iWire] ) {
        if( ic == jc ) continue;
        if( !_clustAvailable[jc] ) continue;
        if( std::count(cand_IDs.begin(),cand_IDs.end(),jc) ) continue;

        //... temporary ... 
        //exclude electrons
        //int eid = clust_edepid[jc];
        //if( eid >= 0 ) {
        //  //std::cout<<"Found an alpha cand that is matched "<<eid<<"   pdg "<<edep_pdg[eid]<<"\n";
        //  if( fabs(edep_pdg[eid]) == 11 ) continue;
        //}
        //
        

        // 
        float dT  = (clust_time[jc]-clust_time[ic])*fSamplePeriod;
        
        // use drift velocity to convert wire offset to 'time'
        // and create a 2D point


        //std::cout<<"   clust "<<jc<<" has dT "<<dT<<"\n";
        
        if (flipdT) {
          dT *= -1.;
          //std::cout<<"      flipped: "<<dT<<"\n";
        }

        // 3D plane-match requirement for alpha
        if( fAlphaReq3D && clust_blipid[jc] < 0 ) continue;
        
        // Count clusters in forward window, and fill some pre-cut histos

        if( dT > 0 && dT < fdT_max ) {
          npileup++;

            // --- alpha charge/nhits cut ---
            if(   clust_charge[jc] > fAlphaCharge_min 
              &&  clust_charge[jc] < fAlphaCharge_max 
              &&  clust_nwires[jc] <= fAlphaWires_max ) {
              
              int blipID = clust_blipid[ic];
              BiPoCandidate c = { blipID, ic, jc, dT, clust_charge[ic], clust_charge[jc]};
              v.push_back(c);
              cand_IDs.push_back(jc);
            }
          
         
        } 
  
      }//<-- end loop over clusters on this wire
    }//endloop over wires
  
    return v;
  
  }

  //###################################################################
  int FindG4Index(int g4id) {
    for(size_t i=0; i<nparticles; i++) {
      if( part_trackID[i] == g4id ) return i;
    }
    return -9;
  }

