  //////////////////////////////////////////////////////////////////
  // 
  //  Analysis ROOT Macro
  //
  ////////////////////////////////////////////////////////////////// 
 
  #include "core/vars.h"
  #include "core/tools.h"
  #include <time.h>
  
  // --- Choose run configuration ---
  int   fConfig  = 2;
  bool  fIsMC    = false;

  // Filenames 
  std::string fFileName[3]= { 
      "BlipAna_BiPo_MC_20220503.root",
      "BlipAna_RadonData_Phase1_20220503.root",
      "BlipAna_RadonData_Phase2_20220503.root"
  };
 

  //TTree* tree = new TTree("bipo_tree","bipo_tree");

  // Time periods in UNIX time
  unsigned int fT0[3] = {  0, 1627415210, 1627594369 };
  unsigned int fT1[3] = {  0, 1627592728, 1627761265 };
  //                            Phase 1     Phase 2
  //                            ~49.3 hrs   ~46.4 hrs

  // --- General macro parameters ---
  bool  fPickyBlipMode    = false;  // Require blips match on all 3 planes
  bool  fRequireHitFit    = false;  // Exclude pulse-train hits
  int   fWireRange        = 1;      // +/- range to look for alpha candidate
  float fBetaCharge_min   = 3.5e3;  // Min charge of beta candidate blip [e-]
  float fBetaCharge_max   = 80e3;   // Max charge of beta candidate blip [e-]
  float fAlphaCharge_min  = 0e3;    // Min charge of alpha candidate cluster [e-]
  float fAlphaCharge_max  = 8e3;    // Max charge of alpha candidate cluster [e-]
  int   fAlphaHits_max    = 2;      // Max hits in alpha candidate cluster
  bool  fAlphaReq3D       = false;  // Require alpha pulse be 3D-matched
  float fdT_binSize       = 20.;    // Bin width for all dT spectra plots [us]
  float fdT_min           = 20.;    // Min dT for looking for candidate [us]
  float fdT_max           = 800.;   // Max dT for looking for candidate [us]
  int   fMaxClustMult     = 999;      // Max number of clusters in time window
  float fSphereRadius     = 25;    
  bool  fFidVolCut        = false;
  float fZmin             = 50; //50;     // Z range (0 to 1037 cm)
  float fZmax             = 985; //1035; //985;    //
  float fYmin             = -80; //-120; //-70;    // Y range (-120 to 120 cm)
  float fYmax             = 80; //120;     //
  
  // --- Detector properties ---
  int   nWiresColl        = 3455;
  float fSamplePeriod     = 0.5; // microseconds
  
  // --- Special switches ---
  bool  fDoWireDiagnostics= false; //true;
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
  
  // File names and T0
  std::string   _fileName = "files/" + fFileName[fConfig];
  unsigned int  _t0       = fT0[fConfig];
  bool          _isMC     = fIsMC;
 
  // Fiducial vol correction factor
  float dz = fZmax-fZmin;
  float dy = fYmax-fYmin;
  float _fidCorFactor = std::max(1., (1037.*230.) / (dz*dy) );
  
  // Counters
  int   _numEvents      = 0;
  int   _numBiPo        = 0;
  int   _numBiPo_6_312  = 0;
  int   _numBiPo_true   = 0;
  int   _numBiPo_6_312_true = 0;

  float _totalLiveTime  = 0;
  std::vector<bool> _blipAvailable;
  std::vector<bool> _clustAvailable;

  // Map of cluster IDs per wire on collection plane
  std::map<int,std::vector<int>> _map_wire_clusters;

  //std::vector<bool> wireIsNoisy(nWiresColl, false);
  //std::vector<bool> wireIsBad(nWiresColl,   false);
  std::unordered_map<int,bool> wireIsNoisy;
  std::unordered_map<int,bool> wireIsBad;

  // Live time
  int   _minTick         = 0;
  int   _maxTick         = 6400 - (int)fdT_max*2;
  float _liveTimePerEvt  = (_maxTick-_minTick)*fSamplePeriod*1e-6; //sec
  
  // Blip fidelity requirements
  int _betaMinPlanes = 2;
  int _betaMaxDiff = 9999;
  
  
  //##########################################################################
  // Structs, functions, histograms
  //##########################################################################
  
  // Struct to hold fit results
  struct FitResult { float rate_signal = -9, rate_signal_err = 0, rate_bg = -9, rate_bg_err = 0, ratio = -9; };
  
  // useful structure to save candidate info in
  struct BiPoCandidate { int blipID, id1, id2; float dT, q1, q2; };
  
  // Functions 
  void                        makePlots();
  void                        makeHistograms();
  std::vector<BiPoCandidate>  FindCandidates(int, int, int, bool, int&);
  std::vector<BiPoCandidate>  FindCandidatesCluster(int, bool, int&, int&);
  FitResult                   fitdT(TH1D*,bool);
  int                         FindG4Index(int);
  
  // ROOT objects
  TTree*      fTree;
  TFile*      fOutFile;

  // Histograms
  TDirectory* tdir_util;
  TDirectory* tdir_plots;
  TH1D* h_blip_charge;
  TH1D* h_blip_nhits;
  TH1D* h_alpha_nhits;
  TH1D* h_nclusts_wire_ave;
  TH1D* h_nclusts_inwindow;
  TH1D* h_ncands_inwindow;
  TH1D* h_cand_dT;
  TH1D* h_cand_dT_shift;
  TH1D* h_cand_dT_sub;
  TH1D* h_clust_dT;
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
  TH2D* h_2D_time_vs_dT_shift;
  TH2D* h_2D_time_vs_dT_sub;

  TH1D* h_time_vs_rate_bipo; 
  TH1D* h_time_vs_rate_bg;
  TH1D* h_time_vs_ratio;
 
  TH1D* h_timeFine_vs_evts;
  TH1D* h_timeFine_vs_rate;
  
  //TH2D* h_time_vs_ph;
  TH2D* h_time_vs_ntrks;
  //TH2D* h_time_vs_trklen;
  
  TH1D* h_alpha_charge;
  TH1D* h_alpha_charge_shift;
  TH1D* h_alpha_charge_sub;
  TH1D* h_beta_charge;
  TH1D* h_beta_charge_shift;
  TH1D* h_beta_charge_sub;

  const int timeBins = 9;
  float timeMax = 45;
  
  TH2D* h_poszy[timeBins];
  TH2D* h_poszy_shift[timeBins];
  TH2D* h_poszy_sub[timeBins];

  float mar_l = 0.12;
  float mar_r = 0.08;
  float mar_b = 0.12;
  float mar_t = 0.08;

  // MC-truth histograms
  TH1D* h_true_alpha_depne;
  TH1D* h_true_alpha_charge;
  TH1D* h_matched_alpha_charge;

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
    h_alpha_nhits       = new TH1D("alpha_nhits","Candidate alphas;Collection Plane Hits",20,0,20);
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
   
    int dTbins = fdT_max / fdT_binSize;
    h_cand_dT       = new TH1D("cand_dT","Selected BiPo Candidates;Time difference [#mus];Candidates per readout", dTbins,0,fdT_max);
    h_cand_dT_shift = (TH1D*)h_cand_dT->Clone("cand_dT_shift"); h_cand_dT_shift ->SetTitle("Shifted-wire region candidates");
    h_cand_dT_sub   = (TH1D*)h_cand_dT->Clone("cand_dT_sub");   h_cand_dT_sub   ->SetTitle("Background-subtracted spectrum");
    
    float alphaQmax = 8e3;
    int   alphaQbins = 40;
    h_alpha_charge      = new TH1D("alpha_charge","Candidate alphas;Collection Plane Charge [e];Events", alphaQbins, 0, alphaQmax);
    h_alpha_charge_shift = (TH1D*)h_alpha_charge->Clone("alpha_charge_shift");
    h_alpha_charge_sub = (TH1D*)h_alpha_charge->Clone("alpha_charge_sub");
    h_alpha_charge_sub  ->SetTitle("Candidate alphas after background subtraction");
    
    float betaQmax = 80e3;
    int   betaQbins = 40;
    h_beta_charge      = new TH1D("beta_charge","Candidate betas;Collection Plane Charge [e];Events", betaQbins, 0, betaQmax);
    h_beta_charge_shift = (TH1D*)h_beta_charge->Clone("beta_charge_shift");
    h_beta_charge_sub = (TH1D*)h_beta_charge->Clone("beta_charge_sub");
    h_beta_charge_sub  ->SetTitle("Candidate betas after background subtraction");
    
    h_time_vs_rate_bipo = new TH1D("time_vs_rate_bipo","BiPo component of dT fit;Event time [hr];Rate per readout",timeBins,0,timeMax);
    h_time_vs_rate_bg   = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_rate_BG");
    h_time_vs_rate_bg   ->SetTitle("Background component");
   
    for(int i=0; i<timeBins;i++){
      h_poszy[i]        = new TH2D(Form("zy_%i",i),      "Bipo candidates;Z [cm]; Y [cm]", Zbins/10,  Zmin, Zmax, Ybins/10,  Ymin, Ymax);
      h_poszy[i]        ->SetOption("colz");
      h_poszy_shift[i]  = new TH2D(Form("zy_shift_%i",i),"Bipo BG;Z [cm]; Y [cm]",         Zbins/10,  Zmin, Zmax, Ybins/10,  Ymin, Ymax);
      h_poszy_shift[i]  ->SetOption("colz");
      h_poszy_sub[i]    = new TH2D(Form("zy_sub_%i",i),  "Post subtraction",               Zbins/10,  Zmin, Zmax, Ybins/10,  Ymin, Ymax);
      h_poszy_sub[i]    ->SetOption("colz");
    }

    if( _isMC ) {
      h_true_alpha_depne  = (TH1D*)h_alpha_charge->Clone("true_alpha_depne");
      h_true_alpha_depne  ->SetTitle("True ionization electrons from alpha");
      h_true_alpha_charge = (TH1D*)h_alpha_charge->Clone("true_alpha_charge");
      h_true_alpha_charge ->SetTitle("True alpha charge at anode");
      h_matched_alpha_charge = (TH1D*)h_alpha_charge->Clone("matched_alpha_charge");
      h_matched_alpha_charge ->SetTitle("Cluster charge matched to alpha");
    }

    // =====================================================
    // Diagnotic and utility histograms
    tdir_util->cd();

    if( fDoWireDiagnostics ) 
    h_nclusts_wire_ave= new TH1D("nclusts_perwire","Average clusters per wire;Average number of clusters per wire",1500,0,3); 
    h_clust_dT        = (TH1D*)h_cand_dT->Clone("clust_dT");  h_clust_dT->SetTitle("Same-wire cluster separations");
    
    h_2D_time_vs_dT      = new TH2D("2D_time_vs_dT",";Event time [hr];Time difference [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_2D_time_vs_dT_shift= (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_shift");
    h_2D_time_vs_dT_sub  = (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_sub");
    
    h_time_vs_N       = new TH1D("time_vs_N",";Event time [hr];Number of entries into dT plot",timeBins,0,timeMax);
    //h_time_vs_rate_bipo = new TH1D("time_vs_rate_bipo","BiPo component of dT fit;Event time [hr];Rate [sec^{-1}]",timeBins,0,timeMax);
    h_time_vs_ratio     = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_ratio");
    h_time_vs_ratio     ->SetTitle("Signal-to-background ratio");
    h_time_vs_ratio     ->GetYaxis()->SetTitle("Signal-to-background ratio");
 
    //timeBins    = 45; timeMax = 45;
    h_timeFine_vs_evts     = new TH1D("timeFine_vs_evts",";Event time [hr];Number of events",timeBins,0,timeMax);
    h_timeFine_vs_rate  = new TH1D("timeFine_vs_rate",";Event time [hr];Number of candidates per readout",timeBins,0,timeMax);
  
    //h_time_vs_ph    = new TH2D("time_vs_ph",";Event time [hr];Average track hit amplitude per evt [ADC]",timeBins,0,timeMax, 100,0,20);
    //h_time_vs_ph    ->SetOption("colz");
    //h_time_vs_ntrks = new TH2D("time_vs_ntrks",";Event time [hr];Number of tracks per evt",timeBins,0,timeMax, 80,0,80);
    //h_time_vs_ntrks    ->SetOption("colz");
    
  }
  
  
  //#################################################################################
  // Primary macro
  //#################################################################################
  void BiPo_macro()
  {
    // ******************************* 
    // Initial configurations
    // *******************************

    // open the file and set up the TTree
    std::cout<<"Reading input file "<<_fileName<<"\n";
    TFile* file = new TFile(_fileName.c_str(),"READ");
    fTree = (TTree*)file->Get("blipanaTrkMask/anatree");
  
    // set branches
    //fTree->SetBranchAddress("event",&event);                         
    //fTree->SetBranchAddress("run",&run);                             
    fTree->SetBranchAddress("timestamp",&timestamp);
    //fTree->SetBranchAddress("ntrks",&ntrks);
    //fTree->SetBranchAddress("trk_length",&trk_length);
    //fTree->SetBranchAddress("ave_trkhit_ph",&ave_trkhit_ph);
    fTree->SetBranchAddress("nclusts",&nclusts);                     
    //fTree->SetBranchAddress("clust_nwires",&clust_nwires);
    fTree->SetBranchAddress("clust_plane",&clust_plane);              
    fTree->SetBranchAddress("clust_wire",&clust_wire);                
    //fTree->SetBranchAddress("clust_startwire",&clust_startwire);                
    //fTree->SetBranchAddress("clust_endwire",&clust_endwire);                
    fTree->SetBranchAddress("clust_nhits",&clust_nhits);              
    fTree->SetBranchAddress("clust_charge",&clust_charge);            
    fTree->SetBranchAddress("clust_time",&clust_time);                
    //fTree->SetBranchAddress("clust_lhit_peakT",&clust_lhit_peakT);         
    //fTree->SetBranchAddress("clust_lhit_gof",&clust_lhit_gof);
    fTree->SetBranchAddress("clust_lhit_isfit",&clust_lhit_isfit);          
    //fTree->SetBranchAddress("clust_lhit_amp",&clust_lhit_amp);
    //fTree->SetBranchAddress("clust_lhit_rms",&clust_lhit_rms);
    fTree->SetBranchAddress("clust_blipid",&clust_blipid);            
    fTree->SetBranchAddress("nblips",&nblips);                       
    fTree->SetBranchAddress("blip_nplanes",&blip_nplanes);
    fTree->SetBranchAddress("blip_maxdiff",&blip_maxdiff);
    //fTree->SetBranchAddress("blip_x",&blip_x);                        
    fTree->SetBranchAddress("blip_y",&blip_y);                        
    fTree->SetBranchAddress("blip_z",&blip_z);                        
    fTree->SetBranchAddress("blip_charge",&blip_charge);              
    fTree->SetBranchAddress("blip_clustid_pl2",&blip_clustid_pl2);
    //fTree->SetBranchAddress("blip_clustid_pl1",&blip_clustid_pl1);
    //fTree->SetBranchAddress("blip_clustid_pl0",&blip_clustid_pl0);
    //fTree->SetBranchAddress("blip_clustid_pl1",&blip_clustid_pl1);
    //fTree->SetBranchAddress("blip_trkdist",&blip_trkdist);
    //fTree->SetBranchAddress("calib_lifetime",&calib_lifetime);
    //fTree->SetBranchAddress("calib_lifetime_err",&calib_lifetime_err);
    if( _isMC ) {
      fTree->SetBranchAddress("nedeps",&nedeps);
      fTree->SetBranchAddress("clust_edepid",&clust_edepid);            
      fTree->SetBranchAddress("nparticles",&nparticles);
      fTree->SetBranchAddress("isPrimary",&isPrimary);
      fTree->SetBranchAddress("trackID",&trackID);
      fTree->SetBranchAddress("pdg",&pdg);
      fTree->SetBranchAddress("startT",&startT);
      fTree->SetBranchAddress("depElectrons",&depElectrons);
      fTree->SetBranchAddress("numElectrons",&numElectrons);
      fTree->SetBranchAddress("edep_g4id",&edep_g4id);
      fTree->SetBranchAddress("edep_pdg",&edep_pdg);
      fTree->SetBranchAddress("edep_depne",&edep_depne);
      fTree->SetBranchAddress("edep_charge",&edep_charge);
      fTree->SetBranchAddress("edep_x",&edep_x);
      fTree->SetBranchAddress("edep_y",&edep_y);
      fTree->SetBranchAddress("edep_z",&edep_z);
    }

    // make output file to store plots
    std::string _outFileName = "output/plots_bipo_" + fFileName[fConfig];
    fOutFile = new TFile(_outFileName.c_str(), "recreate");
    tdir_plots  = fOutFile->mkdir("plots");
    tdir_util   = fOutFile->mkdir("util");
    
    // initialize all histograms
    makeHistograms();
    
    // for picky blip mode, set proper restrictions
    if( fPickyBlipMode ) {
      _betaMinPlanes    = 3;      // Min number of matched planes (must be 2 or 3)
      _betaMaxDiff      = 2;      // Difference in wire intersection points [cm]
    } 

    // Keep list of cluster counts per wire/channel
    std::map<int,int> map_wn;
    //std::map<int,int> map_chn;
    
    // Print-interval counter
    int print_counter = 0;

    for(int iWire=0; iWire<nWiresColl; iWire++){
      wireIsNoisy[iWire]  = false;
      wireIsBad[iWire]    = false;
  		if( std::find(fNoisyWires.begin(), fNoisyWires.end(), iWire ) != fNoisyWires.end() ) wireIsNoisy[iWire] = true;
  		if( std::find(fBadWires.begin(), fBadWires.end(), iWire ) != fBadWires.end() ) wireIsBad[iWire] = true;
    }
    

    // ****************************************************
    // Loop over the events
    // ****************************************************
    std::time_t loopStart = time(0);
    size_t totalEntries = fTree->GetEntries();
    for(size_t iEvent=0; iEvent < totalEntries; iEvent++){
      
      print_counter++;
      if( print_counter > 1000 ) {
        printf("========== EVENT %lu / %lu, BiPo count: %i =====================\n",iEvent,totalEntries,_numBiPo);
        print_counter = 1;
      }
      
      // ..... quick-test options ...........
      //int maxEvt    = 1000; if(  iEvent >= maxEvt ) break;
      //int sparsify  = 10; if(  (iEvent % sparsify) != 0 ) continue; 
      //..................................

      // Retrieve event info
      fTree->GetEntry(iEvent);
      _numEvents++;
      
      // Record the event time relative to start of dataset period
      double eventHour = ( timestamp - _t0 ) / 3600.;
      h_time_vs_N->Fill(eventHour);
      h_timeFine_vs_evts->Fill(eventHour);
      int eventTimeBin = int(eventHour / 5 );
      if( eventTimeBin >= timeBins ) continue;

      // Track-based diagonstics
      //h_time_vs_ntrks->Fill(eventHour,ntrks);
   
      std::vector<int> beta_g4ids;
      std::vector<int> alpha_g4ids;

      //std::vector<int> bipo_trackIDs;
      // Check truth info
      //std::cout<<"Event "<<event<<"\n";
      for(int i=0; i<nparticles; i++){
        
        if( isPrimary[i] && pdg[i] == 11 ) {
          beta_g4ids.push_back(trackID[i]);
        }

        if( isPrimary[i] && pdg[i] == 1000020040 ) {
          alpha_g4ids.push_back(trackID[i]);
          h_true_alpha_depne->Fill(depElectrons[i]);
          h_true_alpha_charge->Fill(numElectrons[i]);
          //std::cout<<numElectrons[i]<<"\n";
          //bipo_trackIDs.push_back(i);
        }
      
      }
      
      std::vector<int> alpha_edeps;
      std::vector<int> beta_edeps;
      
      for(int i=0; i<nedeps; i++){
        
        int g4index = FindG4Index(edep_g4id[i]);

        if( isPrimary[g4index] && edep_pdg[i] == 1000020040 ) {
          h_true_alpha_depne->Fill(edep_depne[i]);
          h_true_alpha_charge->Fill(edep_charge[i]);
          //float T = startT[g4index];
          //if( edep_charge[i] == 0 ) {
          //  std::cout<<"Found 0 charge edep... "<<i<<"   "<<edep_x[i]<<"  "<<edep_y[i]<<"   "<<edep_z[i]<<", T = "<<T/1000<<" us\n";
          //}
          //if( edep_depne[i] == 0 ) std::cout<<"blahblahn\n";
        }
      
      }
      
      std::vector<bool> clust_isAlpha(nclusts,false);
      std::vector<bool> clust_isBeta(nclusts,false);

      // look for collection plane clusts matched to an alpha
      for(int i=0; i < nclusts; i++){
        if( clust_plane[i] != 2 ) continue;
        int eid = clust_edepid[i];
        if( eid < 0 ) continue;
        int g4index = FindG4Index(edep_g4id[eid]);
        if( isPrimary[g4index] && pdg[g4index] == 1000020040 ) {
          h_matched_alpha_charge->Fill( clust_charge[i] );
          clust_isAlpha[i] = true;
        }
        if( isPrimary[g4index] && pdg[g4index] == 11 ) {
          clust_isBeta[i] = true;
        }
      }



      // ====================================================
      // Map of clust IDs per wire on collection plane
      // ====================================================
      _map_wire_clusters.clear();
      for(int i=0; i < nclusts; i++){
        if( clust_plane[i] != 2 ) continue;
        _map_wire_clusters[clust_wire[i]].push_back(i);
        map_wn[clust_wire[i]] += 1;
        //map_chn[clust_chan[i]] += 1;
        h_wt_clusts->Fill(clust_wire[i],clust_time[i]);
      }
    
      // ====================================================
      // Wire diagnostics (time-consuming!)
      // ====================================================
      if( fDoWireDiagnostics ) {

        for(int i=0; i<=nWiresColl; i++){ // plot cluster separations
          if( _map_wire_clusters.find(i) == _map_wire_clusters.end() ) continue;
          
          std::vector<float> v_t;
          for(int j=0; j<(int)_map_wire_clusters[i].size(); j++){
            int id = _map_wire_clusters[i].at(j);
            v_t.push_back(clust_time[id]*fSamplePeriod);
          }

          std::sort(v_t.begin(), v_t.end());
          for(int j=1; j<(int)v_t.size(); j++) 
            h_clust_dT->Fill(v_t.at(j)-v_t.at(j-1));
        
        }
      }
       
      // ===============================================================
      // Clear masks to track blip/cluster availability to avoid double-counts
      // ===============================================================
      _blipAvailable.assign(nblips, true);
      _clustAvailable.assign(nclusts, true);

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
        _blipAvailable[iBlip] = false;
     
        // find associated cluster on collection
        int   ic    = blip_clustid_pl2[iBlip];
        int   ic0    = blip_clustid_pl0[iBlip];
        int   ic1    = blip_clustid_pl1[iBlip];
  
        // plot ZY position
        h_zy_blips->Fill( blip_z[iBlip], blip_y[iBlip] );
        h_wt_blips->Fill( clust_wire[ic], clust_time[ic] );
  
        // skip if this cluster was already included in a BiPo candidate 
        if( !_clustAvailable[ic] ) continue;
        
        // skip pulse train hits
        //if( (fRequireHitFit && clust_lhit_gof[ic] < 0 )) continue;
        if( (fRequireHitFit && !clust_lhit_isfit[ic] )) continue;
      
  			// skip if this is a noisy wire
  		  //if( std::find(fNoisyWires.begin(), fNoisyWires.end(), clust_wire[ic] ) != fNoisyWires.end() ) continue;
  		  //if( std::find(fNoisyWires.begin(), fNoisyWires.end(), clust_startwire[ic] ) != fNoisyWires.end() ) continue;
  		  //if( std::find(fNoisyWires.begin(), fNoisyWires.end(), clust_endwire[ic] ) != fNoisyWires.end() ) continue;
        //
        if( wireIsNoisy[clust_wire[ic]] ) continue;

 
        // 3D match cuts for beta
        bool isWorthy = ( blip_nplanes[iBlip] >= _betaMinPlanes && blip_maxdiff[iBlip] < _betaMaxDiff );
       
        if((blip_nplanes[iBlip] == 3 && blip_maxdiff[iBlip] < 5 ) ) h_zy_blips_picky->Fill( blip_z[iBlip], blip_y[iBlip] );
        
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
  		  //if( std::find(fBadWires.begin(), fBadWires.end(), clust_wire[ic]+1 ) != fBadWires.end() ) continue;
  		  //if( std::find(fBadWires.begin(), fBadWires.end(), clust_wire[ic]-1 ) != fBadWires.end() ) continue;
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
        if( peakT < _minTick ) continue;
        if( peakT > _maxTick ) continue;
        



        //####################################################################
        // Old method using only the lead hit wire
        //####################################################################
        
        int fMaxCandidates = 1;

        // ------------------------------------------------------------
        // Begin search for alpha candidates
        // ------------------------------------------------------------
        int nclusts_inwindow        = 0;
        int nclusts_inwindow_shift  = 0;
        int nclusts_inwindow_shift2  = 0;
  
        // assign wire range
        int refwire = clust_wire[ic] + fRandomWireShift;
        if ( refwire < 0  || refwire > nWiresColl ) 
          refwire = clust_wire[ic] - fRandomWireShift;
        int w0  = refwire - fWireRange;
        int w1  = refwire + fWireRange;
        
        // -------------------------------------------
        // Search for standard candidates
        std::vector<BiPoCandidate> v_cands = FindCandidates(ic, w0, w1, false, nclusts_inwindow);
        h_nclusts_inwindow->Fill(nclusts_inwindow);
        h_ncands_inwindow->Fill((int)v_cands.size());
  
        // -------------------------------------------
        // Search for background candidates (wire shift)
        int shift = 2*fWireRange+1;
        
        // check that in shifted region, there aren't bad wires
        int ref_shift;
        int w0_shift, w1_shift;
        int n_regions_checked = 0;
        
        std::vector<BiPoCandidate> v_cands_shift;
        std::vector<BiPoCandidate> v_cands_shift2;
        
        // 1) positive shift
        ref_shift = refwire + shift;
        w0_shift = w0 + shift;
        w1_shift = w1 + shift;
  		  if( 
         //     std::find(fBadWires.begin(), fBadWires.end(), w0_shift ) == fBadWires.end() 
         // &&  std::find(fBadWires.begin(), fBadWires.end(), w1_shift ) == fBadWires.end() 
         // &&  std::find(fBadWires.begin(), fBadWires.end(), ref_shift ) == fBadWires.end() 
          !wireIsBad[w0_shift] && !wireIsBad[w1_shift] && !wireIsBad[ref_shift]
        ) {
          n_regions_checked++;
          v_cands_shift  = FindCandidates(ic, w0_shift, w1_shift, false, nclusts_inwindow_shift);
        }
        
        // 2) negative shift
        ref_shift = refwire - shift;
        w0_shift = w0 - shift;
        w1_shift = w1 - shift;
  		  if( 
          //    std::find(fBadWires.begin(), fBadWires.end(), w0_shift ) == fBadWires.end() 
          //&&  std::find(fBadWires.begin(), fBadWires.end(), w1_shift ) == fBadWires.end() 
          //&&  std::find(fBadWires.begin(), fBadWires.end(), ref_shift ) == fBadWires.end() 
          !wireIsBad[w0_shift] && !wireIsBad[w1_shift] && !wireIsBad[ref_shift]
        ) {
          n_regions_checked++;
          v_cands_shift2  = FindCandidates(ic, w0_shift, w1_shift, false, nclusts_inwindow_shift2);
        }

        // --------------------------------------------
        // Evaluate standard candidates
        // ---------------------------------------------
        if( v_cands.size()==1 && nclusts_inwindow <=  fMaxClustMult ) {
          
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
            h_timeFine_vs_rate  ->Fill(eventHour);
          }
          if( thisCand.dT >= fdT_min ) {
            _numBiPo++;
            if( clust_isBeta[ic] ) _numBiPo_true++;
            h_beta_charge         ->Fill(thisCand.q1);
            h_alpha_charge        ->Fill(thisCand.q2);
            h_cand_dT             ->Fill(thisCand.dT);
            h_2D_time_vs_dT       ->Fill(eventHour,thisCand.dT);
            h_poszy[eventTimeBin] ->Fill(blip_z[thisCand.blipID],blip_y[thisCand.blipID]);
          }

        }//end evaluation of standard cands
        
  
  
        // --------------------------------------------
        // Evaluate dT-shift candidates
        // ---------------------------------------------
        float weight = 1./float(n_regions_checked);
        //std::cout<<"applying weight "<<weight<<" to background\n";
        if( v_cands_shift.size() && v_cands_shift.size() <= fMaxCandidates && nclusts_inwindow_shift <=  fMaxClustMult ) {
          for(auto& thisCand : v_cands_shift ) {
            if( thisCand.dT >= fdT_min ) {
              h_beta_charge_shift         ->Fill(thisCand.q1,weight);
              h_alpha_charge_shift        ->Fill(thisCand.q2,weight);
              h_cand_dT_shift             ->Fill(thisCand.dT, weight);
              h_2D_time_vs_dT_shift       ->Fill(eventHour,thisCand.dT, weight);
              h_poszy_shift[eventTimeBin] ->Fill(blip_z[thisCand.blipID],blip_y[thisCand.blipID]);
            }
          }//endloop over candidates
        }//end evaluation of wire-shift candidates
        
        if( v_cands_shift2.size() && v_cands_shift2.size() <= fMaxCandidates && nclusts_inwindow_shift2 <=  fMaxClustMult ) {
          for(auto& thisCand : v_cands_shift2 ) {
            if( thisCand.dT >= fdT_min ) {
              h_beta_charge_shift       ->Fill(thisCand.q1,weight);
              h_alpha_charge_shift      ->Fill(thisCand.q2,weight);
              h_cand_dT_shift           ->Fill(thisCand.dT, weight);
              h_2D_time_vs_dT_shift     ->Fill(eventHour,thisCand.dT, weight);
              h_poszy_sub[eventTimeBin] ->Fill(blip_z[thisCand.blipID],blip_y[thisCand.blipID]);
            }
          }//endloop over candidates
        }//end evaluation of wire-shift candidates
      
        //#######################################################################################################


        // ------------------------------------------------------------
        // Begin search for alpha candidates
        // ------------------------------------------------------------
      //  int nclusts_inwindow        = 0;
      //  int nclusts_inwindow_shift  = 0;
      //  int nclusts_inwindow_shift2 = 0;
      //  
      //  int nwiresChecked;
      //  std::vector<BiPoCandidate> v_cands = FindCandidatesCluster(ic, false, nclusts_inwindow, nwiresChecked);

      //  int nwiresChecked_bg;
      //  std::vector<BiPoCandidate> v_cands_shift = FindCandidatesCluster(ic, true, nclusts_inwindow_shift, nwiresChecked_bg);

      //  h_nclusts_inwindow->Fill(nclusts_inwindow);
      //  h_ncands_inwindow->Fill(v_cands.size());

      //  // --------------------------------------------
      //  // Evaluate standard candidates
      //  // ---------------------------------------------
      //  if( v_cands.size() == 1 && nclusts_inwindow <= fMaxClustMult ) {
      //    // plot locations
      //    h_zy_bipos->Fill( blip_z[iBlip], blip_y[iBlip] );
      //    h_wt_bipos->Fill( clust_wire[ic], clust_time[ic] );
      //    // loop through and update histograms/counts
      //    BiPoCandidate thisCand = v_cands.at(0);
      //      _clustAvailable[thisCand.id1] = false;
      //      _clustAvailable[thisCand.id2] = false;
      //      int blipid = clust_blipid[thisCand.id2];
      //      if( blipid >= 0 ) _blipAvailable[blipid] = false;
      //      if( thisCand.dT > 6.25 && thisCand.dT < 312.5 ) {
      //        _numBiPo_6_312++;
      //        h_timeFine_vs_rate  ->Fill(eventHour);
      //      }
      //      if( thisCand.dT > fdT_min ) {
      //        _numBiPo++;
      //        h_beta_charge       ->Fill(thisCand.q1);
      //        h_alpha_charge      ->Fill(thisCand.q2);
      //        h_alpha_nhits       ->Fill(clust_nhits[thisCand.id2]);
      //        h_cand_dT           ->Fill(thisCand.dT);
      //        h_2D_time_vs_dT        ->Fill(eventHour,thisCand.dT);
      //        //if( thisCand.dT > 10 && thisCand.dT < 20 ) {
      //        //  std::cout
      //        //  <<"Found a weird one! dT "<<thisCand.dT
      //        //  <<", RMS " <<clust_lhit_rms[thisCand.id1]<<", "<<clust_lhit_rms[thisCand.id2]
      //        //  <<", amp: "<<clust_lhit_amp[thisCand.id1]<<", "<<clust_lhit_amp[thisCand.id2]<<"\n";
      //        //}
      //      }
      //  }//end evaluation of standard cands
  
      // 
      //  // --------------------------------------------
      //  // Evaluate dT-shift candidates
      //  // ---------------------------------------------
      //  float weight = float(nwiresChecked)/float(nwiresChecked_bg);
      //  //std::cout<<"Applying weight "<<weight<<"\n";
      //  if( v_cands_shift.size() == 1 && nclusts_inwindow_shift <= fMaxClustMult ) {
      //    BiPoCandidate thisCand = v_cands_shift.at(0);
      //      if( thisCand.dT >= fdT_min ) {
      //        h_beta_charge_shift  ->Fill(thisCand.q1,weight);
      //        h_alpha_charge_shift ->Fill(thisCand.q2,weight);
      //        h_cand_dT_shift      ->Fill(thisCand.dT,weight);
      //        h_2D_time_vs_dT_shift ->Fill(eventHour,thisCand.dT,weight);
      //      }
      //  }//end evaluation of dT-shift candidates




      }//end loop over 3D blips

    }//endloop over events
    double loopDuration = ( time(NULL) - loopStart );

    
    // ***************************************************
    // Scale dT plots so they're per readout
    // ***************************************************
    _totalLiveTime = float(_numEvents) * _liveTimePerEvt;
    float scaleFact = 1./( float(_numEvents) );
    h_clust_dT        ->Scale( scaleFact );
    h_cand_dT         ->Scale( scaleFact );
    h_cand_dT_shift   ->Scale( scaleFact );
    h_timeFine_vs_rate->Divide( h_timeFine_vs_evts );

    
    // ***************************************************
    // Histogram subtraction time!
    // ***************************************************
    h_cand_dT_sub     ->Add(h_cand_dT,      1);
    h_2D_time_vs_dT_sub  ->Add(h_2D_time_vs_dT,   1);
    h_alpha_charge_sub->Add(h_alpha_charge, 1);
    h_beta_charge_sub->Add(h_beta_charge, 1);

  
    h_cand_dT_sub     ->Add(h_cand_dT_shift,     -1.);
    h_2D_time_vs_dT_sub  ->Add(h_2D_time_vs_dT_shift,  -1.);
    h_alpha_charge_sub->Add(h_alpha_charge_shift,-1.);
    h_beta_charge_sub->Add(h_beta_charge_shift,-1.);
    
    for(int i=0; i<timeBins; i++){
      h_poszy_sub[i]->Add(h_poszy[i],1);
      h_poszy_sub[i]->Add(h_poszy_shift[i],-1);
    }

   
    // ***************************************************
    // Check for noisy wires if doing wire diagnostics
    // ***************************************************
    float noiseThresh = 0.25;
    if( fDoWireDiagnostics ) {
      std::cout<<"Noisy wires (>"<<noiseThresh<<" cands/evt): \n";
      int nNoisy = 0; 
      // populate histogram of average cluster multiplicity per wire
      if( map_wn.size() ) {
        for(int iwire=0; iwire<=nWiresColl; iwire++){
          float a = 0;
          if( map_wn.find(iwire) != map_wn.end() ) {
            a = map_wn.at(iwire) / (float)_numEvents;
            if( a > noiseThresh ) { std::cout<<iwire<<", "; nNoisy++; } 
          }
          h_nclusts_wire_ave->Fill(a);
        }
        std::cout<<"\n--> Found "<<nNoisy<<" wires with ave clusters/evt > "<<noiseThresh<<"\n";
      
        std::cout<<"Checking for bad/missing wires\n";
        // check for bad/missing wires
        int nBad = 0; 
        for(int iwire=0; iwire<=nWiresColl; iwire++){
          if( map_wn.find(iwire) == map_wn.end() ) {
            std::cout<<iwire<<", "; nBad++;
          }
        }
        std::cout<<"\n--> Found "<<nBad<<" wires with no activity\n";
      }
      
    }
  

    // ***************************************************
    // Write all histos currently in stack
    // ***************************************************
    fOutFile->Write(); 
 

    // ***************************************************
    // Make plots
    // ***************************************************
    makePlots();
  

    // ***************************************************
    // Output summary
    // ***************************************************
    printf("\n*******************************************\n");
    printf("File                : %s\n",      fFileName[fConfig].c_str()); 
    printf("Total events        : %i\n",      _numEvents); 
    printf("Total live time     : %f sec\n",  _totalLiveTime);
    printf("Live time per evt   : %f us\n",   _liveTimePerEvt*1e6);
    printf("dT min/max          : %.2f-%.2f us\n",  fdT_min,fdT_max);
    printf("Ave cands / evt     : %f (%f in 6.25-312.us)\n",h_cand_dT->GetEntries()/(float)_numEvents, _numBiPo_6_312/(float)_numEvents );
    printf("  - truth-matched   : %f (%f in 6.25-312.us)\n",_numBiPo_true/(float)_numEvents, _numBiPo_6_312_true/(float)_numEvents );
    printf("Processing time     : %f sec/evt\n", loopDuration/float(_numEvents));
    printf("Excluding %i noisy wires \n",     (int)fNoisyWires.size());	
    printf("*******************************************\n\n");
    fOutFile->Close();
  
  }
  
  
  
  //#################################################################################
  // Make plots here
  //#################################################################################
  void makePlots()
  {
   
   
    if( !_isMC ) {
  
    // ============================================
    // Rate vs time plots
    // ============================================
    std::cout<<"\nMaking rate vs time plot...\n";
   
    
    h_time_vs_rate_bipo->SetOption("E1");

    TGraphErrors* gr_sig = new TGraphErrors();
    //TGraphErrors* gr_bg;

    TH1D* h_slice;
    for(int i=1; i<=(int)h_time_vs_N->GetXaxis()->GetNbins(); i++){
      h_slice = Make1DSlice( h_2D_time_vs_dT_sub, i, i, Form("dTfit_%i",i) );
      h_slice->Scale( 1. / (float)h_time_vs_N->GetBinContent(i) );
      float t = h_time_vs_N->GetXaxis()->GetBinCenter(i);
      float dt = h_time_vs_N->GetXaxis()->GetBinWidth(i)/sqrt(12);
      std::cout<<"time: "<<t<<" hrs\n";

      FitResult fr = fitdT(h_slice,false);
      if( fr.rate_signal != -9 ) {
        int npt = gr_sig->GetN();
        gr_sig->SetPoint( npt , t, fr.rate_signal );
        gr_sig->SetPointError(  npt,  dt, fr.rate_signal_err);
        h_time_vs_rate_bipo ->SetBinContent(i,fr.rate_signal);
        h_time_vs_rate_bipo ->SetBinError(i, fr.rate_signal_err);
        h_time_vs_rate_bg   ->SetBinContent(i,fr.rate_bg);
        h_time_vs_rate_bg   ->SetBinError(i, fr.rate_bg_err);
        h_time_vs_ratio     ->SetBinContent(i,fr.ratio);
      }
    }
    
    fOutFile->cd();
    h_time_vs_rate_bipo ->Write(0, TObject::kOverwrite );
    h_time_vs_rate_bg   ->Write(0, TObject::kOverwrite );
    h_time_vs_ratio     ->Write(0, TObject::kOverwrite );
  
    TH1D* h1  = h_time_vs_rate_bipo; //(TH1D*)h_time_vs_rate_bipo->Clone("bipo");
    TH1D* h2  = (TH1D*)h_time_vs_rate_bg->Clone("bg");
    FormatTH1D(h1, kBlue, 1, 3 );
    FormatTH1D(h2, kRed, 1, 3 );
    //FormatTH1D(h3, kBlack, 1, 3 );
   
    std::string name = "c_time_vs_rate";
    TCanvas* c = new TCanvas(name.c_str(),name.c_str(),600,450);
    gPad->SetMargin(mar_l, mar_r, mar_b, mar_t ); 
    gPad->SetGridy(1);
    //float max = std::max( GetHistMax(h1), GetHistMax(h2) );
    float min = std::min(GetHistMin(h1),GetHistMin(h2));
    float max = std::max(GetHistMax(h1),GetHistMax(h2));
    float range = (max-min);
    float center = (max+min)/2.;
    float yl = center - 1.2*range/2.;
    float yh = center + 1.2*range/2.;
    h1->GetYaxis()->SetRangeUser(min*1.2,max*1.2);
    h1->DrawCopy();
    h2->DrawCopy("same");
    //h3->DrawCopy("same");
    tdir_plots->cd();
    c->Write();
   
    /*
    name = "c_time_vs_rate_gr";
    TCanvas* c2 = new TCanvas(name.c_str(),name.c_str(),600,450);
    gPad->SetMargin(mar_l, mar_r, mar_b, mar_t ); 
    gr_sig->Draw("AP");
    c2->Write();
    */
    
    }
    

    // ============================================
    // Do final fit on dT spectrum
    // ============================================
    fitdT( h_cand_dT_sub, true );

  } 
  
  
  
  //################################################################################
  // Function that performs the dT fit
  //#################################################################################
  FitResult fitdT(TH1D* h, bool writeCanvas = false ){
  
    FitResult out;
  
    std::cout<<"Fitting dT spectrum "<<h->GetTitle()<<", "<<h->GetEntries()<<"\n";
    std::string label = h->GetName();
    TCanvas* c = new TCanvas(Form("c_fit_%s",label.c_str()),Form("c_fit_%s",label.c_str()),600,500);
    gPad->SetMargin(mar_l, mar_r, mar_b, mar_t ); 
    TH1D* hc = (TH1D*)h->Clone();
    float histMax = GetHistMax(hc);
   
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
    fit->SetParameter(1, histMax );
    fit->SetParLimits(1, 0, 2*histMax );
    fit->FixParameter(2, 164.3 );

    // Draw plot and fit
    c->cd();
    hc->Fit(fit,"QR"); // "WL" = log likelihood w/ weighted bins 
    hc->DrawCopy();
    fit->DrawCopy("same");
    gStyle->SetOptFit(1112);
    //gPad->SetLogy();
    if( writeCanvas ) {
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
    float A       = fit->GetParameter(1);
    float dA      = fit->GetParError(1);
    float B       = fit->GetParameter(2);
    float dB      = fit->GetParError(2);
    float n_bipo  = (A*B)/fdT_binSize;
    float n_bipo_err = sqrt( pow(dA/A,2)+pow(dB/B,2) )*n_bipo;

    // Components within the dT selection window
    float N_total = hc      ->Integral(0,hc->GetXaxis()->GetNbins());
    float N_fit   = fit     ->Integral(fdT_min,fdT_max)/fdT_binSize;
    float N_flat  = f_flat  ->Integral(fdT_min,fdT_max)/fdT_binSize;
    //float N_expbg = f_expBG ->Integral(fdT_min,fdT_max)/fdT_binSize;
    float N_bipo  = f_bipo  ->Integral(fdT_min,fdT_max)/fdT_binSize;
    float N_bg    = N_flat; //N_expbg;
    float N_bg_err = fit->IntegralError(fdT_min,fdT_max)/fdT_binSize;
    float sbratio = N_bipo / N_bg;

    float activity_mBq = 1e3*( (n_bipo*_fidCorFactor) / _liveTimePerEvt ) / 85000.;

    
    printf("================ dT fit =================\n");
    printf("p0                  : %f +/- %f\n", fit->GetParameter(0),fit->GetParError(0));
    printf("p1                  : %f +/- %f\n", fit->GetParameter(1),fit->GetParError(1));
    printf("p2                  : %f +/- %f\n", fit->GetParameter(2),fit->GetParError(2));
    printf("Chi2/ndf            : %f\n",        fit->GetChisquare()/fit->GetNDF());
    printf("Total entries       : %f\n",        hc->GetEntries());
    printf("Total rate          : %f per evt\n",N_total);
    printf(" - BiPo             : %f per evt\n",N_bipo);
    printf(" - Flat             : %f per evt\n",N_flat);
    printf("S/BG ratio          : %f \n",sbratio);
    printf("Extrapolated rate   : %f BiPos per evt\n",n_bipo);
    printf("Fiducial correction : %f BiPos per evt\n",n_bipo*_fidCorFactor);
    printf("Equivalent activity : %f mBq/kg\n",       activity_mBq);
  
    out.rate_signal       = n_bipo;
    out.rate_signal_err   = n_bipo_err;
    out.rate_bg           = N_bg;
    out.rate_bg_err       = N_bg_err;
    out.ratio             = sbratio;
    return out;
  
  }
  

  //################################################################################
  // Function to search for alpha candidates relative to a beta candidate
  //#################################################################################
  std::vector<BiPoCandidate> FindCandidates(int ic, int start, int end, bool flipdT, int& npileup ) {
        
    std::vector<BiPoCandidate> v;
    npileup = 0;
    
    for(int iWire = start; iWire <= end; iWire++){
    
      for(auto& jc : _map_wire_clusters[iWire] ) {
        if( ic == jc ) continue;
        if( !_clustAvailable[jc] ) continue;
        //if( fRequireHitFit && clust_lhit_gof[jc] < 0 ) continue; //MinHitGOF!clust_lhit_isfit[jc]) continue;
        if( fRequireHitFit && !clust_lhit_isfit[jc] ) continue; //MinHitGOF!clust_lhit_isfit[jc]) continue;
        
        //... temporary ... 
        //exclude electrons
        int eid = clust_edepid[jc];
        if( eid >= 0 ) {
          //std::cout<<"Found an alpha cand that is matched "<<eid<<"   pdg "<<edep_pdg[eid]<<"\n";
          if( fabs(edep_pdg[eid]) == 11 ) continue;
        }
        float dT  = (clust_time[jc]-clust_time[ic])*fSamplePeriod;
       
        // 3D plane-match requirement for alpha
        if( fAlphaReq3D && clust_blipid[jc] < 0 ) continue;
        
        // Count clusters in forward window, and fill some pre-cut histos
        if (flipdT) dT *= -1.;
  
        if( dT > 0 && fabs(dT) < fdT_max ) {
          npileup++;
        
            // --- alpha charge/nhits cut ---
            if(   clust_charge[jc] > fAlphaCharge_min 
              &&  clust_charge[jc] < fAlphaCharge_max 
              &&  clust_nhits[jc] <= fAlphaHits_max ) {
              
              int blipID = clust_blipid[ic];
              BiPoCandidate c = { blipID, ic, jc, fabs(dT), clust_charge[ic], clust_charge[jc]};
              v.push_back(c);
            }
          
         
        } 
  
      }//<-- end loop over clusters on this wire
    }//endloop over wires
  
    return v;
  
  }
  

  //################################################################################
  std::vector<BiPoCandidate> FindCandidatesCluster(int ic, bool useSideBands, int& npileup, int& nwiresChecked ) {
        
    // identify which wires to check
    std::set<int> wires;
    
    int w_start = clust_startwire[ic] + fRandomWireShift;
    int w_end   = clust_endwire[ic] + fRandomWireShift;
    if ( w_start < 0  || w_start > nWiresColl ) w_start -= 2*fRandomWireShift;
    if ( w_end < 0  || w_end > nWiresColl ) w_end -= 2*fRandomWireShift;
    
    if( useSideBands ) {
      w_start -= (2*fWireRange+1);
      w_end   += (2*fWireRange+1);
    }

    int w1 = w_start-fWireRange;
    int w2 = w_end+fWireRange;
    for(int w=w1; w<=w2; w++){
      if( abs(w-w_start)<= fWireRange ) wires.insert(w);
      if( abs(w-w_end)  <= fWireRange ) wires.insert(w);
    }
    

    //std::cout<<"Cluster has "<<clust_nwires[ic]<<" --> BG? "<<useSideBands<<"  identified "<<wires.size()<<" wires to check for alpha on\n";
    
    nwiresChecked = wires.size();

    std::vector<BiPoCandidate> v;
    npileup = 0;
   
    for(auto iWire : wires ) {
      //std::cout<<iWire<<"\n";
      for(auto& jc : _map_wire_clusters[iWire] ) {
        if( ic == jc ) continue;
        if( !_clustAvailable[jc] ) continue;
        //if( fRequireHitFit && clust_lhit_gof[jc]<0 ) continue;
        if( fRequireHitFit && !clust_lhit_isfit[jc] ) continue;
        if( fAlphaReq3D && clust_blipid[jc] < 0 ) continue;
        
        // TODO: correct this for wire offset
        float dT  = (clust_time[jc]-clust_time[ic])*fSamplePeriod;
        
        if( dT > 0 && fabs(dT) < fdT_max ) {
          npileup++;
        
            // --- alpha charge/nhits cut ---
            if(   clust_charge[jc] > fAlphaCharge_min 
              &&  clust_charge[jc] < fAlphaCharge_max 
              &&  clust_nhits[jc] <= fAlphaHits_max ) {
          
              int blipID = clust_blipid[ic];
              BiPoCandidate c = { blipID, ic, jc, fabs(dT), clust_charge[ic], clust_charge[jc]};
              v.push_back(c);
            }
        } 
      }
    }//end check on start wire
      
    return v;
  
  }


  //###################################################################
  int FindG4Index(int g4id) {
    for(size_t i=0; i<nparticles; i++) {
      if( trackID[i] == g4id ) return i;
    }
    return -9;
  }

