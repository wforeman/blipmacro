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
  float fMinHr    = -40e9; 
  float fMaxHr    = 40e9; 
  float qscale    = 1.0; //0.905;

  // --- Input files ---
  infile_t fInputFiles[4] = {
  /*0*/   { "BlipAna_20230108_BiPo_OverlayNEST_8ms.root",    "blipanaTrkMask/anatree", true,   0,0},
  /*0*/  //{ "BlipAna_20230108_BiPo_Overlay.root",         "blipanaTrkMask/anatree", true,   0,0},
  /*0*/  //{ "BlipAna_20221215_BiPo_Overlay.root",         "blipanaTrkMask/anatree", true,   0,0},

  /*1*/  { "BlipAna_20230108_Data_RadonDoping_FullFilter.root",   "blipanaTrkMask/anatree", false,  1627415210, 1627592728},
  /*2*/  { "BlipAna_20230108_Data_RadonDoping_FilterBypass.root", "blipanaTrkMask/anatree", false,  1627594380, 1627761265},
  /*3*/  { "BlipAna_20230108_Data_Run3_Unbiased.root",            "blipanaTrkMask/anatree", false,  1528526500, 1532448800}
  };

  // tighter cuts to try for Chao:
  //  - pickybeta
  //  - alphaE < 0.15 MeVee
  //  - dT range: 60-340us

  // --- General selection options ---
  bool  fFidVolCut        = 1;      // Fiducialize beta
  int   fBetaMinPlanes    = 2;      // Min matched planes (3 planes == "picky")
  int   fBetaWires_max    = 4;      // Max wires in beta collPlane cluster
  float fBetaEnergy_min   = 0.5;    // 
  float fBetaEnergy_max   = 3.5;    
  int   fWireRange        = 0;      // +/- range to look for alpha candidate;
  int   fAlphaWires_max   = 2;
  float fAlphaEnergy_min  = 0.00;
  float fAlphaEnergy_max  = 0.24; //0.24; 
  bool  fLinearizeCorr    = 1;      // Linearize detector effects control factor
  float fdT_binSize       = 20.;    // Bin width for all dT spectra plots [us]
  float fdT_min           = 20.;    // Min dT for looking for candidate [us]
  float fdT_max           = 500.;   // Max dT for looking for candidate [us]
  bool  fSkipNoisyWires   = 1;

  // --- MC efficiency for equiv activity calc ---
  double fEfficiencyMC      = 0.078;
  double fEfficiencyMC_err  = 0.030;

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
  
  // --- Noisy wires to skip (collection plane) ---
  /*
  std::vector<int> fNoisyWires{
  0, 1, 2, 3, 6, 7, 83, 384, 1151, 1247, 1535, 1540, 1919, 1920, 2303, 2311, 2334, 2335, 
  2400, 2415, 2687, 2733, 2753, 2783, 2879, 3071, 3072, 3215, 3263, 3274, 3286, 3299, 
  3318, 3327, 3385, 3391, 3408, 3409, 3410, 3411, 3412, 3413, 3414, 3415, 3416, 3417, 
  3418, 3419, 3420, 3421, 3422, 3423, 3424, 3425, 3426, 3427, 3428, 3429, 3430, 3431, 
  3432, 3433, 3434, 3435, 3436, 3437, 3438, 3439, 3440, 3441, 3442, 3443, 3444, 3445, 
  3446, 3447, 3448, 3449, 3451, 3452, 3453, 3454, 3455
  };
  */
  
  std::vector<int> fNoisyWires{
    374, 375, 376, 377, 378, 378, 379, 
    760, 761, 762, 763, 
    1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153,
    1284, 1285, 1286,
    1526, 1527, 1528, 1529, 1530, 1531, 1532, 1533, 1534, 1535,
    1912, 1913, 1914, 1915, 1916, 1917, 1918, 1919, 1920, 1921,
    2111,
    2224, 2225, 2226, 2227, 2228, 2229, 2230,
    2296, 2297, 2298, 2299, 2300, 2301, 2302, 2303, 2304,
    2680, 2681, 2683, 2683, 2684, 2685, 2686, 2687,
    2732, 2733, 2734, 2735, 2781, 2782, 2783,
    2823,2834,  2847, 2851, 
    2974, 2996, 3003, 3041,
    3052, 3070, 3071,
    3215
  };

  


  //#######################################################################
  // Derived parameters
  //######################################################################
  // Fiducial vol correction factor
  float dz = fZlim[1]-fZlim[0];
  float dy = fYlim[1]-fYlim[0];
  float _fiducialFrac = (fFidVolCut) ? std::min(1., (dz*dy)/(1036.*232.) ) : 1.0;
  // Counters / maps / etc
  bool  _isMC           = false;
  int   _numEvents      = 0;
  int   _numBiPo        = 0;
  int   _numBiPo_mcmatch   = 0;
  int   _numBiPo_true     = 0;
  int   _numBiPo_true_perfectReco = 0;
  std::vector<bool> _clustAvailable;
  std::map<int,std::vector<int>> _map_wire_clusters;
  std::vector<bool> wireIsNoisy (nWiresColl,false);

  // Live time
  int   _minTick         = 0;
  int   _maxTick         = 6400 - (int)fdT_max*2;
  float _liveTimePerEvt  = _maxTick*fSamplePeriod*1e-6; //sec
  float _totalLiveTime  = 0;
      
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
  std::vector<BiPoCandidate>  FindCandidates(int, int, int, bool, int&, int&);
  FitResult                   fitdT(TH1D*,bool,bool);
  
  // ROOT objects
  TTree*      fTree;
  TFile*      fOutFile;

  // Histograms
  TDirectory* tdir_util;
  TDirectory* tdir_plots;
  TDirectory* tdir_truth;
 
  TH1D* h_cuts;

  TH1D* h_nclusts_perwire;

  TH1D* h_nclusts_inwindow;
  TH1D* h_ncands_inwindow;
  TH2D* h_wt_clusts;
  TH2D* h_wt_blips;
  TH2D* h_wt_blips_filt;
  TH2D* h_wt_bipos;
  TH2D* h_zy_bipos;
  TH2D* h_zy_bipos_bg;
  TH2D* h_zy_bipos_sub;
  
  TH1D* h_cand_dT;
  TH1D* h_cand_dT_bg;
  TH1D* h_cand_dT_sub;
  
  TH1D* h_control_dT;
  TH1D* h_control_dT_bg;
  TH1D* h_control_dT_ratio;
  TH1D* h_control_ratio;
  
  TH1D* h_time_vs_rate; 
  TH1D* h_time_vs_rate_bg;
  TH1D* h_time_vs_activity; 

  TH1D* h_beta_charge;
  TH1D* h_beta_charge_bg;
  TH1D* h_beta_charge_sub;
  TH1D* h_beta_energy;
  TH1D* h_beta_energy_bg;
  TH1D* h_beta_energy_sub;
  TH1D* h_beta_amp[3];
  TH1D* h_alpha_charge;
  TH1D* h_alpha_charge_bg;
  TH1D* h_alpha_charge_sub;
  TH1D* h_alpha_energy;
  TH1D* h_alpha_energy_bg;
  TH1D* h_alpha_energy_sub;
  TH1D* h_alpha_amp[3];

  TH1D* h_true_alpha_depne;
  TH1D* h_true_alpha_charge;
  TH1D* h_matched_alpha_charge;
  TH1D* h_beta_trueEnergy;
  TH1D* h_beta_trueEnergySum;
  TH1D* h_beta_trueEnergySum_reco;
  TH1D* h_beta_trueEnergySum_recoCuts;
  TH1D* h_beta_nwires;
  TH1D* h_alpha_nwires;
  
  TH1D* h_time_vs_N;        
  TH2D* h_2D_time_vs_dT;
  TH2D* h_2D_time_vs_dT_bg;
  TH2D* h_alpha_energyVsdT;
  
  //##########################################################################
  // Initialize histograms
  //##########################################################################
  void makeHistograms()
  {
    fOutFile->cd();
    tdir_plots  = fOutFile->mkdir("plots");
    tdir_util   = fOutFile->mkdir("util");
    tdir_truth  = fOutFile->mkdir("truth");

    
    h_nclusts_perwire = new TH1D("nclusts_perwire","Collection plane;Wire number",3456,0,3456);
    h_nclusts_inwindow  = new TH1D("nclusts_inwindow","Mean clusters per wire in time window following Bi-candidate",50,0,5);
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
    h_cand_dT_bg  = (TH1D*)h_cand_dT->Clone("cand_dT_bg");  h_cand_dT_bg  ->SetTitle("Opposite dT candidates");
    h_cand_dT_sub   = (TH1D*)h_cand_dT->Clone("cand_dT_sub");   h_cand_dT_sub   ->SetTitle("Background-subtracted spectrum");
   
    h_control_dT        = new TH1D("control_dT","OFFSET REGION;Time difference [#mus];Number of candidates", dTbins,0.,fdT_max);
    h_control_dT_bg   = new TH1D("control_dT_bg","OFFSET REGION;Time difference [#mus];Number of candidates", dTbins,0.,fdT_max);
    h_control_dT_ratio  = new TH1D("control_dT_ratio","OFFSET REGION;Time difference [#mus];Number of candidates", dTbins,0.,fdT_max);
    
    float alphaQmax   = 6e3;
    int   alphaQbins  = alphaQmax/200.;
    float betaQmax    = 90e3;
    int   betaQbins   = betaQmax/2000;
    float alphaEmax   = fAlphaEnergy_max;
    int   alphaEbins  = fAlphaEnergy_max/0.01;
    float betaEmax    = fBetaEnergy_max;
    int   betaEbins   = fBetaEnergy_max/0.10;
    
    h_beta_energy       = new TH1D("beta_energy","Candidate betas;Energy [MeV];Entries per second", betaEbins, 0, betaEmax);
    h_beta_energy_bg    = (TH1D*)h_beta_energy->Clone("beta_energy_bg");
    h_beta_energy_sub   = (TH1D*)h_beta_energy->Clone("beta_energy_sub");
    h_alpha_energy      = new TH1D("alpha_energy","Candidate alphas;Electron-equivalent energy [MeVee];Entries per second", alphaEbins, 0, alphaEmax);
    h_alpha_energy_bg   = (TH1D*)h_alpha_energy->Clone("alpha_energy_bg");
    h_alpha_energy_sub  = (TH1D*)h_alpha_energy->Clone("alpha_energy_sub");
    //h_alpha_energy_sub  ->SetTitle("Candidate alphas after background subtraction");
    //h_beta_energy_sub   ->SetTitle("Candidate betas after background subtraction");

    h_beta_charge       = new TH1D("beta_charge","Candidate betas;Collected charge [e^{-}];Events", betaQbins, 0, betaQmax);
    h_beta_charge_bg    = (TH1D*)h_beta_charge->Clone("beta_charge_bg");
    h_beta_charge_sub   = (TH1D*)h_beta_charge->Clone("beta_charge_sub");
    h_alpha_charge      = new TH1D("alpha_charge","Candidate alphas;Collected charge [e^{-}];Entries per second", alphaQbins, 0, alphaQmax);
    h_alpha_charge_bg   = (TH1D*)h_alpha_charge->Clone("alpha_charge_bg");
    h_alpha_charge_sub  = (TH1D*)h_alpha_charge->Clone("alpha_charge_sub");
    //h_alpha_charge_sub  ->SetTitle("Candidate alphas after background subtraction");
    //h_beta_charge_sub   ->SetTitle("Candidate betas after background subtraction");

    h_beta_amp[0]   = new TH1D("beta_amp","Candidate betas;Hit amplitude [ADC];Entries per second", 40,0,20);
    h_beta_amp[1]   = (TH1D*)h_beta_amp[0]->Clone("beta_amp_bg");
    h_beta_amp[2]   = (TH1D*)h_beta_amp[0]->Clone("beta_amp_sub");
    h_alpha_amp[0]  = new TH1D("alpha_amp","Candidate alphas;Hit amplitude [ADC];Entries per second", 50,0,5);
    h_alpha_amp[1]  = (TH1D*)h_alpha_amp[0]->Clone("alpha_amp_bg");
    h_alpha_amp[2]  = (TH1D*)h_alpha_amp[0]->Clone("alpha_amp_sub");

    

    h_alpha_energyVsdT = new TH2D("alpha_energyVsdT","Alpha candidates;#DeltaT [#mus];Energy [MeVee]",dTbins,0.,fdT_max,alphaEbins,0,alphaEmax);
    h_alpha_energyVsdT->SetOption("colz");

    const float binPeriodHrs[4] = { 2,  2,   2,   48 };
    const float binPeriodMax[4] = { 44, 44, 44,  1104 };
    int   timeBins  = binPeriodMax[fConfig]/binPeriodHrs[fConfig];
    float timeMax   = binPeriodMax[fConfig];
    h_time_vs_rate      = new TH1D("time_vs_rate",";Time [hr];Rate per 3.2 ms readout",timeBins,0,timeMax);
    h_time_vs_activity  = (TH1D*)h_time_vs_rate->Clone("time_vs_activity");
    h_time_vs_activity  ->GetYaxis()->SetTitle("Equivalent activity [mBq/kg]");
    h_time_vs_rate_bg   = (TH1D*)h_time_vs_rate->Clone("time_vs_rate_BG");
    h_time_vs_rate_bg   ->SetTitle("Background component");
   // if( _isMC ) {
   // }

    // =====================================================
    // Diagnotic and utility histograms
    tdir_util->cd();
    h_2D_time_vs_dT   = new TH2D("2D_time_vs_dT",";Time [hr];#DeltaT [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_2D_time_vs_dT_bg= (TH2D*)h_2D_time_vs_dT->Clone("2D_time_vs_dT_bg");
    h_time_vs_N       = new TH1D("time_vs_N",";Time [hr];Number of entries into dT plot",timeBins,0,timeMax);
    
    // =====================================================
    // MC-truth based histograms
    tdir_truth->cd();
    h_beta_nwires = new TH1D("beta_nwires","True beta cluster;Number of collection plane wires;Entries",10,0,10);
    h_alpha_nwires = new TH1D("alpha_nwires","True alpha cluster;Number of collection plane wires;Entries",10,0,10);
    h_beta_trueEnergy             = new TH1D("beta_trueEnergy",             "All true decays;Electron true energy [MeV];Entries per bin",70,0,betaEmax);
    h_beta_trueEnergySum          = new TH1D("beta_trueEnergySum",          "All true decays;Decay vertex true energy [MeV];Entries per bin",70,0,betaEmax);
    h_beta_trueEnergySum_reco     = new TH1D("beta_trueEnergySum_reco",     "Reco'd on coll plane;Decay vertex true energy [MeV];Entries per bin",70,0,betaEmax);
    h_beta_trueEnergySum_recoCuts = new TH1D("beta_trueEnergySum_recoCuts", "blip cuts;Decay vertex true energy [MeV];Entries per bin",70,0,betaEmax);
    h_true_alpha_depne      = (TH1D*)h_alpha_charge->Clone("true_alpha_depne");
    h_true_alpha_depne      ->SetTitle("True ionization electrons from alpha");
    h_true_alpha_charge     = (TH1D*)h_alpha_charge->Clone("true_alpha_charge");
    h_true_alpha_charge     ->SetTitle("True alpha charge at anode");
    h_matched_alpha_charge  = (TH1D*)h_alpha_charge->Clone("matched_alpha_charge");
    h_matched_alpha_charge  ->SetTitle("Cluster charge matched to alpha");
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
    if( !_isMC ) qscale = 1.0;
    if( qscale < 1. ) printf("WARNING: scalinig charge by %f\n",qscale);

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
    //fTree->SetBranchAddress("clust_starttime",&clust_starttime);
    //fTree->SetBranchAddress("clust_endtime",&clust_endtime);
    fTree->SetBranchAddress("clust_nhits",&clust_nhits);              
    fTree->SetBranchAddress("clust_charge",&clust_charge);            
    fTree->SetBranchAddress("clust_time",&clust_time);                
    fTree->SetBranchAddress("clust_blipid",&clust_blipid);            
    fTree->SetBranchAddress("clust_deadwiresep",&clust_deadwiresep);
    fTree->SetBranchAddress("clust_amp",&clust_amp);
    fTree->SetBranchAddress("nblips",&nblips);                       
    fTree->SetBranchAddress("blip_nplanes",&blip_nplanes);
    fTree->SetBranchAddress("blip_energy",&blip_energy);
    fTree->SetBranchAddress("blip_y",&blip_y);                        
    fTree->SetBranchAddress("blip_z",&blip_z);                        
    fTree->SetBranchAddress("blip_charge",&blip_charge);              
    fTree->SetBranchAddress("blip_pl0_clustid",&blip_clustid[0]);
    fTree->SetBranchAddress("blip_pl1_clustid",&blip_clustid[1]);
    fTree->SetBranchAddress("blip_pl2_clustid",&blip_clustid[2]);
    fTree->SetBranchAddress("blip_yzcorr",&blip_yzcorr);
    if( _isMC ) {
      fTree->SetBranchAddress("nedeps",&nedeps);
      fTree->SetBranchAddress("clust_edepid",&clust_edepid);            
      fTree->SetBranchAddress("nparticles",&nparticles);
      fTree->SetBranchAddress("part_isPrimary",&part_isPrimary);
      fTree->SetBranchAddress("part_startT",part_startT);
      fTree->SetBranchAddress("part_mother",&part_mother);
      fTree->SetBranchAddress("part_pdg",&part_pdg);
      fTree->SetBranchAddress("part_KE",&part_KE);
      fTree->SetBranchAddress("edep_isPrimary",&edep_isPrimary);
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

    // initialize all histograms
    setRootStyle();
    makeHistograms();

    // (find somewhere better to put these)
    _minTick          = (int)fdT_max*2;
    _liveTimePerEvt   = (_maxTick-_minTick)*fSamplePeriod*1e-6; //sec
    if( fSkipNoisyWires == false ) fNoisyWires.clear();
    for(auto iwire : fNoisyWires )  wireIsNoisy[iwire] = true;



    // ****************************************************
    // BEGIN EVENT LOOP
    // ****************************************************
    int print_counter = 0;
    std::time_t loopStart = time(0);
    for(size_t iEvent=0; iEvent < fTree->GetEntries(); iEvent++){

      print_counter++;
      if( print_counter > 1000 ) {
        print_counter = 1;
        printf("========== EVENT %lu / %lu, %6.2f %%, BiPo count: %i =====================\n",
          iEvent,(size_t)fTree->GetEntries(),100*iEvent/float(fTree->GetEntries()),_numBiPo);
      }
     
      
      // ..... quick-test options ...........
      //int maxEvt    = 1000; if(  iEvent >= maxEvt ) break;
      //int sparsify  = 100; if(  (iEvent % sparsify) != 0 ) continue; 
      //..................................

      // Retrieve event info
      // Record the event time relative to start of dataset period
      fTree->GetEntry(iEvent);
      double eventHr = ( timestamp - inFile.t0 ) / 3600.;
      if( !_isMC && (eventHr < fMinHr || eventHr > fMaxHr) ) continue;
      h_time_vs_N->Fill(eventHr);
      _numEvents++;

      
      // ====================================================
      // Map of clust IDs per wire on collection plane, and 
      // some simple cluster cuts and charge scaling
      // ====================================================
      _clustAvailable.assign(nclusts, true);
      _map_wire_clusters.clear();
      for(int i=0; i < nclusts; i++){
        clust_charge[i] *= qscale;
        if( clust_deadwiresep[i] < 0 ) clust_deadwiresep[i]=999; // TEMP (fixed in blipana as of 1/5/23)
        if( clust_deadwiresep[i] < 1 ) _clustAvailable[i]=false;
        if( clust_plane[i] != 2 ) continue;
        h_wt_clusts->Fill(clust_startwire[i],clust_time[i]);
        for(int j=clust_startwire[i]; j<=clust_endwire[i]; j++){
          _map_wire_clusters[j].push_back(i);
          h_nclusts_perwire->Fill(j);
        }
      }
  

      // ======================================
      // Check truth info
      // ======================================
      clust_isAlpha.assign(nclusts, false); 
      clust_isBeta.assign(nclusts, false); 
      clust_isGamma.assign(nclusts, false); 
      if( _isMC ) {
        
        int alphaPDG = 1000020040;
        
        // -------------------------------------------------
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
          // skip any secondaries for now
          if( !part_isPrimary[i] || part_mother[i] > 0 ) continue;
          // Bi214 betas
          if( part_pdg[i]==11 ) {
            nPrimaryElectrons++;
            sumBetaEnergy += KE;
            //h_beta_trueEnergy->Fill(KE);
          }  
          // when we reach a Po214 alpha, reset the counters so
          // on the next loop we can find the next beta
          if( part_pdg[i] == alphaPDG ) {
            h_beta_trueEnergySum->Fill(sumBetaEnergy);
            nPrimaryElectrons = 0;
            sumBetaEnergy = 0;
          }
        }//endloop over MCParticles
        

        //--------------------------------------------------
        // loop over the energy deposits ("true blips")
        for(int i=0; i<nedeps; i++){
          if( !edep_isPrimary[i] ) continue;
          int q_dep   = edep_electrons[i];
          int q_drift = edep_charge[i];
          if( edep_pdg[i] == alphaPDG ) {
            _numBiPo_true++;
            if( q_drift > 0 ) h_true_alpha_charge ->Fill(q_drift/0.826);
            if( q_dep > 0   ) h_true_alpha_depne  ->Fill(q_dep);
            float driftTime   = edep_tdrift[i];
            float readoutTime = driftTime + part_startT[i];
            int   readoutTick = readoutTime/fSamplePeriod;
            if( readoutTick > _minTick && readoutTick < _maxTick ) 
              _numBiPo_true_perfectReco++;
          }
        }
        
        //----------------------------------------------
        // look for collection plane clusts matched to an alpha, beta, or gamma
        for(int i=0; i < nclusts; i++){
          //if( clust_plane[i] != 2 ) continue;
          int eid = clust_edepid[i];
          if( eid < 0 ) {
            if( fIgnoreNonMC ) _clustAvailable[i] = false;
            continue;
          }
          int g4index = edep_g4id[eid];
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
          
          if( clust_plane[i] != 2 ) continue;
          if( clust_isBeta[i] ) {
            h_beta_nwires->Fill(clust_nwires[i]);
            h_beta_trueEnergySum_reco->Fill(edep_energy[eid]);
          }
          if( clust_isAlpha[i] ) h_alpha_nwires->Fill(clust_nwires[i]);

        }//clust loop
      }//isMC

      
      // ==============================================================
      // Create list of blip IDs sorted by charge
      // ==============================================================
      /*
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
        float beta_charge = blip_charge[iBlip]; //*corrFact;
        float beta_energy = blip_energy[iBlip]; //( beta_charge / fRecomb ) * (23.6e-6);
        int   beta_wires  = clust_nwires[ic];

        // blip plane-matching requirements
        if( blip_nplanes[iBlip] < fBetaMinPlanes ) continue;
        
        // evaluate if in fiducial volume
        if( fFidVolCut ) {
          if( blip_z[iBlip] < fZlim[0] || blip_z[iBlip] > fZlim[1] ) continue;
          if( blip_y[iBlip] < fYlim[0] || blip_y[iBlip] > fYlim[1] ) continue;
        }
        
        // apply charge/size cuts on beta
        if(   beta_wires  > fBetaWires_max  ||
              //beta_energy < fBetaEnergy_min || 
              beta_energy > fBetaEnergy_max ) 
          continue;


  			// skip if any wire in cluster is from a noisy wire
        // skip if this cluster is at the very edge of the wireplane
        int w_start = std::max(0,clust_startwire[ic]-1);
        int w_end   = std::min(clust_endwire[ic]+1,(int)nWiresColl-1);
        if( w_start < 10 || w_end > int(nWiresColl-10) ) continue;
        bool flag = false;
        for(int iwire = w_start; iwire <= w_end; iwire++){
          if( wireIsNoisy[iwire] ) {
            flag = true;
            break;
          }
        }
        if( flag ) continue;
        
        
        
        // require semi-isolation (no other clusters within +/-20us)
        for(int iwire = w_start; iwire <= w_end; iwire++){
          for(auto& jc : _map_wire_clusters[iwire] ) {
            if( ic == jc ) continue;
            float dt = fSamplePeriod*fabs(clust_time[jc]-clust_time[ic]);
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
        


        // ------------------------------------------------------------
        // Begin search for alpha candidates
        // ------------------------------------------------------------
        int nclusts_inwindow        = 0;
        int nwires = 0;
        std::vector<BiPoCandidate> v_cands = FindCandidates(ic, 0, fWireRange, false, nclusts_inwindow,nwires);
        h_nclusts_inwindow  ->Fill(float(nclusts_inwindow)/nwires);
        h_ncands_inwindow   ->Fill((int)v_cands.size());
       
        bool passes = false;
        for(auto& thisCand : v_cands) {
          //_clustAvailable[ic] = false;
          //_clustAvailable[thisCand.id2] = false;
          
          if( thisCand.dT < fdT_min ) continue;
          h_beta_charge   ->Fill(thisCand.q1);
          h_beta_energy   ->Fill(thisCand.e1);
          h_beta_amp[0]   ->Fill(clust_amp[ic]);

          if( thisCand.e1 < fBetaEnergy_min ) continue;
          h_alpha_charge  ->Fill(thisCand.q2);
          h_alpha_energy  ->Fill(thisCand.e2);
          h_alpha_amp[0]  ->Fill(clust_amp[thisCand.id2]);
          h_alpha_energyVsdT->Fill(thisCand.dT,thisCand.e2);
          if( false && thisCand.dT < 40 ){
          std::cout<<"Evt "<<iEvent<<", hour "<<eventHr<<", beta clust/wire: "<<ic<<"   "<<clust_startwire[ic]<<" - "<<clust_endwire[ic]
            <<";   alpha clust/wire: "<<thisCand.id2<<"   "<<clust_startwire[thisCand.id2]<<"; dT "<<thisCand.dT<<"\n";
            std::cout<<"Beta amp: "<<clust_amp[ic]<<",   alpha amp: "<<clust_amp[thisCand.id2]<<"\n";
            int beta_start = clust_starttime[ic];
            int beta_end  = clust_endtime[ic];
            int alpha_start = clust_starttime[thisCand.id2];
            int alpha_end  = clust_endtime[thisCand.id2];
            std::cout<<"Beta start/end: "<<beta_start<<" - "<<beta_end<<"       alpha start/end: "<<alpha_start<<" - "<<alpha_end<<"\n";
          }
          passes=true;
          h_zy_bipos    ->Fill( blip_z[iBlip], blip_y[iBlip]);
          h_wt_bipos    ->Fill( clust_startwire[ic], clust_time[ic]);
          h_cand_dT     ->Fill(thisCand.dT);
          h_2D_time_vs_dT->Fill(eventHr,thisCand.dT);
          //h_alpha_nwires        ->Fill(clust_nwires[thisCand.id2]);
        }
        if( passes ) { 
          _numBiPo++;
          if( clust_isBeta[ic] ) _numBiPo_mcmatch++;
        }
        
        // --------------------------------------------
        // Evaluate dT-flip candidates
        // ---------------------------------------------
        std::vector<BiPoCandidate> v_cands_bg = FindCandidates(ic,0, fWireRange, true, nclusts_inwindow, nwires);
        for(auto& thisCand : v_cands_bg) {
          //_clustAvailable[ic] = false;
          //_clustAvailable[thisCand.id2] = false;
          if( thisCand.dT < fdT_min ) continue;
          h_beta_charge_bg  ->Fill(thisCand.q1);
          h_beta_energy_bg  ->Fill(thisCand.e1);
          h_beta_amp[1]     ->Fill(clust_amp[ic]);
          if( thisCand.e1 < fBetaEnergy_min ) continue;
          h_alpha_charge_bg ->Fill(thisCand.q2);
          h_alpha_energy_bg ->Fill(thisCand.e2);
          h_alpha_amp[1]    ->Fill(clust_amp[thisCand.id2]);
          h_zy_bipos_bg     ->Fill( blip_z[iBlip], blip_y[iBlip]);
          h_cand_dT_bg    ->Fill(thisCand.dT);
          h_2D_time_vs_dT_bg->Fill(eventHr,thisCand.dT);
        }
        
        // --------------------------------------------
        // Evaluate 'control region' for detector effects
        // (look in region shifted + 5 wires away)
        // ---------------------------------------------
        std::vector<BiPoCandidate> v_control      = FindCandidates(ic, 5, 2, false, nclusts_inwindow, nwires);
        std::vector<BiPoCandidate> v_control_bg = FindCandidates(ic, 5, 2, true,  nclusts_inwindow, nwires);
        for( auto& vc : v_control     ) {
          if(vc.dT < fdT_min || vc.e1 < fBetaEnergy_min ) continue;
          h_control_dT->Fill( vc.dT );
        }
        for( auto& vc : v_control_bg  ) {
          if(vc.dT < fdT_min || vc.e1 < fBetaEnergy_min ) continue;
          h_control_dT_bg  ->Fill( vc.dT );
        }



      }//end loop over 3D blips

    }//endloop over events
    double loopDuration = ( time(NULL) - loopStart );


    // ***************************************************
    // Scale everything so it's 'per second in full AV'
    // ***************************************************
    _totalLiveTime = float(_numEvents) * _liveTimePerEvt;
    float scaleFact = (1./_totalLiveTime)*(1./_fiducialFrac);
    h_cand_dT           ->Scale( scaleFact); //, "width");
    h_cand_dT_bg      ->Scale( scaleFact); //, "width");
    h_beta_charge       ->Scale( scaleFact );
    h_beta_charge_bg    ->Scale( scaleFact );
    h_beta_energy       ->Scale( scaleFact );
    h_beta_energy_bg    ->Scale( scaleFact );
    h_beta_amp[0]       ->Scale( scaleFact );
    h_beta_amp[1]       ->Scale( scaleFact );
    h_alpha_charge      ->Scale( scaleFact );
    h_alpha_charge_bg   ->Scale( scaleFact );
    h_alpha_energy      ->Scale( scaleFact );
    h_alpha_energy_bg   ->Scale( scaleFact );
    h_alpha_amp[0]      ->Scale( scaleFact );
    h_alpha_amp[1]      ->Scale( scaleFact );

    // Keep sum squares
    h_control_dT        ->Sumw2();
    h_control_dT_bg     ->Sumw2();
    h_control_dT_ratio  ->Divide(h_control_dT,h_control_dT_bg);
    h_control_dT_ratio  ->Fit(f_backward_corr,"R");

    if( fLinearizeCorr )  h_cand_dT_bg->Multiply( f_backward_corr );
    else                  h_cand_dT_bg->Multiply( h_control_dT_ratio );
   
    // ***************************************************
    // Histogram subtraction time!
    // ***************************************************
    h_cand_dT_sub         ->Add(h_cand_dT, h_cand_dT_bg, 1, -1);
    h_alpha_charge_sub    ->Add(h_alpha_charge,   1);
    h_beta_charge_sub     ->Add(h_beta_charge,    1);
    h_alpha_energy_sub    ->Add(h_alpha_energy,   1);
    h_beta_energy_sub     ->Add(h_beta_energy,    1);
    h_zy_bipos_sub         ->Add(h_zy_bipos,      1);
    //h_beta_ratio_sub      ->Add(h_beta_ratio,     1);

    //h_cand_dT_sub->Add(h_cand_dT_bg,-1.); 

    h_alpha_charge_sub->Add(h_alpha_charge_bg,-1.);
    h_beta_charge_sub->Add(h_beta_charge_bg,-1.);
    h_alpha_energy_sub->Add(h_alpha_energy_bg,-1.);
    h_beta_energy_sub->Add(h_beta_energy_bg,-1.);
    h_zy_bipos_sub->Add(h_zy_bipos_bg,-1);
    //h_beta_ratio_sub->Add(h_beta_ratio_bg,-1.);
    
    h_beta_amp[2]   ->Add(h_beta_amp[0],  h_beta_amp[1],  1,  -1);
    h_alpha_amp[2]  ->Add(h_alpha_amp[0], h_alpha_amp[1], 1,  -1);
  
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
    printf("Processing time     : %f sec (%f sec/evt)\n", loopDuration, loopDuration/float(_numEvents));
    printf("Excluded %i noisy wires \n",     (int)fNoisyWires.size());	
    printf("*******************************************\n\n");
    fOutFile->Close();
  
    std::cout<<"Backward correction at 500: "<<f_backward_corr->Eval(500)<<"\n";

    if( _isMC ) {
      float rate_Bq       = _numBiPo_true_perfectReco/_totalLiveTime;
      float rate_readout  = _numBiPo_true_perfectReco/float(_numEvents);
      float rate_sim      = _numBiPo_true/(_numEvents*0.0056);
      std::cout<<"Simulated rate: "<<rate_sim<<" per sec --> "<<rate_sim*0.0032<<" per 3.2ms\n";
      std::cout<<"BiPos occuring in dT window: "<<_numBiPo_true_perfectReco<<" --> "<<rate_readout<<" per 3.2ms\n)";
    }
    
  }
  
  
  
  //#################################################################################
  // Make plots here
  //#################################################################################
  void makePlots()
  {
    // Histograms needed:
    //   - h_cand_dT_sub
    //   - h_2D_time_vs_dT
    //   - h_2D_time_vs_dT_bg
    //   - h_time_vs_N
    //   - h_alpha_charge_sub
    //   - h_beta_charge_sub

    std::string name;
    float       range, min, max;
    
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
        h_bg        = Make1DSlice( h_2D_time_vs_dT_bg, i, i, Form("time_vs_dT_bg_%i",i) );
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
      
      std::cout<<"\nRate over time: \n";
      std::cout<<" 2 hours: "<<h_time_vs_rate->Interpolate(2)<<"\n";
      std::cout<<" 4 hours: "<<h_time_vs_rate->Interpolate(4)<<"\n";
      std::cout<<" 6 hours: "<<h_time_vs_rate->Interpolate(6)<<"\n";
      std::cout<<" 8 hours: "<<h_time_vs_rate->Interpolate(8)<<"\n";
      std::cout<<" 10 hours: "<<h_time_vs_rate->Interpolate(10)<<"\n\n";


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
      TCanvas* canvas = new TCanvas(name.c_str(),name.c_str(),500,380);
      gStyle->SetOptStat(0);
      gPad->SetGridy(1);
      h_p0->GetYaxis()->SetRangeUser(min-0.3*range, max+0.5*range);
      h_p0->GetYaxis()->SetTitleOffset(1.1); 
      h_p0->GetYaxis()->SetTitle("Fit component integral");
      h_p0->DrawCopy();
      h_p1->DrawCopy("same");
      canvas->Write();


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
    
    std::cout<<"\n\nFitting "<<h->GetTitle()<<", "<<h->GetEntries()<<"\n";
    FitResult out;
    std::string label = h->GetName();
    
    TCanvas* canvas = new TCanvas(Form("c_fit_%s",label.c_str()),Form("c_fit_%s",label.c_str()),500,420);
    float histMax = GetHistMax(h);
    float histMin = GetHistMin(h);
    float range   = (histMax-histMin);
   
    TGraphErrors* gr = MakeGraph(h);
    
    // Full fit to un-subtracted dT spectra
    //TF1* fit = new TF1("FullFit","[0] + [1]*exp(-x/[2]) + [3]*exp(-x/[4])",fdT_min,fdT_max);
    //fit->SetParameter(0, histMax/5 );
    //fit->SetParLimits(0, 0, histMax );
    //fit->SetParameter(3, histMax/5 );
    //fit->SetParLimits(3, 0, histMax );
    //fit->SetParameter(4, 15 );
    //fit->SetParLimits(4, 5, 40);
    //fit->SetParameter(1, histMax );
    //fit->SetParLimits(1, 0, histMax*2 );
    //fit->FixParameter(2, 164.5 );

    // Define single exp fit + flat BG
    TF1* fit = new TF1("FullFit","[0] + [1]*exp(-x/[2])",fdT_min,fdT_max);
    fit->SetParameter(0, (histMax+histMin)/2 );
    fit->SetParameter(1, (histMax-histMin) );
    fit->FixParameter(2, 164.3 );
    if( constrainNorm ) fit->SetParLimits(1, 0., histMax*10);

    // Draw plot and fit
    canvas->cd();
    gr->Fit(fit,"REQ");
  
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


    gr->Draw("AP");
    gStyle->SetOptFit(1112);
    if( writeCanvas ) {
      tdir_util->cd();
      gr->Write(label.c_str());
      tdir_plots->cd();
      canvas->Write();
      fOutFile->cd();
    }
    

    // Full extrapolated BiPo component; integral of A*exp(-x/B) from 0-infinity is A*B
    double p0       = fit->GetParameter(0);
    double p0_err   = fit->GetParError(0);
    double p1       = fit->GetParameter(1);
    double p1_err   = fit->GetParError(1);
    double tau      = fit->GetParameter(2);
    double tau_err  = fit->GetParError(2); // tau is fixed so this should be =0
    double chi2     = fit->GetChisquare();
    int    ndf      = fit->GetNDF();
    double n_bipo   = (p1*tau)/fdT_binSize;
    double n_bipo_err = (p1_err*tau)/fdT_binSize;

    // Factor out the signal and BG components into separate TF1 functions
    TF1* f_bipo = new TF1("bipo","[0]*exp(-x/[1])");
    TF1* f_flat = new TF1("flat","[0]");
    f_flat->FixParameter(0,p0);   f_flat->SetParError(0,p0_err);
    f_bipo->SetParameter(0,p1);   f_bipo->SetParError(0,p1_err);
    f_bipo->SetParameter(1,tau);  f_bipo->SetParError(1,tau_err);
  
    // Components within dT selection window (to be used for S/N metrics)
    //double N_total    = h     ->Integral(1,h->GetXaxis()->GetNbins());
    double N_total    = fit   ->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_bipo     = f_bipo->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_bg       = f_flat->Integral(fdT_min,fdT_max)/fdT_binSize;
    double N_bipo_err = fabs(N_bipo*(n_bipo_err/n_bipo));
    double N_bg_err   = fabs(N_bg*(p0_err/p0));
    double SNratio    = N_bipo / N_bg;

    // As extra systematic, fit the flat comp to zero
    fit->FixParameter(0,0); 
    gr->Fit(fit,"QN");
    double n2 = (fit->GetParameter(1)*fit->GetParameter(2))/fdT_binSize;
    double n_bipo_syst = fabs(n_bipo-n2);
   
    // efficiency factor + error
    float effMC        = (fEfficiencyMC>0) ? fEfficiencyMC : 1.;
    float effMC_syst    = (fEfficiencyMC>0) ? std::max(0.,fEfficiencyMC_err) : 0.;
  
    // Fractional error terms for specific activity calculation
    // - err1:  stat error from fit
    // - err2a: syst error from fit
    // - err2b: syst error from efficiency
    double err1   = (n_bipo!=0) ? fabs(n_bipo_err/n_bipo) : 0.;
    double err2a  = (n_bipo!=0) ? fabs(n_bipo_syst/n_bipo) : 0.;
    double err2b  = (effMC>0) ? effMC_syst/effMC : 0.;

    double activity       = 1e3 * n_bipo / effMC / 85000.;
    double activity_err   = fabs(activity) * err1;
    double activity_syst  = fabs(activity) * sqrt( pow(err2a,2) + pow(err2b,2) );
    double activity_errTot = sqrt( pow(activity_err,2) + pow(activity_syst,2) );
    
    // normalize "rate" to be per 3.2ms readout
    n_bipo      *= 0.0032;
    n_bipo_err  *= 0.0032;
    n_bipo_syst *= 0.0032;
    
    printf("================ dT fit results =================\n");
    printf("p0          : %f +/- %f\n", p0, p0_err);
    printf("p1          : %f +/- %f\n", p1, p1_err);
    printf("Chi2/ndf    : %f / %i = %f\n",chi2,ndf,chi2/ndf);
    printf("Total rate  : %f per sec\n",N_total);
    printf(" - BiPo     : %f per sec\n",N_bipo);
    printf(" - Flat     : %f per sec\n",N_bg);
    printf("S/N         : %f\n",SNratio);
    printf("Rate        = %f +/- %f (stat) +/- %f (syst) BiPos per 3.2 ms\n",n_bipo,n_bipo_err,n_bipo_syst);
    printf("            = %f +/- %f BiPos per 3.2 ms\n",n_bipo,sqrt(pow(n_bipo_err,2)+pow(n_bipo_syst,2)));
    printf("Activity*   = %f +/- %f (stat) +/- %f (syst) mBq/kg\n",activity,activity_err,activity_syst);
    printf("            = %f +/- %f mBq/kg\n",activity,activity_errTot);
    printf("              *efficiency: %.3f +/- %.3f\n",effMC,effMC_syst);

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
  std::vector<BiPoCandidate> FindCandidates(int ic, int wire_shift, int range, bool flipdT, int& npileup, int& nwires ) {
   
    //std::cout<<"looking for cands; beta clust "<<ic<<", beta range "<<clust_startwire[ic]<<", "<<clust_endwire[ic]<<"... shift "<<wire_shift<<"    range "<<range<<"\n";
    std::vector<BiPoCandidate> v;
    std::vector<int> cand_IDs;
    npileup = 0;
    nwires = 0;
    
    int blipID = clust_blipid[ic];

    std::set<int> wires;
    int w1  = clust_startwire[ic] - wire_shift;
    int w2  = clust_endwire[ic] + wire_shift;
    
    //std::cout<<"w1 w2 "<<w1<<"  "<<w2<<"\n";
     

    for(int i=w1-range; i<=w1+range; i++) wires.insert(i);
    for(int i=w2-range; i<=w2+range; i++) wires.insert(i);

    nwires = 0;

    for(auto iWire : wires ) {
      //std::cout<<"checking wire "<<iWire<<"\n";

      if( wireIsNoisy[iWire] ) continue;
      nwires++;
      

      //std::cout<<"There are "<<_map_wire_clusters[iWire].size()<<" clusts on this wire\n";
      for(auto& jc : _map_wire_clusters[iWire] ) {
        
        if( ic == jc ) continue;
        if( !_clustAvailable[jc] ) continue;
        if( std::count(cand_IDs.begin(),cand_IDs.end(),jc) ) continue;

        float dT        = (clust_time[jc]-clust_time[ic])*fSamplePeriod;
        if (flipdT) dT  *= -1.;
        
        if( dT < 0 || dT > fdT_max ) continue;
        npileup++;
        //if( dT < fdT_min ) continue;
        
        float beta_q = qscale*blip_charge[blipID];
        float beta_E = qscale*blip_energy[blipID];
        float alpha_q = clust_charge[jc]*blip_yzcorr[blipID];
        float alpha_E = charge_to_energy(alpha_q,fRecomb);
        if( clust_nwires[jc]  > fAlphaWires_max   ) continue;
        if( alpha_E           < fAlphaEnergy_min  ) continue;
        if( alpha_E           > fAlphaEnergy_max  ) continue;
        //if( alpha_E           > beta_E )            continue;
        
        BiPoCandidate c = { clust_blipid[ic], ic, jc, dT, beta_q, alpha_q, beta_E, alpha_E};
        v.push_back(c);
        cand_IDs.push_back(jc);
        
      }//<-- end loop over clusters on this wire
    }//<-- endloop over wires
  
    return v;
  
  }

  //###################################################################
  void configureMacro(){
  
  }


  //###################################################################
  //int FindG4Index(int g4id) {
  //  for(size_t i=0; i<nparticles; i++) {
  //    if( part_trackID[i] == g4id ) return i;
  //  }
  //  return -9;
  //}


