  //////////////////////////////////////////////////////////////////
  // 
  //  Analysis ROOT Macro
  //
  ////////////////////////////////////////////////////////////////// 
  
  #include "core/vars.h"
  #include "core/tools.h"
  #include <time.h>
  
  // --- Choose dataset ---
  int   fPhase  = 2;
  
  // --- Input files ---
  std::string fFileName[2]= { "files/BlipAna_20220405_RadonData_Phase1.root",
                              "files/BlipAna_20220405_RadonData_Phase2.root" };
  
  std::string fTreeName = "blipanaTrkMask/anatree";

  // --- Time periods in UNIX time ---
  unsigned int fPhase_T0[2] = {  1627415210, 1627594369 };
  unsigned int fPhase_T1[2] = {  1627592728, 1627761265 };
  //                            Phase 1     Phase 2
  //                            ~49.3 hrs   ~46.4 hrs
  
  // --- General macro parameters ---
  int   fSubtractionMode   = 2;     // 0= none; 1= subtract flip; 2= subtract shift 
  bool  fPickyBlipMode    = false;   // Require blips match on all 3 planes
  int   fMinTick          = 0;      // min waveform tick to start search (0-6400 ticks)
  int   fWireRange        = 1;      // +/- range to look for alpha candidate
  float fMaxClusterSpan   = 30;     // Veto clusters longer than this [ticks]
  float fMinGOF           = -999;   // pulse train hits have gof < 0
  float fBetaCharge_min   = 2e3;    // Min charge of beta candidate blip [e-]
  float fBetaCharge_max   = 100e3;  // Max charge of beta candidate blip [e-]
  float fAlphaCharge_min  = 0e3;    // Min charge of alpha candidate cluster [e-]
  float fAlphaCharge_max  = 5e3;   // Max charge of alpha candidate cluster [e-]
  int   fAlphaHits_max    = 999;    // Max hits in alpha candidate cluster
  bool  fAlphaReq3D       = false;  // Require alpha pulse be 3D-matched
  float fdT_binSize       = 20.;    // Bin width for all dT spectra plots [us]
  float fdT_min           = 20.;    // Min dT for looking for candidate [us]
  float fdT_max           = 800.;   // Max dT for looking for candidate [us]
  int   fMaxClustMult     = 5; //5      // Max number of clusters in time window
  int   fMaxCandidates    = 1; //3      // Max number of alpha-like candidates
  float fSphereRadius     = 25;     
  int   fMinSphereMult3D  = -1; //0;
  int   fMaxSphereMult3D  = -1; //10;
  float fZmin             = 50; //50;     // Z range (0 to 1037 cm)
  float fZmax             = 985; //1035; //985;    //
  float fYmin             = -70; //-120; //-70;    // Y range (-120 to 120 cm)
  float fYmax             = 70; //120;     //
  
  // --- Detector properties ---
  int   nWiresColl        = 3455;
  float fSamplePeriod     = 0.5; // microseconds
  
  // --- Special switches ---
  bool  fDoWireDiagnostics= false;
  int   fRandomWireShift  = 0; 
  
  // --- Noisy wires to skip (collection plane) ---
  std::vector<int> fNoisyWires{
  //0, 1, 2, 83, 384, 1247, 1535, 1540, 1919, 1920, 2303, 2334, 2335, 2400, 2415, 2465, 
  //2687, 3071, 3072, 3408, 3412, 3422, 3423, 3424, 3425, 3426, 3428, 3429, 3432, 3433, 
  //3435, 3437, 3438, 3439, 3442, 3444, 3445, 3446, 3448, 3449, 3452, 3453, 3455
  //};
  0, 1, 2, 3, 6, 7, 9, 10, 11, 16, 72, 80, 83, 85, 384, 1151, 1152, 1247, 1535, 1540, 
  1716, 1919, 1920, 2111, 2229, 2230, 2297, 2298, 2300, 2301, 2303, 2304, 2309, 2311, 
  2327, 2329, 2330, 2333, 2334, 2335, 2400, 2415, 2464, 2465, 2681, 2686, 2687, 2705, 
  2732, 2733, 2753, 2764, 2783, 2789, 2879, 3071, 3072, 3215, 3251, 3263, 3269, 3274, 
  3280, 3286, 3289, 3294, 3295, 3299, 3313, 3318, 3327, 3335, 3336, 3338, 3339, 3342, 
  3345, 3347, 3348, 3351, 3353, 3359, 3360, 3361, 3378, 3385, 3391, 3408, 3409, 3410, 
  3411, 3412, 3413, 3415, 3416, 3418, 3419, 3420, 3421, 3422, 3423, 3424, 3425, 3426, 
  3427, 3428, 3429, 3430, 3431, 3432, 3433, 3434, 3435, 3436, 3437, 3438, 3439, 3441, 
  3442, 3443, 3444, 3445, 3446, 3447, 3448, 3449, 3451, 3452, 3453, 3454, 3455
  };
  
  //#######################################################################
  // Derived parameters
  //######################################################################
  
  // File names and T0
  std::string   _fileName = fFileName[fPhase-1];
  unsigned int  _t0       = fPhase_T0[fPhase-1];
  
  // Fiducial vol correction factor
  float dz = fZmax-fZmin;
  float dy = fYmax-fYmin;
  float fFidCorFactor = std::max(1., (1037.*240.) / (dz*dy) );
  
  // Counters
  int   _numEvents      = 0;
  float _numBiPo        = 0;
  float _numBiPo_6_312  = 0;
  float _totalLiveTime  = 0;
  std::vector<bool> _blipAvailable(kMaxBlips, true);
  std::vector<bool> _clustAvailable(kMaxHits, true);

  // Map of cluster IDs per wire on collection plane
  std::map<int,std::vector<int>> _map_wire_clusters;
  
  
  // Live time
  int   _minTick         = fMinTick;
  int   _maxTick         = 6400 - (int)fdT_max*2;
  float _liveTimePerEvt  = (_maxTick-_minTick)*fSamplePeriod*1e-6; //sec
  
  // Blip fidelity requirements
  int _betaMinPlanes = 2;
  int _betaMaxDiff = 9999;
  
  
  //##########################################################################
  // Structs, functions, histograms
  //##########################################################################
  
  // Struct to hold fit results
  struct FitResult { float rate_signal = -9, rate_bg = -9, ratio = -9; };
  
  // useful structure to save candidate info in
  struct BiPoCandidate { int id1, id2; float dT, q1, q2; };
  
  // Functions 
  void                        makePlots();
  void                        makeHistograms();
  std::vector<BiPoCandidate>  FindCandidates(int, int, int, bool, int&);
  FitResult                   fitdT(TH1D*,bool);
  
  // ROOT objects
  TTree*      fTree;
  TFile*      fOutFile;
  TDirectory* fOutFile_plots;
  
  // Histograms
  TH1D* h_blip_charge;
  TH1D* h_blip_nhits;
  TH1D* h_alpha_nhits;
  TH1D* h_clust_timespan;
  TH1D* h_nclusts_wire_ave;
  TH1D* h_nclusts_inwindow;
  TH1D* h_ncands_inwindow;
  TH1D* h_cand_dT;
  TH1D* h_cand_dT_flip;
  TH1D* h_cand_dT_shift;
  TH1D* h_cand_dT_sub;
  TH1D* h_clust_dT;
  TH1D* h_hit_gof_2D;
  TH1D* h_hit_gof_3D;
  TH2D* h_zy_blips;
  TH2D* h_zy_blips_filt;
  TH2D* h_zy_bipos;
  TH2D* h_wt_clusts;
  TH2D* h_wt_blips;
  TH2D* h_wt_blips_filt;
  TH2D* h_wt_bipos;
  TH1D* h_blip_sphere_mult;
  
  
  TH1D* h_time_vs_N;        // basic event count
  
  TH2D* h_time_vs_dT;       // dT for candidate region
  TH2D* h_time_vs_dT_flip;  // dT for background region
  TH2D* h_time_vs_dT_shift;
  TH2D* h_time_vs_dT_sub;

  TH1D* h_time_vs_rate_bipo; 
  TH1D* h_time_vs_rate_bg;
  TH1D* h_time_vs_rate_sum;
  TH1D* h_time_vs_ratio;
  
  TH1D* h_timeFine_vs_evts;
  TH1D* h_timeFine_vs_rate;
  
  //TH2D* h_time_vs_ph;
  TH2D* h_time_vs_ntrks;
  //TH2D* h_time_vs_trklen;
  
  TH1D* h_alpha_charge;
  TH1D* h_alpha_charge_flip;
  TH1D* h_alpha_charge_shift;
  TH1D* h_alpha_charge_sub;
  TH1D* h_beta_charge;
  TH1D* h_beta_charge_flip;
  TH1D* h_beta_charge_shift;
  TH1D* h_beta_charge_sub;
  
  //##########################################################################
  // Initialize histograms
  //##########################################################################
  void makeHistograms()
  {
    if( fDoWireDiagnostics )
      h_nclusts_wire_ave= new TH1D("nclusts_perwire","Average clusters per wire;Average number of clusters per wire",1500,0,3); 
    h_clust_timespan    = new TH1D("clust_timespan","Cluster timespan;Ticks",200,0,200);
    h_blip_charge       = new TH1D("blip_charge","3D Blips;Collection Plane Charge [e];Events", 500,0,50e3);
    h_blip_nhits        = new TH1D("blip_nhits","Candidate blips;Collection Plane Hits",20,0,20);
    h_alpha_nhits       = new TH1D("alpha_nhits","Candidate alphas;Collection Plane Hits",20,0,20);
    h_nclusts_inwindow  = new TH1D("nclusts_inwindow","Number of clusters in time window following Bi-candidate",20,0,20);
    h_ncands_inwindow   = new TH1D("ncands_inwindow","Number of Po candidates in time window following Bi-candidate",10,0,10);
    h_blip_sphere_mult  = new TH1D("blip_sphere_mult",Form("Sphere radius %f cm;3D blip multiplicity",fSphereRadius),50,0,50);
  
    float Zmin = -100;  float Zmax = 1100;  int Zbins = 300;
    float Ymin = -150;  float Ymax = 150;   int Ybins = 75;
    float Tmin = 0;     float Tmax = 6400;  int Tbins = 3200/2; 
    float Wmin = -100;  float Wmax = 3500;  int Wbins = 1800; 
    h_zy_blips      = new TH2D("zy_blips","3D blips;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_blips_filt = new TH2D("zy_blips_filt","3D blips (quality cuts);Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_bipos      = new TH2D("zy_bipos","BiPo candidates;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
    h_zy_blips      ->SetOption("colz");
    h_zy_blips_filt ->SetOption("colz");
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
    //h_cand_dT       = new TH1D("cand_dT","Selected BiPo Candidates;Time Difference [#mus];Candidates per second", dTbins,0,fdT_max);
    h_cand_dT       = new TH1D("cand_dT","Selected BiPo Candidates;Time difference [#mus];Candidates per readout", dTbins,0,fdT_max);
    h_cand_dT_flip  = (TH1D*)h_cand_dT->Clone("cand_dT_flip");  h_cand_dT_flip  ->SetTitle("-dT region candidates");
    h_cand_dT_shift = (TH1D*)h_cand_dT->Clone("cand_dT_shift"); h_cand_dT_shift ->SetTitle("Shifted-wire region candidates");
    h_cand_dT_sub   = (TH1D*)h_cand_dT->Clone("cand_dT_sub");   h_cand_dT_sub   ->SetTitle("Background-subtracted spectrum");
    h_clust_dT      = (TH1D*)h_cand_dT->Clone("clust_dT");      h_clust_dT      ->SetTitle("Same-wire cluster separations");
    
    h_alpha_charge      = new TH1D("alpha_charge","Candidate alphas;Collection Plane Charge [e];Events", 50,0,5e3);
    h_alpha_charge_flip = (TH1D*)h_alpha_charge->Clone("alpha_charge_flip");
    h_alpha_charge_shift = (TH1D*)h_alpha_charge->Clone("alpha_charge_shift");
    h_alpha_charge_sub = (TH1D*)h_alpha_charge->Clone("alpha_charge_sub");
    h_alpha_charge_sub  ->SetTitle("Candidate alphas after background subtraction");
    
    h_beta_charge      = new TH1D("beta_charge","Candidate betas;Collection Plane Charge [e];Events", 500,0,50e3);
    h_beta_charge_flip = (TH1D*)h_beta_charge->Clone("beta_charge_flip");
    h_beta_charge_shift = (TH1D*)h_beta_charge->Clone("beta_charge_shift");
    h_beta_charge_sub = (TH1D*)h_beta_charge->Clone("beta_charge_sub");
    h_beta_charge_sub  ->SetTitle("Candidate betas after background subtraction");
    
    int timeBins    = 9; float timeMax = 45;
    h_time_vs_dT      = new TH2D("time_vs_dT",";Event time [hr];Time difference [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_time_vs_dT      ->SetOption("colz");
    h_time_vs_dT_flip = (TH2D*)h_time_vs_dT->Clone("time_vs_dT_flip");
    h_time_vs_dT_shift= (TH2D*)h_time_vs_dT->Clone("time_vs_dT_flip");
    h_time_vs_dT_sub  = (TH2D*)h_time_vs_dT->Clone("time_vs_dT_flip");
    /*
    h_time_vs_dT_flip = new TH2D("time_vs_dT_flip",";Event time [hr];Time difference [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_time_vs_dT_flip ->SetOption("colz");
    h_time_vs_dT_shift= new TH2D("time_vs_dT_shift",";Event time [hr];Time difference [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_time_vs_dT_shift->SetOption("colz");
    h_time_vs_dT_sub= new TH2D("time_vs_dT_shift",";Event time [hr];Time difference [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
    h_time_vs_dT_shift->SetOption("colz");
    */
    h_time_vs_N       = new TH1D("time_vs_N",";Event time [hr];Number of entries into dT plot",timeBins,0,timeMax);
    //h_time_vs_rate_bipo = new TH1D("time_vs_rate_bipo","BiPo component of dT fit;Event time [hr];Rate [sec^{-1}]",timeBins,0,timeMax);
    h_time_vs_rate_bipo = new TH1D("time_vs_rate_bipo","BiPo component of dT fit;Event time [hr];Rate per readout",timeBins,0,timeMax);
    h_time_vs_rate_bg   = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_rate_BG");
    h_time_vs_rate_bg   ->SetTitle("Background component");
    h_time_vs_ratio     = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_ratio");
    h_time_vs_ratio     ->SetTitle("Signal-to-background ratio");
    h_time_vs_ratio     ->GetYaxis()->SetTitle("Signal-to-background ratio");
    h_time_vs_rate_sum   = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_rate_sum");
    h_time_vs_rate_sum   ->SetTitle("Summed components");
   
    timeBins    = 45; timeMax = 45;
    h_timeFine_vs_evts     = new TH1D("timeFine_vs_evts",";Event time [hr];Number of events",timeBins,0,timeMax);
    h_timeFine_vs_rate  = new TH1D("timeFine_vs_rate",";Event time [hr];Number of candidates per readout",timeBins,0,timeMax);
  
    //h_time_vs_ph    = new TH2D("time_vs_ph",";Event time [hr];Average track hit amplitude per evt [ADC]",timeBins,0,timeMax, 100,0,20);
    //h_time_vs_ph    ->SetOption("colz");
    h_time_vs_ntrks = new TH2D("time_vs_ntrks",";Event time [hr];Number of tracks per evt",timeBins,0,timeMax, 80,0,80);
    h_time_vs_ntrks    ->SetOption("colz");
    //h_time_vs_trklen= new TH2D("time_vs_trklen",";Event time [hr];Track lengths [cm]",timeBins,0,timeMax, 300,0,300);
    //h_time_vs_trklen    ->SetOption("colz");
  
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
    TFile* file = new TFile(_fileName.c_str(),"READ");
    fTree = (TTree*)file->Get(fTreeName.c_str());
  
    // set branches
    fTree->SetBranchAddress("event",&event);                         
    fTree->SetBranchAddress("run",&run);                             
    fTree->SetBranchAddress("timestamp",&timestamp);
    fTree->SetBranchAddress("ntrks",&ntrks);
    fTree->SetBranchAddress("trk_length",&trk_length);
    //fTree->SetBranchAddress("ave_trkhit_ph",&ave_trkhit_ph);
    fTree->SetBranchAddress("nclusts",&nclusts);                     
    fTree->SetBranchAddress("clust_plane",&clust_plane);              
    fTree->SetBranchAddress("clust_wire",&clust_wire);                
    fTree->SetBranchAddress("clust_nhits",&clust_nhits);              
    fTree->SetBranchAddress("clust_charge",&clust_charge);            
    fTree->SetBranchAddress("clust_time",&clust_time);                
    fTree->SetBranchAddress("clust_lhit_peakT",&clust_lhit_peakT);          
    fTree->SetBranchAddress("clust_lhit_gof",&clust_lhit_gof);          
    fTree->SetBranchAddress("clust_timespan",&clust_timespan);          
    fTree->SetBranchAddress("clust_blipid",&clust_blipid);            
    fTree->SetBranchAddress("nblips",&nblips);                       
    fTree->SetBranchAddress("blip_nplanes",&blip_nplanes);
    fTree->SetBranchAddress("blip_maxdiff",&blip_maxdiff);
    fTree->SetBranchAddress("blip_x",&blip_x);                        
    fTree->SetBranchAddress("blip_y",&blip_y);                        
    fTree->SetBranchAddress("blip_z",&blip_z);                        
    fTree->SetBranchAddress("blip_charge",&blip_charge);              
    fTree->SetBranchAddress("blip_clustid_pl2",&blip_clustid_pl2);
    fTree->SetBranchAddress("blip_trkdist",&blip_trkdist);
  
    // make output file to store plots
    fOutFile = new TFile("output/plots_bipo.root", "recreate");
    fOutFile_plots  = fOutFile->mkdir("plots");
    
    // initialize all histograms
    makeHistograms();
 
    // if we're doing dT-flip subtraction, make sure there's enough room
    if( fSubtractionMode == 1 ) _minTick         = std::max(fMinTick, int(fdT_max)*2);
    
    // for picky blip mode, set proper restrictions
    if( fPickyBlipMode ) {
      _betaMinPlanes    = 3;      // Min number of matched planes (must be 2 or 3)
      _betaMaxDiff      = 2;      // Difference in wire intersection points [cm]
    } 

    // Keep list of cluster counts per wire/channel
    std::map<int,int> map_wn;
    //std::map<int,int> map_chn;
    
    // Record start-time
    std::time_t loopStart = time(0);
    
    // Print-interval counter
    int print_counter = 0;
 


    // ****************************************************
    // Loop over the events
    // ****************************************************
    size_t totalEntries = fTree->GetEntries();
    for(size_t iEvent=0; iEvent < totalEntries; iEvent++){
      
      print_counter++;
      if( print_counter > 1000 ) {
        printf("========== EVENT %lu / %lu, BiPo count: %.0f =====================\n",iEvent,totalEntries,_numBiPo);
        print_counter = 1;
      }
      
      // ..... quick-test options ...........
      //int maxEvt    = 1000; if(  iEvent >= maxEvt ) break;
      //int sparsify  = 100; if(  (iEvent % sparsify) != 0 ) continue; 
      //..................................

      // Retrieve event info
      fTree->GetEntry(iEvent);
      _numEvents++;

      // Record the event time relative to start of dataset period
      double eventHour = ( timestamp - _t0 ) / 3600.;
      h_time_vs_N->Fill(eventHour);
      h_timeFine_vs_evts->Fill(eventHour);
  
      // Track-based diagonstics
      h_time_vs_ntrks->Fill(eventHour,ntrks);
     
      // ====================================================
      // Map of clust IDs per wire on collection plane
      // ====================================================
      _map_wire_clusters.clear();
      for(int i=0; i < nclusts; i++){
        if( clust_plane[i] != 2 ) continue;
        _map_wire_clusters[clust_wire[i]].push_back(i);
        map_wn[clust_wire[i]] += 1;
        //map_chn[clust_chan[i]] += 1;
        h_clust_timespan->Fill(clust_timespan[i]);
        h_wt_clusts->Fill(clust_wire[i],clust_lhit_peakT[i]);
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
      std::fill(_blipAvailable.begin(),_blipAvailable.end(), true);
      std::fill(_clustAvailable.begin(),_clustAvailable.end(), true);
  
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
  
        // plot ZY position
        h_zy_blips->Fill( blip_z[iBlip], blip_y[iBlip] );
        h_wt_blips->Fill( clust_wire[ic], clust_lhit_peakT[ic] );
  
        // skip if this cluster was already included in a BiPo candidate 
        if( !_clustAvailable[ic] ) continue;
        
        // skip pulse train hits
        if( clust_lhit_gof[ic] < fMinGOF  || clust_timespan[ic] > fMaxClusterSpan ) continue;
  
  			// skip if this is a noisy wire
  		  if( std::find(fNoisyWires.begin(), fNoisyWires.end(), clust_wire[ic] ) != fNoisyWires.end() ) continue;
  
        // 3D match cuts for beta
        if( blip_nplanes[iBlip] < _betaMinPlanes || blip_maxdiff[iBlip] > _betaMaxDiff ) continue;
        
        // fill some blip/cluster location histograms
        h_zy_blips_filt->Fill( blip_z[iBlip], blip_y[iBlip] );
        h_wt_blips_filt->Fill( clust_wire[ic], clust_lhit_peakT[ic] );
        h_blip_charge->Fill(clust_charge[ic]);
        h_blip_nhits->Fill(clust_nhits[ic]);
  
        // apply charge/size cuts on beta
        if(   clust_charge[ic] < fBetaCharge_min 
          ||  clust_charge[ic] > fBetaCharge_max ) continue;
        
        // evaluate if in fiducial volume
        if( blip_z[iBlip] < fZmin || blip_z[iBlip] > fZmax ) continue;
        if( blip_y[iBlip] < fYmin || blip_y[iBlip] > fYmax ) continue;
      
        // check for nearby 3D blips
        int nMultSphere3D = 0;
        for(int ii=0; ii<nblips; ii++){
          if( ii == iBlip ) continue;
          if( blip_nplanes[ii] < _betaMinPlanes || blip_maxdiff[ii] > _betaMaxDiff ) continue;
          TVector3 p1(blip_x[iBlip],blip_y[iBlip],blip_z[iBlip]);
          TVector3 p2(blip_x[ii],blip_y[ii],blip_z[ii]);
          float ds = (p2-p1).Mag();
          if(ds<fSphereRadius) nMultSphere3D++; 
        }
        h_blip_sphere_mult->Fill(nMultSphere3D);
        if( fMinSphereMult3D > 0 && nMultSphere3D < fMinSphereMult3D ) continue;
        if( fMaxSphereMult3D > 0 && nMultSphere3D > fMaxSphereMult3D ) continue;
      
        // skip if we are near end of wire (account for 400us/800tick trigger offset)
        if( (clust_lhit_peakT[ic]) < _minTick ) continue;
        if( (clust_lhit_peakT[ic]) > _maxTick ) continue;
      

        // ------------------------------------------------------------
        // Begin search for alpha candidates
        // ------------------------------------------------------------
        int nclusts_inwindow        = 0;
        int nclusts_inwindow_flip   = 0;
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
        // Search for background candidates (dT flip)
        std::vector<BiPoCandidate> v_cands_flip = FindCandidates(ic, w0, w1, true, nclusts_inwindow_flip);
  
        // -------------------------------------------
        // Search for background candidates (wire shift)
        int shift = 2*fWireRange+1;
        std::vector<BiPoCandidate> v_cands_shift  = FindCandidates(ic, w0+shift, w1+shift, false, nclusts_inwindow_shift);
        std::vector<BiPoCandidate> v_cands_shift2 = FindCandidates(ic, w0-shift, w1-shift, false, nclusts_inwindow_shift2);
  
        
        // --------------------------------------------
        // Evaluate standard candidates
        // ---------------------------------------------
        if( v_cands.size() && v_cands.size() <= fMaxCandidates && nclusts_inwindow <=  fMaxClustMult ) {
          // plot locations
          h_zy_bipos->Fill( blip_z[iBlip], blip_y[iBlip] );
          h_wt_bipos->Fill( clust_wire[ic], clust_time[ic] );
          // loop through and update histograms/counts
          float weight = 1./v_cands.size();
          for(auto& thisCand : v_cands ) {
            _clustAvailable[thisCand.id1] = false;
            _clustAvailable[thisCand.id2] = false;
            int blipid = clust_blipid[thisCand.id2];
            if( blipid >= 0 ) _blipAvailable[blipid] = false;
            if( thisCand.dT > 6.25 && thisCand.dT < 312.5 ) {
              _numBiPo_6_312 += weight;
              h_timeFine_vs_rate  ->Fill(eventHour,weight);
            }
            if( thisCand.dT >= fdT_min ) {
              _numBiPo += weight;
              h_beta_charge       ->Fill(thisCand.q1,weight);
              h_alpha_charge      ->Fill(thisCand.q2,weight);
              h_cand_dT           ->Fill(thisCand.dT,weight);
              h_time_vs_dT        ->Fill(eventHour,thisCand.dT,weight);
            }
          }//endloop over candidates
        }//end evaluation of standard cands
        
  
        // --------------------------------------------
        // Evaluate dT-flip candidates
        // ---------------------------------------------
        if( v_cands_flip.size() && v_cands_flip.size() <= fMaxCandidates && nclusts_inwindow_flip <=  fMaxClustMult ) {
          // loop through and update histograms/counts
          float weight = 1./v_cands_flip.size();
          for(auto& thisCand : v_cands_flip ) {
            if( thisCand.dT >= fdT_min ) {
              h_beta_charge_flip  ->Fill(thisCand.q1,weight);
              h_alpha_charge_flip ->Fill(thisCand.q2,weight);
              h_cand_dT_flip      ->Fill(thisCand.dT,weight);
              h_time_vs_dT_flip   ->Fill(eventHour,thisCand.dT,weight);
            }
          }//endloop over candidates
        }//end evaluation of dT-flip candidates
        
  
        // --------------------------------------------
        // Evaluate dT-shift candidates
        // ---------------------------------------------
        if( v_cands_shift.size() && v_cands_shift.size() <= fMaxCandidates && nclusts_inwindow_shift <=  fMaxClustMult ) {
          float weight = 0.5/v_cands_shift.size();
          for(auto& thisCand : v_cands_shift ) {
            if( thisCand.dT >= fdT_min ) {
              h_beta_charge_shift ->Fill(thisCand.q1,weight);
              h_alpha_charge_shift->Fill(thisCand.q2,weight);
              h_cand_dT_shift     ->Fill(thisCand.dT, weight);
              h_time_vs_dT_shift  ->Fill(eventHour,thisCand.dT, weight);
            }
          }//endloop over candidates
        }//end evaluation of wire-shift candidates
        
        if( v_cands_shift2.size() && v_cands_shift2.size() <= fMaxCandidates && nclusts_inwindow_shift2 <=  fMaxClustMult ) {
          float weight = 0.5/v_cands_shift2.size();
          for(auto& thisCand : v_cands_shift2 ) {
            if( thisCand.dT >= fdT_min ) {
              h_beta_charge_shift ->Fill(thisCand.q1,weight);
              h_alpha_charge_shift->Fill(thisCand.q2,weight);
              h_cand_dT_shift     ->Fill(thisCand.dT, weight);
              h_time_vs_dT_shift  ->Fill(eventHour,thisCand.dT, weight);
            }
          }//endloop over candidates
        }//end evaluation of wire-shift candidates
 

  
      }//end loop over 3D blips      
  
    }//endloop over events
    
  
    
    // ****************************************************
    // Record the time it took to complete the loop
    // ****************************************************
    double loopDuration = ( time(NULL) - loopStart );
    
    h_cand_dT->Sumw2();
    h_cand_dT_flip->Sumw2();
    h_cand_dT_shift->Sumw2();

    // ***************************************************
    // Scale dT plots so they're per readout
    // ***************************************************
    _totalLiveTime = float(_numEvents) * _liveTimePerEvt;
    float scaleFact = 1./( float(_numEvents) );
    h_clust_dT        ->Scale( scaleFact );
    h_cand_dT         ->Scale( scaleFact );
    h_cand_dT_flip    ->Scale( scaleFact );
    h_cand_dT_shift   ->Scale( scaleFact );
    h_timeFine_vs_rate->Divide( h_timeFine_vs_evts );


    // ***************************************************
    // Histogram subtraction time!
    // ***************************************************
    h_cand_dT_sub     ->Add(h_cand_dT,      1);
    h_time_vs_dT_sub  ->Add(h_time_vs_dT,   1);
    h_alpha_charge_sub->Add(h_alpha_charge, 1);
    h_beta_charge_sub->Add(h_beta_charge, 1);
    
    switch(fSubtractionMode) {
      case 1: 
        h_cand_dT_sub     ->Add(h_cand_dT_flip,     -1);
        h_time_vs_dT_sub  ->Add(h_time_vs_dT_flip,  -1);
        h_alpha_charge_sub->Add(h_alpha_charge_flip,-1);
        h_beta_charge_sub->Add(h_beta_charge_flip,-1);
        break;
      case 2:
        h_cand_dT_sub     ->Add(h_cand_dT_shift,     -1);
        h_time_vs_dT_sub  ->Add(h_time_vs_dT_shift,  -1);
        h_alpha_charge_sub->Add(h_alpha_charge_shift,-1);
        h_beta_charge_sub->Add(h_beta_charge_shift,-1);
        break;
    }
    

   
    // ***************************************************
    // Check for noisy wires if doing wire diagnostics
    // ***************************************************
    float noiseThresh = 0.15;
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
      }
      std::cout<<"\n--> Found "<<nNoisy<<" wires with ave clusters/evt > "<<noiseThresh<<"\n";
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
    printf("File                : %s\n",      _fileName.c_str()); 
    printf("Total events        : %i\n",      _numEvents); 
    printf("Total live time     : %f sec\n",  _totalLiveTime);
    printf("Live time per evt   : %f us\n",   _liveTimePerEvt*1e6);
    printf("dT min/max          : %.2f-%.2f us\n",  fdT_min,fdT_max);
    printf("Ave cands / evt     : %f (%f in 6.25-312.us)\n",h_cand_dT->GetEntries()/(float)_numEvents, _numBiPo_6_312/(float)_numEvents );
    printf("Processing time     : %f sec/evt\n", loopDuration/float(_numEvents));
    printf("Excluding %i noisy wires \n",     (int)fNoisyWires.size());	
    printf("\n*******************************************\n\n");
    fOutFile->Close();
  
  }
  
  
  
  //#################################################################################
  // Make plots here
  //#################################################################################
  void makePlots()
  {
    
    // ============================================
    // Do final fit on dT spectrum
    // ============================================
    fitdT( h_cand_dT_sub, true );
   
  
    // ============================================
    // Rate vs time plots
    // ============================================
    std::cout<<"\nMaking rate vs time plot...\n";
    TH1D* h_slice;
    for(int i=1; i<=(int)h_time_vs_N->GetXaxis()->GetNbins(); i++){
      h_slice = Make1DSlice( h_time_vs_dT_sub, i, i, Form("dTfit_%i",i) );
      //float N = h_time_vs_N->GetBinContent(i);
      h_slice->Scale( 1. / (float)h_time_vs_N->GetBinContent(i) );
      std::cout<<"time: "<<h_time_vs_N->GetXaxis()->GetBinCenter(i)<<" hrs\n";
      FitResult fr = fitdT(h_slice,false);
      if( fr.rate_signal > 0 ) {
        h_time_vs_rate_bipo ->SetBinContent(i,fr.rate_signal);
        h_time_vs_rate_bg   ->SetBinContent(i,fr.rate_bg);
        h_time_vs_rate_sum  ->SetBinContent(i,fr.rate_bg+fr.rate_signal);
        h_time_vs_ratio     ->SetBinContent(i,fr.ratio);
      }
    }
    
    fOutFile->cd();
    h_time_vs_rate_bipo ->Write(0, TObject::kOverwrite );
    h_time_vs_rate_bg   ->Write(0, TObject::kOverwrite );
    h_time_vs_rate_sum  ->Write(0, TObject::kOverwrite );
    h_time_vs_ratio     ->Write(0, TObject::kOverwrite );
  
    TH1D* h1  = h_time_vs_rate_bipo; //(TH1D*)h_time_vs_rate_bipo->Clone("bipo");
    TH1D* h2  = (TH1D*)h_time_vs_rate_bg->Clone("bg");
    //TH1D* h3  = (TH1D*)h_time_vs_rate_sum->Clone("sum");
    FormatTH1D(h1, kBlue, 1, 3 );
    FormatTH1D(h2, kRed, 1, 3 );
    //FormatTH1D(h3, kBlack, 1, 3 );
    
    std::string name = "c_time_vs_rate";
    TCanvas* c = new TCanvas(name.c_str(),name.c_str(),600,500);
    //float max = std::max( GetHistMax(h1), GetHistMax(h2) );
    float max = std::max(GetHistMax(h1),GetHistMax(h2));
    h1->GetYaxis()->SetRangeUser(0,max*1.2);
    h1->DrawCopy();
    h2->DrawCopy("same");
    //h3->DrawCopy("same");
    fOutFile_plots->cd();
    c->Write();
    
  } 
  
  
  
  //################################################################################
  // Function that performs the dT fit
  //#################################################################################
  FitResult fitdT(TH1D* h, bool writeCanvas = false ){
  
    FitResult out;
  
    std::cout<<"Fitting dT spectrum "<<h->GetTitle()<<", "<<h->GetEntries()<<"\n";
    std::string label = h->GetName();
    TCanvas* c = new TCanvas(Form("c_fit_%s",label.c_str()),Form("c_fit_%s",label.c_str()),600,500);
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
    fit->SetParameter(0, histMax/5 );
    fit->SetParLimits(0, 0, histMax );
    fit->SetParameter(1, histMax );
    fit->SetParLimits(1, 0, histMax*2 );
    fit->FixParameter(2, 164.5 );

    // Draw plot and fit
    c->cd();
    hc->Fit(fit,"QR"); 
    hc->DrawCopy();
    fit->DrawCopy("same");
    gStyle->SetOptFit(1112);
    //gPad->SetLogy();
    if( writeCanvas ) {
      fOutFile_plots->cd();
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
    float n_bipo  = (f_bipo->GetParameter(0)*f_bipo->GetParameter(1))/fdT_binSize;
    
    // Components within the dT selection window
    float N_total = hc      ->Integral(0,hc->GetXaxis()->GetNbins());
    float N_fit   = fit     ->Integral(fdT_min,fdT_max)/fdT_binSize;
    float N_flat  = f_flat  ->Integral(fdT_min,fdT_max)/fdT_binSize;
    //float N_expbg = f_expBG ->Integral(fdT_min,fdT_max)/fdT_binSize;
    float N_bipo  = f_bipo  ->Integral(fdT_min,fdT_max)/fdT_binSize;
    float N_bg    = N_flat; //N_expbg;
    float sbratio = N_bipo / N_bg;
    
    printf("================ dT fit =================\n");
    printf("p0                  : %f +/- %f\n", fit->GetParameter(0),fit->GetParError(0));
    printf("p1                  : %f +/- %f\n", fit->GetParameter(1),fit->GetParError(1));
    printf("p2                  : %f +/- %f\n", fit->GetParameter(2),fit->GetParError(2));
    printf("Chi2/ndf            : %f\n",        fit->GetChisquare()/fit->GetNDF());
    printf("Total entries       : %f\n",        hc->GetEntries());
    printf("Total rate          : %f per evt\n",N_total);
    printf(" - BiPo component   : %f per evt\n",N_bipo);
    //printf(" - ExpBG component  : %f per evt\n",N_expbg);
    printf(" - Flat component   : %f per evt\n",N_flat);
    printf("S/BG ratio          : %f \n",sbratio);
    printf("\n");
    printf("Fiducial-corrected  : %f BiPos per evt\n\n",n_bipo*fFidCorFactor);
  
    out.rate_signal = n_bipo;
    out.rate_bg     = N_bg;
    out.ratio       = sbratio;
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
        if( clust_lhit_gof[jc] < fMinGOF ) continue;
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
          
              BiPoCandidate c = { ic, jc, fabs(dT), clust_charge[ic], clust_charge[jc]};
              v.push_back(c);
            }
          
         
        } 
  
      }//<-- end loop over clusters on this wire
    }//endloop over wires
  
    return v;
  
  }
  
  
