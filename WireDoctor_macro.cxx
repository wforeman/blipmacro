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
  
  std::string fFileName[3] =  {    
                                //"files/BlipAna_BiPo_MC_Overlay_AlphaQY_20220422.root",
                                "files/BlipAna_BiPo_MC_Overlay_AlphaQY_20220422.root",
                                "files/BlipAna_RadonData_Phase1_20220422.root",
                                "files/BlipAna_RadonData_Phase2_OldSettings_20220502.root"
//                                "files/BlipAna_RadonDAta_Phase2_20220422.root"
                              };
  

  // --- Output histogram file ---
  std::string fOutFileName  = "output/plots_wire.root";
  std::string fTreeName = "blipanaTrkMask/anatree";
  
  // --- Detector properties ---
  int   nWiresColl        = 3455;
  float fSamplePeriod     = 0.5; // microseconds
  
  // --- Special switches ---
  bool  fDoWireDiagnostics= true; //true;
  
  std::string   _fileName = fFileName[fConfig];
  
  // Counters
  int   _numEvents      = 0;
  
  // Map of cluster IDs per wire on collection plane
  std::map<int,std::vector<int>> _map_wire_clusters;
  
  // Functions 
  void                        makeHistograms();
  
  // ROOT objects
  TTree*      fTree;
  TFile*      fOutFile;

  // Histograms
  TH1D* h_nclusts_wire_ave;
  TH1D* h_clust_dT;
  
  TH1D* h_wire_vs_aveN;

  //##########################################################################
  // Initialize histograms
  //##########################################################################
  void makeHistograms()
  {
    
    h_nclusts_wire_ave= new TH1D("nclusts_perwire","Average clusters per wire;Average number of clusters per wire",1500,0,3);
    h_wire_vs_aveN = new TH1D("wire_vs_aveN",";Collection plane wire;Average number of hits",nWiresColl,0,nWiresColl);
    //h_clust_dT        = (TH1D*)h_cand_dT->Clone("clust_dT");  h_clust_dT->SetTitle("Same-wire cluster separations");
    
  }
  
  
  //#################################################################################
  // Primary macro
  //#################################################################################
  void WireDoctor_macro()
  {
    std::cout<<"Reading input file "<<_fileName<<"\n";
    TFile* file = new TFile(_fileName.c_str(),"READ");
    fTree = (TTree*)file->Get(fTreeName.c_str());
  
    fTree->SetBranchAddress("timestamp",&timestamp);
    fTree->SetBranchAddress("nclusts",&nclusts);                     
    fTree->SetBranchAddress("clust_plane",&clust_plane);              
    fTree->SetBranchAddress("clust_wire",&clust_wire);                
    fTree->SetBranchAddress("clust_startwire",&clust_startwire);                
    fTree->SetBranchAddress("clust_endwire",&clust_endwire);                
    fTree->SetBranchAddress("clust_nhits",&clust_nhits);              
    fTree->SetBranchAddress("clust_charge",&clust_charge);            
    fTree->SetBranchAddress("clust_time",&clust_time);                
    fTree->SetBranchAddress("clust_lhit_isfit",&clust_lhit_isfit);          
    fTree->SetBranchAddress("clust_timespan",&clust_timespan);          
    fTree->SetBranchAddress("clust_blipid",&clust_blipid);            

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");
    
    makeHistograms();
 
    std::map<int,int> map_wn;
    
    int print_counter = 0;

    // ****************************************************
    // Loop over the events
    // ****************************************************
    size_t totalEntries = fTree->GetEntries();
    for(size_t iEvent=0; iEvent < totalEntries; iEvent++){
      
      print_counter++;
      if( print_counter > 1000 ) {
        printf("========== EVENT %lu / %lu =====================\n",iEvent,totalEntries);
        print_counter = 1;
      }
      
      // Retrieve event info
      fTree->GetEntry(iEvent);
      _numEvents++;

      // ====================================================
      // Map of clust IDs per wire on collection plane
      // ====================================================
      _map_wire_clusters.clear();
      for(int i=0; i < nclusts; i++){
        if( clust_plane[i] != 2 ) continue;
        _map_wire_clusters[clust_wire[i]].push_back(i);
        map_wn[clust_wire[i]] += 1;
        //map_chn[clust_chan[i]] += 1;
        //h_clust_timespan->Fill(clust_timespan[i]);
        //h_wt_clusts->Fill(clust_wire[i],clust_time[i]);
      }
    
      // ====================================================
      // Wire diagnostics (time-consuming!)
      // ====================================================
      /*
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
            //h_clust_dT->Fill(v_t.at(j)-v_t.at(j-1));
        
        }
      }
      */
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
 

  }
  
