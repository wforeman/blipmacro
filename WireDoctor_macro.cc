  //////////////////////////////////////////////////////////////////
  // 
  // Diagnose which wires are noisy
  //
  ////////////////////////////////////////////////////////////////// 
  
  #include "core/vars.h"
  #include "core/tools.h"
  
  std::string fFileName     = "BlipAna_RadonData_Phase1_20220710.root";
  std::string fTreeName     = "blipanaTrkMask/anatree";
  std::string fOutFileName  = "output/plots_wirediagnostics_" + fFileName;
  
  float noiseThresh = 0.15;
  
  // Detector properties
  int   nWiresColl        = 3455;
  
  // Counters, maps, etc
  int _numEvents      = 0;
  std::map<int,std::vector<int>> _map_wire_clusters;
  
  // Histograms
  TH1D* h_nclusts_wire_ave;
  TH1D* h_clust_dT;
  TH1D* h_wire_vs_aveN;
  
  
  //#################################################################################
  // Primary macro
  //#################################################################################
  void WireDoctor_macro()
  {
    // ****************************************************
    // Initial config
    // ****************************************************
    std::cout<<"Reading input file "<<fFileName<<"\n";
    TFile* file = new TFile(("files/"+fFileName).c_str(),"READ");
    TTree* fTree = (TTree*)file->Get(fTreeName.c_str());
  
    //fTree->SetBranchAddress("timestamp",&timestamp);
    fTree->SetBranchAddress("nclusts",&nclusts);                     
    fTree->SetBranchAddress("clust_plane",&clust_plane);              
    fTree->SetBranchAddress("clust_wire",&clust_wire);                
    fTree->SetBranchAddress("clust_startwire",&clust_startwire);                
    fTree->SetBranchAddress("clust_endwire",&clust_endwire);                

    TFile* fOutFile = new TFile(fOutFileName.c_str(), "recreate");
    
    // ****************************************************
    // Make histograms
    // ****************************************************
    h_nclusts_wire_ave= new TH1D("nclusts_perwire","Average clusters per wire;Average number of clusters per wire",1500,0,3);
    h_wire_vs_aveN = new TH1D("wire_vs_aveN",";Collection plane wire;Average number of hits",nWiresColl,0,nWiresColl);
 
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
      }
    
    }//endloop on events
  

    // ***************************************************
    // Check for noisy wires if doing wire diagnostics
    // ***************************************************
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
  

    // ***************************************************
    // Write all histos currently in stack
    // ***************************************************
    fOutFile->Write(); 
 

  }
  
