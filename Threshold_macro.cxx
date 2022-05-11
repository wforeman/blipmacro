//////////////////////////////////////////////////////////////////
// 
//  Analysis ROOT Macro
//
////////////////////////////////////////////////////////////////// 

#include "core/vars.h"
#include "core/tools.h"

// -- General macro parameters --
std::string   fFileName     = "files/BlipAna_Eminus_0to5MeV_Overlay_Standard_20220422.root";
std::string   fTreeName     = "blipanaCosFilt/anatree";

// --- Histograms ---
TH1D*   h_thresh_true_energy;
TH1D*   h_thresh_reco[3];
TH1D*   h_thresh_reco3D;
TH1D*   h_thresh[3];
TH1D*   h_thresh3D;
TH2D*   h_blip_zy;
TH2D*   h_blip_offset;
TH1D*   h_blip_xres;

// ROOT objects
TTree*  fTree;
TFile*  fOutFile;

std::string _outName = "output/plots_threshold_" + fFileName;
fOutFile = new TFile(_outFileName.c_str(), "recreate");

// Macro Functions
void    makeHistograms();
void    makePlots();
float   GetThreshold(TGraph*,float);

//#################################################################################
void makeHistograms(){
  
  // Blip spatial offset plots
  float Zmin = -100;  float Zmax = 1100;  int Zbins = 300;
  float Ymin = -150;  float Ymax = 150;   int Ybins = 75;
  h_blip_zy     = new TH2D("blip_zy","Reconstructed blip location;Z [cm]; Y [cm]", Zbins,Zmin,Zmax, Ybins, Ymin, Ymax);
  h_blip_zy     ->SetOption("colz");
  h_blip_offset = new TH2D("blip_offset","Blip spatial offset;Z offset [cm];Y offset [cm]",200,-10,10,200,-10,10);
  h_blip_offset ->SetOption("colz");
  
  h_blip_xres   = new TH1D("blip_xres","Blip X-offset;X offset [cm]",200,-10,10);

  // Energy threshold plots
  float emax = 1.5; // MeV
  int   ebins = 30; 
  h_thresh_true_energy  = new TH1D("true_energy","True blip energy;True electron energy dep [MeV];",ebins,0,emax);
  h_thresh_reco3D       = new TH1D("reco3D","Reconstructed 3D blips;True electron energy dep [MeV]",ebins,0,emax);
  h_thresh3D            = new TH1D("3Dthresh","Reconstructed threshold 3D;True electron energy dep [MeV]",ebins,0,emax);
  for(int i=0; i<kNplanes; i++){
    h_thresh_reco[i] = new TH1D(Form("pl%i_reco",i),Form("Reconstructed clusters, plane %i;True electron energy dep [MeV]",i),ebins,0,emax);
    h_thresh[i] = new TH1D(Form("pl%i_threshold",i),Form("Reconstruction threshold, plane %i;True electron energy dep [MeV]",i),ebins,0,emax);
  }
}
  


//#################################################################################
void Threshold_macro(){
  
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  
  // set branches
  fTree->SetBranchAddress("event",&event);                         
  fTree->SetBranchAddress("total_depEnergy",&total_depEnergy);     
  fTree->SetBranchAddress("nparticles",&nparticles);               
  fTree->SetBranchAddress("depEnergy",&depEnergy);                  
  //TBranch *br = 0; fTree->SetBranchAddress("process",&process,&br);
  fTree->SetBranchAddress("nedeps",&nedeps);                       
  fTree->SetBranchAddress("edep_g4id",&edep_g4id);                  
  fTree->SetBranchAddress("edep_blipid",&edep_blipid);              
  fTree->SetBranchAddress("edep_clustid",&edep_clustid);            
  fTree->SetBranchAddress("edep_energy",&edep_energy);              
  fTree->SetBranchAddress("edep_charge",&edep_charge);              
  fTree->SetBranchAddress("edep_x",&edep_x);                        
  fTree->SetBranchAddress("edep_y",&edep_y);                        
  fTree->SetBranchAddress("edep_z",&edep_z);                        
  fTree->SetBranchAddress("nhits",&nhits);
  /*
  fTree->SetBranchAddress("hit_plane",&hit_plane);                  
  fTree->SetBranchAddress("hit_wire",&hit_wire);                    
  fTree->SetBranchAddress("hit_channel",&hit_channel);              
  fTree->SetBranchAddress("hit_peakT",&hit_peakT);                  
  fTree->SetBranchAddress("hit_time",&hit_time);                    
  fTree->SetBranchAddress("hit_rms",&hit_rms);                      
  fTree->SetBranchAddress("hit_ph",&hit_amp);                        
  fTree->SetBranchAddress("hit_area",&hit_area);                    
  fTree->SetBranchAddress("hit_charge",&hit_charge);                
  fTree->SetBranchAddress("hit_isreal",&hit_isreal);                
  fTree->SetBranchAddress("hit_trkid",&hit_trkid);                  
  fTree->SetBranchAddress("hit_g4id",&hit_g4id);                    
  fTree->SetBranchAddress("hit_g4frac",&hit_g4frac);                
  fTree->SetBranchAddress("hit_g4energy",&hit_g4energy);            
  fTree->SetBranchAddress("hit_g4charge",&hit_g4charge);            
  fTree->SetBranchAddress("hit_clustid",&hit_clustid);              
  */
  fTree->SetBranchAddress("nclusts",&nclusts);                     
  fTree->SetBranchAddress("clust_ismatch",&clust_ismatch);
  fTree->SetBranchAddress("clust_plane",&clust_plane);              
  fTree->SetBranchAddress("clust_wire",&clust_wire);                
  fTree->SetBranchAddress("clust_charge",&clust_charge);            
  fTree->SetBranchAddress("clust_time",&clust_time);                
  fTree->SetBranchAddress("clust_g4charge",&clust_g4charge);        
  fTree->SetBranchAddress("clust_g4energy",&clust_g4energy);        
  fTree->SetBranchAddress("clust_edepid",&clust_edepid);            
  fTree->SetBranchAddress("clust_blipid",&clust_blipid);            
  fTree->SetBranchAddress("nblips",&nblips);                       
  fTree->SetBranchAddress("blip_tpc",&blip_tpc);                    
  fTree->SetBranchAddress("blip_x",&blip_x);                        
  fTree->SetBranchAddress("blip_y",&blip_y);                        
  fTree->SetBranchAddress("blip_z",&blip_z);                        
  fTree->SetBranchAddress("blip_charge",&blip_charge);              
  fTree->SetBranchAddress("blip_edepid",&blip_edepid);              
  
  // make output file to store plots
  fOutFile = new TFile(fOutName.c_str(),"recreate");
  
  // initialize histograms
  makeHistograms();
  
  // Configure histograms and TFile
  configure();

  // Loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
      
    //int sparsify  = 10; if(  (iEvent % sparsify) != 0 ) continue; 
    
    fTree->GetEntry(iEvent);
    std::cout<<"========== EVENT "<<iEvent<<" ========================\n";

    // ===========================================
    // Event-wide values
    float true_charge           = 0;
    float true_charge_75keV     = 0;
    float true_charge_150keV     = 0;
    float true_charge_300keV    = 0;
    float true_charge_perf      = 0;  
    float reco_charge_pl[3]     = {0};
    int   numhits[kNplanes][2]  = {0};


    // ===========================================
    // Look at truth information first
    //std::cout<<"MCParticles: "<<nparticles<<"\n";
    h_totalVisE->Fill(total_depEnergy);
    for(int i=0; i<nparticles; i++){
      // ignore dummy entries, nuclear fragments, and neutrons
      //if( pdg[i] < -9999 || pdg[i] > 100000 || pdg[i] == 2112 ) continue;
      //printG4Particle(i);
    }
    

    // ===========================================
    // Look at true blips ("edeps")
    //std::cout<<"True energy depositions: "<<nedeps<<"\n";
    for(int i=0; i<nedeps; i++){
      //std::cout<<"  "<<i<<"  E= "<<edep_energy[i]<<"   X= "<<edep_x[i]<<"\n";
      h_thresh_true_energy->Fill(edep_energy[i]);
      h_edep_charge->Fill(edep_charge[i]/1e3);
      true_charge += edep_charge[i];
      if(edep_energy[i]>0.075) true_charge_75keV += edep_charge[i];
      if(edep_energy[i]>0.150) true_charge_150keV += edep_charge[i];
      if(edep_energy[i]>0.300) true_charge_300keV += edep_charge[i];
      if(edep_blipid[i]>=0)    true_charge_perf   += edep_charge[i];
    }


    // ===========================================
    // Loop over the hits
    /*
    for(int i=0; i<nhits;i++){
      int plane = hit_plane[i];
      int isReal // = (int)hit_isreal[i];
                = (int)(hit_g4id[i] >= 0 );
      numhits[plane][isReal]++;
      h_hitamp[plane][isReal]->Fill(hit_amp[i]);
      h_hitrms[plane][isReal]->Fill(hit_rms[i]);
      h_hitratio[plane][isReal]->Fill(hit_amp[i]/hit_rms[i]);
      h_hitrms_vs_amp[plane][isReal]->Fill(hit_rms[i],hit_amp[i]);
      if( isReal ) {
        h_nelec_TrueVsReco[plane]->Fill(hit_charge[i]/1e3,hit_g4charge[i]/1e3);
        h_nelec_Resolution[plane]->Fill((hit_charge[i]-hit_g4charge[i])/hit_g4charge[i]);
      }
    }
    */

  
    // ==========================================
    // Loop over the reconstructed hitclusts and blips
    for(int i=0; i<nclusts; i++){
      int plane = clust_plane[i];
      int blipid = clust_blipid[i];
      int edepid = clust_edepid[i];
      float x = clust_time[i];
      //std::cout<<"Cluster "<<i<<" on plane "<<plane<<", nwires = "<<clust_nwires[i]<<", maps to edepid "<<edepid<<"   X= "<<clust_xpos[i]<<"\n";
      if( edepid >= 0 ) {
        h_thresh_reco[plane]  -> Fill(edep_energy[edepid]);
        //std::cout<<"   true energy: "<<edep_energy[edepid]*1000<<"\n";
      }

      // Look at time separation between planes
      /*
      for(int j=i+1; j<nclusts; j++){
        int plane_j = clust_plane[j];
        float Ti = clust_time[i];
        float Tj = clust_time[j];
        float wid_i = (clust_endTime[i]-clust_startTime[i])/2.;
        float wid_j = (clust_endTime[j]-clust_startTime[j])/2.;
        if( plane != plane_j ) {
          h_clust_dT->Fill( Tj - Ti );
          h_clust_dTfrac->Fill( (Tj-Ti)/std::max(wid_i,wid_j) );
        }
      }
      */

      // If a clust is considered good if both:
      //  (1) matched to at least one other clust in time
      //  (2) part of a blip object, meaning at least one valid wire-crossing was found
      bool requireMatch = true;
      bool requireBlip  = true;
      if( (!requireMatch || clust_ismatch[i]) && (!requireBlip || clust_blipid[i] >=0 )) {
        reco_charge_pl[plane] += clust_charge[i];
        //std::cout<<"  cluster "<<i<<" on plane "<<plane<<" at T= "<<clust_time[i]
        //<<" has charge "<<clust_charge[i]<<"   "<<reco_charge_pl[plane]<<", "<<clust_nwires[i]<<" wires, and is mapped to edep "<<edepid<<"\n";
      }

    }//endloop over clusts
  
    // Loop over BLIPS 
    for(int i=0; i<nblips; i++){
      h_blip_zy->Fill(blip_z[i],blip_y[i]);
      int eid = blip_edepid[i];
      if( eid >= 0 ){
        h_thresh_reco3D-> Fill(edep_energy[eid]);
      }
      float dx = blip_x[i]-edep_x[eid];
      float dy = blip_y[i]-edep_y[eid];
      float dz = blip_z[i]-edep_z[eid];
      h_blip_xres->Fill(dx);
      h_blip_offset->Fill(dz,dy);
    }


    // =========================================
    // Fill event-wide histograms
    h_true_charge->Fill(true_charge/1e3);
    h_true_charge_75keV ->Fill(true_charge_75keV/1e3);
    h_true_charge_150keV->Fill(true_charge_150keV/1e3);
    h_true_charge_300keV->Fill(true_charge_300keV/1e3);
    //h_true_charge_perf  ->Fill(true_charge_perf/1e3);
    for(int iPl=0; iPl<kNplanes; iPl++){
      //std::cout<<"Reco charge plane "<<iPl<<" "<<reco_charge_pl[iPl]<<"\n";
      if( reco_charge_pl[iPl] > 0 ) {
        h_reco_charge_plane[iPl]  ->Fill(reco_charge_pl[iPl]/1e3);
        if( iPl==2 ) h_reco_charge->Fill(reco_charge_pl[iPl]/1e3);
      }

      for(int iReal=0; iReal<2; iReal++){
        //std::cout<<"plane "<<iPl<<"   isReal "<<iReal<<"  numhits "<<numhits[iPl][iReal]<<"\n";
        h_nhits[iPl][iReal]->Fill(numhits[iPl][iReal]);
      }
    }

  }//endloop over events
  

  makePlots();


  fOutFile->Write(); 
  fOutFile->Close();

}






//#################################################################################
void makePlots(){
  
  bool makeThresholdPlot  = true;
  bool makeHitNumPlot     = true;
  bool makeHitMetricPlot  = true;

  if(makeThresholdPlot){
    // Make threshold plot
    TGraphErrors* gr_thresh[3];
    TGraphErrors* gr_thresh3D;
    TLegend* leg_thresh;

    size_t Nbins = h_thresh_true_energy->GetXaxis()->GetNbins();

    // Hitclust efficiency (2D)
    for(int i=0; i<kNplanes;i++) {
      gr_thresh[i] = new TGraphErrors();
      for(size_t j=0; j<Nbins; j++){
        float num = h_thresh_reco[i]->GetBinContent(j);
        float denom = h_thresh_true_energy->GetBinContent(j);
        float ratio = 0;
        if( denom > 0 ) ratio = num/float(denom);
        if( denom > 0 && ratio < 1 ) {
          h_thresh[i]->SetBinContent(j,ratio);
          gr_thresh[i]->SetPoint(gr_thresh[i]->GetN(),h_thresh_reco[i]->GetBinCenter(j),ratio);
        }
      }
      std::cout<<"Threshold graph for plane "<<i<<" has "<<gr_thresh[i]->GetN()<<" points.\n";
    }

    // Blip efficiency (3D)
    gr_thresh3D = new TGraphErrors();  
    for(size_t j=0; j<Nbins; j++){
      float num = h_thresh_reco3D->GetBinContent(j);
      float denom = h_thresh_true_energy->GetBinContent(j);
      float ratio = 0;
      if( denom > 0 ) ratio = num/float(denom);
      if( denom > 0 && ratio < 1 ) {
        h_thresh3D->SetBinContent(j,ratio);
        gr_thresh3D->SetPoint(gr_thresh3D->GetN(),h_thresh_reco3D->GetBinCenter(j),ratio);
      }
    }
    
    
    gr_thresh[2]->GetXaxis()->SetTitle("Energy Deposited by Electron [MeV]");
    gr_thresh[2]->GetYaxis()->SetTitle("Reconstruction Efficiency");
    gr_thresh[2]->GetXaxis()->SetRangeUser(0,h_thresh_true_energy->GetXaxis()->GetXmax());
    FormatTGraph(gr_thresh[2],kBlue,kBlue,20,1,0.7,1);
    FormatTGraph(gr_thresh[1],kGreen+2,kGreen+2,  22,1,0.7,1);
    FormatTGraph(gr_thresh[0],kRed,kRed,    21,1,0.7,1);
    FormatTGraph(gr_thresh3D,kBlack,kBlack, 4,2,0.7,2);

    std::cout<<"A\n";

    TCanvas* c_th = new TCanvas("threshold","threshold",500,500);
    std::cout<<"Made canvas\n";
    gr_thresh[2]->Draw("APL");
    gr_thresh[1]->Draw("PL same");
    gr_thresh[0]->Draw("PL same");
    gr_thresh3D->Draw("PL same");
    
    leg_thresh = MakeLegend(0.5,0.3,0.035,4,0.25);
    leg_thresh->AddEntry(gr_thresh[2],"Collection Plane","PL");
    leg_thresh->AddEntry(gr_thresh[1],"Induction V Plane","PL");
    leg_thresh->AddEntry(gr_thresh[0],"Induction U Plane","PL");
    leg_thresh->AddEntry(gr_thresh3D, "3D Blip Reconstruction","PL");
    leg_thresh->Draw();
    
    std::cout<<"B\n";

    c_th->Write();
   
    // Print out 50% crossings
    float effthresh;

    effthresh = 0.2;
    for(size_t i=0; i<kNplanes;i++){
      float thresh = GetThreshold(gr_thresh[i],effthresh);
      std::cout<<effthresh*100.<<"% threshold crossing, plane "<<i<<" --> "<<thresh<<" keV.\n";
    }
    std::cout<<100*effthresh<<"% threshold crossing for 3D --> "<<GetThreshold(gr_thresh3D,effthresh)<<" keV.\n\n";
    
    effthresh = 0.5;
    for(size_t i=0; i<kNplanes;i++){
      float thresh = GetThreshold(gr_thresh[i],effthresh);
      std::cout<<effthresh*100.<<"% threshold crossing, plane "<<i<<" --> "<<thresh<<" keV.\n";
    }
    std::cout<<effthresh*100.<<"% threshold crossing for 3D --> "<<GetThreshold(gr_thresh3D,effthresh)<<" keV.\n\n";
 
    
    effthresh = 0.8;
    for(size_t i=0; i<kNplanes;i++){
      float thresh = GetThreshold(gr_thresh[i],effthresh);
      std::cout<<effthresh*100.<<"% threshold crossing, plane "<<i<<" --> "<<thresh<<" keV.\n";
    }
    std::cout<<effthresh*100.<<"% threshold crossing for 3D --> "<<GetThreshold(gr_thresh3D,effthresh)<<" keV.\n\n";

  }

  //******************************************************    
  // Make noise hit number plots
  if(makeHitNumPlot) {
    
    TCanvas* c_num = new TCanvas("hits","hits",500,500);
    TLegend* leg;
    TH1D* h[3];
    for(int i=0; i<3; i++) h[i] = Clone(h_nhits[i][0]);
    float max = std::max(GetHistMax(h[0]),GetHistMax(h[1]));
    max = std::max(max,GetHistMax(h[2]));
    h[2]->SetMaximum(max*1.1);
    h[2]->GetXaxis()->SetRangeUser(0,150);
    h[2]->GetXaxis()->SetTitle("Number of noise hits per readout");
    FormatTH1D(h[2],kBlack,1,2);
    FormatTH1D(h[1],kBlue,1,2);
    FormatTH1D(h[0],kRed,1,2);
    h[2]->DrawCopy("hist");
    h[1]->DrawCopy("hist same");
    h[0]->DrawCopy("hist same");
   
    leg = MakeLegend(0.5,0.8,0.035,3,0.25);
    leg->AddEntry(h[2],"Collection Plane","L");
    leg->AddEntry(h[1],"Induction V Plane","L");
    leg->AddEntry(h[0],"Induction U Plane","L");
    leg->Draw();
  }

  
  //*****************************************************
  // Make hit metric plots


  float mar_l = 0.1;
  float mar_r = 0.03;
  float mar_b = 0.1;
  float mar_t = 0.1;
  
  // Hit widths
  TCanvas* c_hw = new TCanvas("hit_widths","hit_widths",1200,400);
  c_hw->Divide(3,1);
  for(int i=0; i<kNplanes;i++) {
    c_hw->cd(i+1);
    gPad->SetMargin(mar_l,mar_r,mar_b,mar_t);
    gStyle->SetOptStat(0);
    ScaleHist(h_hitrms[i][0],"integral");
    ScaleHist(h_hitrms[i][1],"integral");
    TH1D* h0 = Clone(h_hitrms[i][0]);
    TH1D* h1 = Clone(h_hitrms[i][1]);

    float max = std::max(GetHistMax(h0),GetHistMax(h1));
    max *= 1.1;
    h0->SetLineWidth(2);
    FormatAxes(h0, 0.05, 0.045, 1.0, 1.4);
    CopyHistoFormat(h0,h1);
    h0->SetMaximum(max);
    h0->SetTitle(Form("Plane %i",i));
    h0->SetLineColor(kRed);
    h1->SetLineColor(kBlack);
    h0->DrawCopy("hist");
    h1->DrawCopy("hist same");
  }
  c_hw->Write();
  
  // Hit width vs PH scatterplots
  TCanvas* c_hwamp = new TCanvas("hit_widths_vs_amp","hit_widths_vs_amp",1200,400);
  c_hwamp->Divide(3,1);
  for(int i=0; i<kNplanes;i++) {
    c_hwamp->cd(i+1);
    gPad->SetMargin(0.15,0.03,0.1,0.1);
    gStyle->SetOptStat(0);
    TH2D* h0 = Clone(h_hitrms_vs_amp[i][0]);
    TH2D* h1 = Clone(h_hitrms_vs_amp[i][1]);
    FormatAxes(h0, 0.05, 0.045, 1.0, 1.4);
    CopyHistoFormat(h0,h1);
    h0->SetTitle(Form("Plane %i",i));
    h0->SetMarkerColorAlpha(kRed,1.0);
    h1->SetMarkerColorAlpha(kBlack,1.0);
    h0->SetMarkerStyle(20);
    h0->SetMarkerSize(0.1);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(0.1);
    h0->DrawCopy("hist");
    h1->DrawCopy("hist same");
  }
}

//#######################################################
float GetThreshold(TGraph* g, float thresh){
  for(size_t e=0; e<5000; e+=1){
    float eval = g->Eval(e*1e-3);
    if( eval >= thresh ) return e;
  }

  return -1;

}
