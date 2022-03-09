//////////////////////////////////////////////////////////////////
// 
//  SBND Analysis ROOT Macro
//
////////////////////////////////////////////////////////////////// 

#include "core/vars.h"
#include "core/tools.h"

// General macro parameters
//std::string   fFileName     = "files/hist_BlipAna_defaultReco_220125.root";
std::string   fFileName     = "files/hist_BlipAna_lowThreshReco_220128.root";
std::string   fOutName      = "output/plots.root";
std::string   fTreeName     = "blipana/anatree";

// Histograms
TH1D*   h_totalVisE;
TH1D*   h_edep_charge;
TH1D*   h_true_charge;
TH1D*   h_true_charge_75keV;
TH1D*   h_true_charge_150keV;
TH1D*   h_true_charge_300keV;
//TH1D*   h_true_charge_perf;
TH1D*   h_reco_charge;
TH1D*   h_reco_charge_plane[kNplanes];
TH1D*   h_charge_frac[kNplanes];    // fraction of true charge reconstructed
TH1D*   h_hit_purity[kNplanes];     // (sel #hits)/(total true #hits)
TH1D*   h_hit_comp[kNplanes];       // (sel true #hits)/(total true #hits)
TH1D*   h_nhits[kNplanes][2];       // --- hit diagnostics ---
TH1D*   h_hitph[kNplanes][2];
TH1D*   h_hitwidth[kNplanes][2];
TH1D*   h_hitphratio[kNplanes][2];
TH1D*   h_hitneighbors[kNplanes][2];
TH2D*   h_hitwidth_vs_ph[kNplanes][2];
TH2D*   h_nelec_TrueVsReco[kNplanes];
TH1D*   h_nelec_Resolution[kNplanes];
TH1D*   h_clust_dT;                 // --- interplane timing ---
TH1D*   h_clust_dTfrac;
TH1D*   h_thresh_true_energy;       // --- for threshold calc ---
TH1D*   h_thresh_reco3D;
TH1D*   h_thresh_reco[3];
TH1D*   h_thresh[3];
TH1D*   h_thresh3D;
TH2D*   h_blip_zy;
TH2D*   h_blip_offset;
TH1D*   h_blip_xres;

TDirectory* dir_diag;
TDirectory* dir_thresh;

//#################################################################################
void configure(){
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches();
  // make output file to store plots
  fOutFile = new TFile(fOutName.c_str(),"recreate");
  // initialize all histograms
  makeHistograms();
  fOutFile->cd();
}


//#################################################################################
void BlipAna_macro(){
  
  // Configure histograms and TFile
  configure();

  // Loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
     
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
    std::cout<<"MCParticles: "<<nparticles<<"\n";
    h_totalVisE->Fill(total_depEnergy);
    for(int i=0; i<nparticles; i++){
      // ignore dummy entries, nuclear fragments, and neutrons
      if( pdg[i] < -9999 || pdg[i] > 100000 || pdg[i] == 2112 ) continue;
      printG4Particle(i);
    }
    

    // ===========================================
    // Look at true blips ("edeps")
    //std::cout<<"True energy depositions: "<<nedeps<<"\n";
    for(int i=0; i<nedeps; i++){
      //std::cout<<"  "<<i<<"  E= "<<edep_energy[i]<<"   X= "<<edep_x[i]<<"\n";
      h_thresh_true_energy->Fill(edep_energy[i]*1000);
      h_edep_charge->Fill(edep_charge[i]/1e3);
      true_charge += edep_charge[i];
      if(edep_energy[i]>0.075) true_charge_75keV += edep_charge[i];
      if(edep_energy[i]>0.150) true_charge_150keV += edep_charge[i];
      if(edep_energy[i]>0.300) true_charge_300keV += edep_charge[i];
      if(edep_blipid[i]>=0)    true_charge_perf   += edep_charge[i];
    }


    // ===========================================
    // Loop over the hits 
    for(int i=0; i<nhits;i++){
      int plane = hit_plane[i];
      int isReal // = (int)hit_isreal[i];
                = (int)(hit_g4id[i] >= 0 );
      numhits[plane][isReal]++;
      h_hitph[plane][isReal]->Fill(hit_ph[i]);
      h_hitwidth[plane][isReal]->Fill(hit_rms[i]);
      h_hitphratio[plane][isReal]->Fill(hit_ph[i]/hit_rms[i]);
      h_hitwidth_vs_ph[plane][isReal]->Fill(hit_rms[i],hit_ph[i]);
      if( isReal ) {
        h_nelec_TrueVsReco[plane]->Fill(hit_charge[i]/1e3,hit_g4charge[i]/1e3);
        h_nelec_Resolution[plane]->Fill((hit_charge[i]-hit_g4charge[i])/hit_g4charge[i]);
      }
    }

  
    // ==========================================
    // Loop over the reconstructed hitclusts and blips
    for(int i=0; i<nclusts; i++){
      int plane = clust_plane[i];
      int blipid = clust_blipid[i];
      int edepid = clust_edepid[i];
      float x = clust_time[i];
      //std::cout<<"Cluster "<<i<<" on plane "<<plane<<", nwires = "<<clust_nwires[i]<<", maps to edepid "<<edepid<<"   X= "<<clust_xpos[i]<<"\n";
      if( edepid >= 0 ) {
        h_thresh_reco[plane]  -> Fill(edep_energy[edepid]*1000);
        //std::cout<<"   true energy: "<<edep_energy[edepid]*1000<<"\n";
      }

      // Look at time separation between planes
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

      // If a clust is considered good if both:
      //  (1) matched to at least one other clust in time
      //  (2) part of a blip object, meaning at least one valid wire-crossing was found
      bool requireMatch = true;
      bool requireBlip  = true;
      if( (!requireMatch || clust_ismatched[i]) && (!requireBlip || clust_blipid[i] >=0 )) {
        reco_charge_pl[plane] += clust_charge[i];
        std::cout<<"  cluster "<<i<<" on plane "<<plane<<" at T= "<<clust_time[i]
        <<" has charge "<<clust_charge[i]<<"   "<<reco_charge_pl[plane]<<", "<<clust_nwires[i]<<" wires, and is mapped to edep "<<edepid<<"\n";
      }

    }//endloop over clusts
  
    // Loop over BLIPS 
    for(int i=0; i<nblips; i++){
      h_blip_zy->Fill(blip_z[i],blip_y[i]);
      int eid = blip_edepid[i];
      if( eid >= 0 ){
        h_thresh_reco3D-> Fill(edep_energy[eid]*1000);
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
      std::cout<<"Reco charge plane "<<iPl<<" "<<reco_charge_pl[iPl]<<"\n";
      if( reco_charge_pl[iPl] > 0 ) {
        h_reco_charge_plane[iPl]  ->Fill(reco_charge_pl[iPl]/1e3);
        if( iPl==2 ) h_reco_charge->Fill(reco_charge_pl[iPl]/1e3);
      }

      for(int iReal=0; iReal<2; iReal++){
        std::cout<<"plane "<<iPl<<"   isReal "<<iReal<<"  numhits "<<numhits[iPl][iReal]<<"\n";
        h_nhits[iPl][iReal]->Fill(numhits[iPl][iReal]);
      }
    }
    
    

  }//endloop over events
  

  makePlots();


  fOutFile->Write(); 
  fOutFile->Close();

}



//#################################################################################
void makeHistograms(){
  // make the histograms
  float chargeMax = 1000;   int chargeBins= 200; // using units of x10^3 electrons
  h_totalVisE = new TH1D("totalVisE","Total deposited energy;Energy deposited [MeV];Events",  500,0,50.);
  h_edep_charge = new TH1D("edep_charge","Energy deposition charge;Charge [ #times10^{3} electrons ];Events",chargeBins,0,chargeMax);
  h_true_charge = new TH1D("true_charge","Total true charge;Charge [ #times10^{3} electrons ];Events",chargeBins,0,chargeMax);
  h_true_charge_75keV = new TH1D("true_charge_75keV","Total true charge (75 keV threshold);Charge [ #times10^{3} electrons ];Events",chargeBins,0,chargeMax);
  h_true_charge_150keV = new TH1D("true_charge_150keV","Total true charge (150 keV threshold);Charge [ #times10^{3} electrons ];Events",chargeBins,0,chargeMax);
  h_true_charge_300keV = new TH1D("true_charge_300keV","Total true charge (300 keV threshold);Charge [ #times10^{3} electrons ];Events",chargeBins,0,chargeMax);
  //h_true_charge_perfReco = new TH1D("true_charge_perfectReco","Total true charge of reconstructed blips;Charge [ #times10^{3} electrons ];Events",chargeBins,0,chargeMax);
  h_reco_charge = new TH1D("total_reco_charge","Total reconstructed charge, best plane;Charge [ #times10^{3} electrons ];Events",  chargeBins,0,chargeMax);
  h_clust_dT      = new TH1D("clust_dT","Inter-plane hit time separation;dT [ticks]",200,-20,20);
  h_clust_dTfrac  = new TH1D("clust_dTfrac","Inter-plane hit time separation; dT/width",120,-3,3);
  for(int i=0; i<kNplanes; i++) {
    h_reco_charge_plane[i]  = new TH1D(Form("pl%i_reco_charge",i),Form("Reco charge, plane %i;Charge [ #times10^{3} electrons ];Events",i),chargeBins,0,chargeMax);
    h_charge_frac[i]        = new TH1D(Form("pl%i_charge_frac",i),Form("Fraction of true charge reco'd, plane %i;Reco charge fraction;Events",i),100,0,2);
    h_hit_purity[i]         = new TH1D(Form("pl%i_hit_purity",i), Form("Hit purity, plane %i;Purity;Events",i),120,0,1.2);
    h_hit_comp[i]           = new TH1D(Form("pl%i_hit_comp",i),   Form("Hit completeness, plane %i;Completeness;Events",i),100,0,1.2);
  }
  
  // Blip spatial offset
  h_blip_zy     = new TH2D("blip_zy","Reconstructed blip location;Z [cm]; Y [cm]",70,-100,600,60,-300,300);
  h_blip_zy     ->SetOption("colz");
  h_blip_offset = new TH2D("blip_offset","Blip spatial offset;Z offset [cm];Y offset [cm]",200,-10,10,200,-10,10);
  h_blip_offset->SetOption("colz");
  
  h_blip_xres   = new TH1D("blip_xres","Blip X-offset;X offset [cm]",200,-5,5);

  // Plane hit diagnostics
  // make subdirectories for histograms
  dir_diag = fOutFile->mkdir("PlaneDiagnostics");
  float hitMax  = 150;    int hitBins   = 150;
  float phMax = 20;      int phBins    = 100;
  float widthMax = 10;    int widthBins = 200;
  float ratioMax = 30;    int ratioBins = 150;
  dir_diag->cd();
  for(int i=0; i<kNplanes; i++) {
    std::vector<std::string> label = { "noise", "real" };
    for(int j=0; j<2; j++){
      h_nhits[i][j]     = new TH1D(Form("pl%i_%s_nhits",i,label[j].c_str()),       Form("Plane %i, %s hits;Number of hits",i,label[j].c_str()),hitBins,0,hitMax);
      h_hitph[i][j]     = new TH1D(Form("pl%i_%s_hit_ph",i,label[j].c_str()),      Form("Plane %i, %s hits;Pulse height [ADC]",i,label[j].c_str()),phBins,0,phMax);
      h_hitwidth[i][j]  = new TH1D(Form("pl%i_%s_hit_width",i,label[j].c_str()),   Form("Plane %i, %s hits;Hit width [ticks]",i,label[j].c_str()),widthBins,0,widthMax);
      h_hitphratio[i][j]= new TH1D(Form("pl%i_%s_hit_phratio",i,label[j].c_str()), Form("Plane %i, %s hits;Height-to-width ratio",i,label[j].c_str()),  ratioBins,0,ratioMax);
      h_hitwidth_vs_ph[i][j]= new TH2D(Form("pl%i_%s_hit_width_vs_ph",i,label[j].c_str()), Form("Plane %i, %s hits;Hit width [ticks];Pulse height [ADC]",i,label[j].c_str()),  
        widthBins/2,0,widthMax,
        phBins/2,0,phMax);
    }
    h_nelec_TrueVsReco[i] = new TH2D( Form("pl%i_nelec_TrueVsReco",i),
      Form("Plane %i;True hit charge [ #times 10^{3} electrons ];Reconstructed hit charge [ #times 10^{3} electrons ]",i),60,0,30, 60,0,30);
    h_nelec_TrueVsReco[i]->SetOption("colz");
    h_nelec_Resolution[i] = new TH1D( Form("pl%i_nelec_res",i),Form("Plane %i;Hit charge resolution: (reco-true)/true",i),200,-2,2);
  }
  
  // Energy threshold
  dir_thresh = fOutFile->mkdir("EnergyThreshold");
  dir_thresh->cd();
  float emax = 4000.; // keV
  int   ebins = 40; 
  h_thresh_true_energy = new TH1D("true_energy","True blip energy;True electron energy dep [keV];",ebins,0,emax);
  h_thresh_reco3D = new TH1D("reco3D","Reconstructed 3D blips;True electron energy dep [keV]",ebins,0,emax);
  h_thresh3D = new TH1D("3Dthresh","Reconstructed threshold 3D;True electron energy dep [keV]",ebins,0,emax);
  for(int i=0; i<kNplanes; i++){
    h_thresh_reco[i] = new TH1D(Form("pl%i_reco",i),Form("Reconstructed clusters, plane %i;True electron energy dep [keV]",i),ebins,0,emax);
    h_thresh[i] = new TH1D(Form("pl%i_threshold",i),Form("Reconstruction threshold, plane %i;True electron energy dep [keV]",i),ebins,0,emax);
  }
}



//#################################################################################
void makePlots(){
  
  bool makeThresholdPlot  = true;
  bool makeHitNumPlot     = true;
  bool makeHitMetricPlot  = true;

  if(makeThresholdPlot){
    // Make threshold plot
    dir_thresh->cd();
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
    
    
    gr_thresh[2]->GetXaxis()->SetTitle("Deposited Energy [keV]");
    gr_thresh[2]->GetYaxis()->SetTitle("Reconstruction Efficiency");
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
    
    float thresh_value[3]={};
    for(size_t i=0; i<kNplanes;i++){
      for(size_t e=0; e<5000; e+=10){
        float eval = gr_thresh[i]->Eval(e);
        if( thresh_value[i]==0 && eval >= 0.5 ) {
          thresh_value[i] = e;
          std::cout<<"50% threshold crossing for plane "<<i<<" is "<<e<<" keV.\n";
          break;
        }
      }
    }
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
    h[2]->GetXaxis()->SetRangeUser(0,100);
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

  dir_diag->cd();

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
    ScaleHist(h_hitwidth[i][0],"integral");
    ScaleHist(h_hitwidth[i][1],"integral");
    TH1D* h0 = Clone(h_hitwidth[i][0]);
    TH1D* h1 = Clone(h_hitwidth[i][1]);

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
  TCanvas* c_hwph = new TCanvas("hit_widths_vs_ph","hit_widths_vs_ph",1200,400);
  c_hwph->Divide(3,1);
  for(int i=0; i<kNplanes;i++) {
    c_hwph->cd(i+1);
    gPad->SetMargin(0.15,0.03,0.1,0.1);
    gStyle->SetOptStat(0);
    TH2D* h0 = Clone(h_hitwidth_vs_ph[i][0]);
    TH2D* h1 = Clone(h_hitwidth_vs_ph[i][1]);
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
