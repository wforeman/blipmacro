//////////////////////////////////////////////////////////////////
// 
//  Analysis ROOT Macro
//
////////////////////////////////////////////////////////////////// 

#include "core/vars.h"
#include "core/tools.h"


// General macro parameters
std::string   fFileName     = "BlipAna_20220627_Eminus_0to5MeV_Overlay_Standard.root";
std::string   fTreeName     = "blipana/anatree";
float         thresh_emax   = 3.0; // MeV
int           thresh_ebins  = 30; // # bins
int           clust_max     = 1000;
int           clust_bins    = 1000;
int           frac_bins     = 50;

// Functions
void    makeHistograms();
void    makePlots();
float   GetThreshold(TGraph*,float);

// ROOT objects
TFile*      fInputFile;
TFile*      fOutFile;
TDirectory* dir_thresh;
TH1D*       h_thresh_true_energy;
TH1D*       h_thresh_reco_2D[kNplanes];
TH1D*       h_thresh_reco_3D;
TH1D*       h_thresh_reco_3D_3plane;
TDirectory* dir_noise;
TH1D*       h_nclusts_total[kNplanes];
TH1D*       h_nclusts_noise[kNplanes];
TH1D*       h_nclusts_unmatched[kNplanes];
TH1D*       h_matchfrac[kNplanes];
TH1D*       h_blip_goodfrac;
TH1D*       h_nmatches[kNplanes];

//#################################################################################
void makeHistograms(){
  
  // Energy threshold
  dir_thresh = fOutFile->mkdir("EnergyThreshold");
  dir_thresh->cd();
  h_thresh_true_energy    = new TH1D("thresh_true_energy"   ,"True blips            ;Energy Deposited by Electron [MeV]",                          thresh_ebins,0,thresh_emax);
  h_thresh_reco_3D        = new TH1D("thresh_reco_3D"       ,"3D blips (2-3 planes) ;Energy Deposited by Electron [MeV];Reconstruction Efficiency",thresh_ebins,0,thresh_emax);
  h_thresh_reco_3D_3plane = new TH1D("thresh_reco_3D_3plane","3D blips (3 planes)   ;Energy Deposited by Electron [MeV];Reconstruction Efficiency",thresh_ebins,0,thresh_emax);
  for(int i=0; i<kNplanes; i++){
    h_thresh_reco_2D[i]   = (TH1D*)h_thresh_true_energy->Clone(Form("thresh_reco_2D_pl%i",i));
    h_thresh_reco_2D[i] ->SetTitle(Form("2D clusters, plane %i;Energy Deposited by Electron;Reconstruction Efficiency [MeV]",i));
  }

  // Noise metrics
  dir_noise = fOutFile->mkdir("NoiseMetrics");
  dir_noise->cd();
  for(int i=0; i<kNplanes; i++){
    h_nclusts_total[i]      = new TH1D(Form("pl%i_nclusts",i),          Form("2D clusters, plane %i;Number of clusters",i),                   clust_max,0,clust_bins);   
    h_nclusts_noise[i]      = new TH1D(Form("pl%i_nclusts_noise",i),    Form("Non-truth-matched 2D clusters, plane %i;Number of clusters",i), clust_max,0,clust_bins);
    h_nclusts_unmatched[i]  = new TH1D(Form("pl%i_nclusts_unmatched",i),Form("Non-plane-matched 2D clusters, plane %i;Number of clusters",i), clust_max,0,clust_bins); 
    h_matchfrac[i]          = new TH1D(Form("pl%i_matchfrac",i),        Form("Fraction of clusters plane-matched, plane %i",i),               frac_bins,0,1.+1./frac_bins);
    h_nmatches[i]           = (TH1D*)fInputFile->Get(Form("blipana/BlipRecoAlg/pl%i_nmatches",i));
  }
  h_blip_goodfrac      = new TH1D("blip_goodfrac","Fraction of perfectly truth-matched blips",101,0,1.01);

}


//#################################################################################
void BlipAna_macro(){
 
  //===================================================
  // open the file and set up the TTree
  //===================================================
  fInputFile = new TFile(("files/"+fFileName).c_str(),"READ");
  TTree* fTree = (TTree*)fInputFile->Get(fTreeName.c_str());
  fTree->SetBranchAddress("event",&event);                         
  fTree->SetBranchAddress("nedeps",&nedeps);                       
  fTree->SetBranchAddress("edep_g4id",&edep_g4id);                  
  fTree->SetBranchAddress("edep_blipid",&edep_blipid);              
  fTree->SetBranchAddress("edep_clustid",&edep_clustid);            
  fTree->SetBranchAddress("edep_energy",&edep_energy);              
  fTree->SetBranchAddress("edep_charge",&edep_charge);              
  fTree->SetBranchAddress("edep_x",&edep_x);                        
  fTree->SetBranchAddress("edep_y",&edep_y);                        
  fTree->SetBranchAddress("edep_z",&edep_z);                        
  fTree->SetBranchAddress("nclusts",&nclusts);                     
  fTree->SetBranchAddress("clust_ismatch",&clust_ismatch);
  fTree->SetBranchAddress("clust_plane",&clust_plane);              
  fTree->SetBranchAddress("clust_wire",&clust_wire);                
  fTree->SetBranchAddress("clust_startwire",&clust_startwire);                
  fTree->SetBranchAddress("clust_endwire",&clust_endwire);                
  fTree->SetBranchAddress("clust_charge",&clust_charge);            
  fTree->SetBranchAddress("clust_time",&clust_time);                
  fTree->SetBranchAddress("clust_edepid",&clust_edepid);            
  fTree->SetBranchAddress("clust_blipid",&clust_blipid);            
  fTree->SetBranchAddress("nblips",&nblips);                       
  fTree->SetBranchAddress("blip_x",&blip_x);                        
  fTree->SetBranchAddress("blip_y",&blip_y);                        
  fTree->SetBranchAddress("blip_z",&blip_z);                       
  fTree->SetBranchAddress("blip_energy",&blip_energy);
  fTree->SetBranchAddress("blip_charge",&blip_charge);              
  fTree->SetBranchAddress("blip_edepid",&blip_edepid);              
  fTree->SetBranchAddress("blip_nplanes",&blip_nplanes);
  fTree->SetBranchAddress("blip_sigmayz",&blip_sigmayz);
  fTree->SetBranchAddress("blip_pl0_clustid",&blip_clustid[0]);
  fTree->SetBranchAddress("blip_pl1_clustid",&blip_clustid[1]);
  fTree->SetBranchAddress("blip_pl2_clustid",&blip_clustid[2]);
  
  // make output file to store plots
  fOutFile = new TFile(("output/plots_"+fFileName).c_str(),"recreate");
  makeHistograms();
  
  // Maps, counters, etc
  std::map<int, std::set<int>> map_plane_livewires;
  std::vector<int> maxWire(kNplanes, 0);
  int nblips_truthmatched         = 0;
  int nblips_truthmatched_perfect = 0;
  int print_counter               = 0;
  
  //===================================================
  // Loop over the events
  //===================================================
  int totalEntries = fTree->GetEntries();
  for(int iEvent=0; iEvent<totalEntries; iEvent++){
    fTree->GetEntry(iEvent);
    
    print_counter++;
    if( print_counter > 1000 ) {
      printf("========== EVENT %i / %i =====================\n",iEvent,totalEntries);
      print_counter = 1;
    }
   
    // ------------------------------------------
    // Look at true blips ("edeps")
    // ------------------------------------------
    for(int i=0; i<nedeps; i++) {
      h_thresh_true_energy  ->Fill(edep_energy[i]);
    }

    // ------------------------------------------
    // Loop over the reconstructed clusters
    // ------------------------------------------
    std::vector<int> numclusts(kNplanes,0);
    std::vector<int> numclusts_noise(kNplanes,0);
    std::vector<int> numclusts_matched(kNplanes,0);
    std::vector<int> numclusts_unmatched(kNplanes,0);
    
    for(int i=0; i<nclusts; i++){
      
      int plane = clust_plane[i];
      int blipid = clust_blipid[i];
      int edepid = clust_edepid[i];
      
      numclusts[plane]++;
      if( clust_ismatch[i] ) numclusts_matched[plane]++;
      else                   numclusts_unmatched[plane]++;

      if( clust_endwire[i] > maxWire[plane] ) maxWire[plane] = clust_endwire[i];
      for(int iwire = clust_startwire[i]; iwire <= clust_endwire[i]; iwire++)
      map_plane_livewires[plane].insert(iwire);
      
      if( edepid >= 0 ) h_thresh_reco_2D[plane]  -> Fill(edep_energy[clust_edepid[i]]);
      else              numclusts_noise[plane]++;
      
    }//endloop over clusts

    // ------------------------------------------
    // Loop over 3D blips
    // ------------------------------------------
    for(int i=0; i<nblips; i++){

      int nplanes = blip_nplanes[i];
      int eid     = blip_edepid[i];
      
      if( eid >= 0 ){

        nblips_truthmatched++;

        // only count the blip as "good" if there are no falsely-matched
        // clusters within it, and all clusters match to the same true particle
        bool isGood = true;
        for(size_t ipl = 0; ipl < kNplanes; ipl++){
          if( blip_clustid[ipl][i] <= 0 ) continue;
          int clust_eid = clust_edepid[blip_clustid[ipl][i]];
          if( clust_eid != eid ) {
            isGood = false;
            break;
          }
        }
        
        if( isGood ) {
          nblips_truthmatched_perfect++;
          float true_energy = edep_energy[eid];
          h_thresh_reco_3D-> Fill(true_energy);
          if( nplanes == 3 ) h_thresh_reco_3D_3plane->Fill(true_energy);
        }
      
      }
      
    }//endloop over blips
    
  // Fill histograms
  for(int ipl=0; ipl<kNplanes; ipl++ ) {
    if( !numclusts[ipl] ) continue;
    h_nclusts_total[ipl]    ->Fill(numclusts[ipl]);
    h_nclusts_noise[ipl]    ->Fill(numclusts_noise[ipl]);
    h_nclusts_unmatched[ipl]->Fill(numclusts_unmatched[ipl]);
    h_matchfrac[ipl]        ->Fill(numclusts_matched[ipl]/float(numclusts[ipl]));
  }
  
  if( nblips_truthmatched )
  h_blip_goodfrac->Fill( nblips_truthmatched_perfect / float(nblips_truthmatched) );


  }//endloop over events


  // ========================================
  // Normalize efficiency plots
  // ========================================
  for(int ipl=0; ipl<kNplanes; ipl++)
  h_thresh_reco_2D[ipl]   ->Divide( h_thresh_true_energy );
  h_thresh_reco_3D        ->Divide( h_thresh_true_energy );
  h_thresh_reco_3D_3plane ->Divide( h_thresh_true_energy );
  for(int ipl=0; ipl<kNplanes; ipl++)
  h_matchfrac[ipl]        ->Scale( 1./ totalEntries );

  std::cout<<"\nLive wires per plane:\n";
  std::cout<<map_plane_livewires[0].size()<<" (max wire: "<<maxWire[0]<<")\n";
  std::cout<<map_plane_livewires[1].size()<<" (max wire: "<<maxWire[1]<<")\n";
  std::cout<<map_plane_livewires[2].size()<<" (max wire: "<<maxWire[2]<<")\n";

  fOutFile->Write(); 
  makePlots();
  fOutFile->Close();

}






//#################################################################################
void makePlots(){
  
  TCanvas* c_comb     = new TCanvas("combined","combined",1000,500);
  c_comb->Divide(2,1);

  std::vector< Color_t > colors { kRed, kGreen+2, kBlue, kBlack, kMagenta };

  // =======================================================
  // Threshold/efficiency plot
  // =======================================================
  dir_thresh->cd();
  TGraphErrors* gr_thresh_2D[3];
  TGraphErrors* gr_thresh_3D;
  TGraphErrors* gr_thresh_3D_3plane;
  TLegend* leg_thresh;
  
  for(int i=0; i<kNplanes;i++) 
  gr_thresh_2D[i]     = new TGraphErrors();
  gr_thresh_3D        = new TGraphErrors();
  gr_thresh_3D_3plane = new TGraphErrors();
  
  // Loop over true energy bins
  size_t Nbins = h_thresh_true_energy->GetXaxis()->GetNbins();
  for(size_t j=0; j<Nbins; j++){
    float energy = h_thresh_true_energy->GetBinCenter(j);
    // 2D clusters
    for(int i=0; i<kNplanes;i++)
    gr_thresh_2D[i] ->SetPoint(gr_thresh_2D[i]->GetN(),energy,h_thresh_reco_2D[i]->GetBinContent(j));
    // 3D blips
    gr_thresh_3D        ->SetPoint(gr_thresh_3D       ->GetN(), energy,h_thresh_reco_3D       ->GetBinContent(j));
    gr_thresh_3D_3plane ->SetPoint(gr_thresh_3D_3plane->GetN(), energy,h_thresh_reco_3D_3plane->GetBinContent(j));
  }

  gr_thresh_2D[2]->GetXaxis()->SetTitle("Energy deposited by electron [MeV]");
  gr_thresh_2D[2]->GetYaxis()->SetTitle("Reconstruction efficiency");
  gr_thresh_2D[2]->GetXaxis()->SetRangeUser(0,h_thresh_true_energy->GetXaxis()->GetXmax());
  gr_thresh_2D[2]->GetXaxis()->SetTitleOffset(1.2);
  FormatTGraph(gr_thresh_2D[2],     kBlue,kBlue,20,1,0.7,1);
  FormatTGraph(gr_thresh_2D[1],     kGreen+2, kGreen+2,   22, 1,  0.7,  1);
  FormatTGraph(gr_thresh_2D[0],     kRed,     kRed,       21, 1,  0.7,  1);
  FormatTGraph(gr_thresh_3D,        kBlack,   kBlack,     4,  2,  0.7,  2);
  FormatTGraph(gr_thresh_3D_3plane, kMagenta, kMagenta,   4,  2,  0.7,  2);

  //TCanvas* c_th     = new TCanvas("threshold","threshold",500,500);
  c_comb->cd(1);
  gPad->SetMargin(0.15, 0.05,0.10, 0.10);
  gr_thresh_2D[2]     ->Draw("APL");
  gr_thresh_2D[1]     ->Draw("PL same");
  gr_thresh_2D[0]     ->Draw("PL same");
  gr_thresh_3D        ->Draw("PL same");
  gr_thresh_3D_3plane ->Draw("PL same");
  leg_thresh  = MakeLegend(0.5,0.35,0.035,5,0.25);
  leg_thresh  ->AddEntry(gr_thresh_2D[0],     "Induction U Plane","PL");
  leg_thresh  ->AddEntry(gr_thresh_2D[1],     "Induction V Plane","PL");
  leg_thresh  ->AddEntry(gr_thresh_2D[2],     "Collection Plane","PL");
  leg_thresh  ->AddEntry(gr_thresh_3D,        "3D (2-3 planes)","PL");
  leg_thresh  ->AddEntry(gr_thresh_3D_3plane, "3D (3 planes)","PL");
  leg_thresh  ->Draw();
  //c_th        ->Write();


  // =======================================================
  // Cluster matching metrics (proxy for noise levels)
  // =======================================================
  dir_noise->cd();
  TGraph*   g_mf[kNplanes];
  TLegend*  leg_matchfrac;
  
  c_comb  ->cd(2);
  gPad    ->SetMargin(0.15, 0.05,0.10, 0.10);
  gStyle  ->SetOptStat(0);

  for(int i=0; i<kNplanes; i++){
    g_mf[i] = HistToGraph(h_matchfrac[i]); 
    g_mf[i] ->SetTitle("");
    g_mf[i] ->SetLineColor(colors[i]);
    g_mf[i] ->SetLineWidth(2);
    g_mf[i] ->SetLineStyle(1);
  }
  
  float max = 1.2*std::max( GetHistMax(h_matchfrac[0]), std::max( GetHistMax(h_matchfrac[1]),GetHistMax(h_matchfrac[2])));
  g_mf[0] ->GetXaxis()->SetTitle("Fraction of plane-matched clusters");
  g_mf[0] ->GetYaxis()->SetTitle("Number of events [normalized]");
  g_mf[0] ->GetXaxis()->SetTitleOffset(1.2);
  g_mf[0] ->GetYaxis()->SetRangeUser(0,max);
  g_mf[0] ->Draw("APL");
  g_mf[1] ->Draw("PL same");
  g_mf[2] ->Draw("PL same");
  leg_matchfrac  = MakeLegend(0.45,0.88,0.035,3,0.25);
  leg_matchfrac  ->AddEntry(g_mf[0], Form("Induction U Plane: #sigma = %.2f",h_matchfrac[0]->GetMean()),"L");
  leg_matchfrac  ->AddEntry(g_mf[1], Form("Induction V Plane: #sigma = %.2f",h_matchfrac[1]->GetMean()),"L");
  leg_matchfrac  ->AddEntry(g_mf[2], Form("Collection Plane: #sigma = %.2f",h_matchfrac[2]->GetMean()),"L");
  leg_matchfrac  ->Draw();
  
  
  fOutFile->cd();  
  c_comb->Write();
  

  // ================================================================
  // Summary
  // ===============================================================
  // Live wire fraction per plane:
  //  plane 0: 0.831
  //  plane 1: 0.963
  //  plane 2: 0.900
  //  P( 2 + (1||0) ) = 0.894
  //  P( 0 + 1 + 2 )  = 0.720
  
  std::cout<<"\n Input file: "<<fFileName<<"\n";

  float effthresh = 0.50;
  std::cout<<"\n 50% threshold crossing points\n";
  for(int i=0; i<kNplanes;i++)
  printf(" plane %i       : %4.0f keV\n",  i,         GetThreshold(gr_thresh_2D[i],     effthresh));
  printf(" 3D (2-3 plns) : %4.0f keV (%4.0f keV*)\n",GetThreshold(gr_thresh_3D,        effthresh), GetThreshold(gr_thresh_3D,effthresh*0.894));
  printf(" 3D (3 plns)   : %4.0f keV (%4.0f keV*)\n",GetThreshold(gr_thresh_3D_3plane, effthresh), GetThreshold(gr_thresh_3D_3plane,effthresh*0.720));
  printf(" * = scaled based on max achievable\n");
  
  std::cout<<"\n Plane-match fraction:\n";
  for(int i=0; i<kNplanes;i++)
  printf(" plane %i  : %.3f\n",  i, h_matchfrac[i]->GetMean());

  std::cout<<"\n Match candidate multiplicity per cluster:\n";
  for(int i=0; i<kNplanes && i != 2;i++)
  printf(" plane %i  : %f\n", i, h_nmatches[i]->GetMean());
 
  std::cout<<"\n Blip quality rate: "<<h_blip_goodfrac->GetMean()<<"\n\n";


}

//#######################################################
float GetThreshold(TGraph* g, float thresh){
  for(size_t e=0; e<5000; e+=1){
    float eval = g->Eval(e*1e-3);
    if( eval >= thresh ) return e;
  }

  return -1;

}
