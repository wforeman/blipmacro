//////////////////////////////////////////////////////////////////
// 
//  Analysis ROOT Macro
//
////////////////////////////////////////////////////////////////// 

#include "core/vars_bipo.h"
#include "core/tools.h"
#include <algorithm>
#include <time.h>

// Choose dataset to look at
int   fPhase            = 2;

// Input files
std::string   fileDir   = "files";
std::string   outDir    = "output";
std::string   files[2]  = { "BlipAna_RadonData_Phase1_20220217.root",
                            "BlipAna_RadonData_Phase2_20220217.root" };

// Special switches
bool  fDoWireDiagnostics= false; //true;
int   fRandomWireShift  = 0; //1000; 

// General macro parameters
int   fMinTick          = 0;
float fMaxClusterSpan   = 30;     // Veto clusters longer than this [ticks]
int   fBetaMinPlanes    = 2;      // Min number of matched planes (must be 2 or 3)
float fMinHitPH         = -999;   // Don't consider clusters with lead hits
int   fWireRange        = 1;      // +/- range to look for alpha candidate
float fBetaMaxDiff      = 3;      // Difference in wire intersection points [cm]
float fBetaCharge_min   = 0;      // Min charge of beta candidate blip [e-]
float fBetaCharge_max   = 60e3;   // Max charge of beta candidate blip [e-]
int   fBetaHits_max     = 999;
float fAlphaCharge_min  = 0;      // Min charge of alpha candidate cluster [e-]
float fAlphaCharge_max  = 15e3;   // Max charge of alpha candidate cluster [e-]
int   fAlphaHits_max    = 999;      // Max hits in alpha candidate cluster
bool  fAlphaReq3D       = false;  // Require alpha pulse be 3D-matched
float fdT_binSize       = 20.;    // Bin width for all dT spectra plots [us]
float fdT_min           = 20.;    // Min dT for looking for candidate [us]
float fdT_max           = 800.;   // Max dT for looking for candidate [us]
int   fMaxClustMult     = 2;      // Max number of clusters in time window
int   fMaxCandidates    = 1;      // Max number of alpha-like candidates
float fSphereRadius     = 25;     
int   fMinSphereMult2D  = -999; //0;
int   fMaxSphereMult2D  =  999; //10;
int   fMinSphereMult3D  = -999; //0;
int   fMaxSphereMult3D  =  999; //10;
float fZmin             = 50;     // Z range (0 to 1037 cm)
float fZmax             = 985;    //
float fYmin             = -70;    // Y range (-120 to 120 cm)
float fYmax             = 70;     //

// Detector properties
int   nWiresColl        = 3455;
float fSamplePeriod     = 0.5; // microseconds

// Time periods in UNIX time 
//                            Phase 1     Phase 2
//                            ~49.3 hrs   ~46.4 hrs
unsigned int phase_T0[2] = {  1627415210, 1627594369 };
unsigned int phase_T1[2] = {  1627592728, 1627761265 };

// Output options
int printInterval       = 1000;

//##########################################################################3

unsigned int T0             = phase_T0[fPhase-1];
std::string   fFileName     = files[fPhase-1];
std::string   fInFilePath   = fileDir + "/" + fFileName;
std::string   fOutFilePath  = outDir + "/" + "plots_bipo.root";
std::string   fTreeName     = "blipana/anatree";

float dz = fZmax-fZmin;
float dy = fYmax-fYmin;
float fFidCorFactor = (1037.*240.) / (dz*dy);

// Counters
int numEvents       = 0;
int numBiPo         = 0;
int numBiPo_150_300 = 0;
float totalLiveTime = 0;
unsigned int minUnixTime  = 0;
unsigned int maxUnixTime  = 0;
int   maxWirePlane[3]     = { -999, -999, -999 };

// Derived parameters
int   minTick         = fMinTick;
int   maxTick         = 6400 - (int)fdT_max*2;
float liveTimePerEvt  = (maxTick-minTick)*fSamplePeriod*1e-6; //sec

// Vector of noisy wires to skip (collection plane)
std::vector<int> noisyWires{
/*
17, 21, 83, 352, 1247, 1540, 1577, 2335, 2400, 2415, 2464, 2704, 3391, 3408, 3409, 3410, 3411, 3412, 3413, 3414, 3415, 3416, 3417, 3418, 3419, 3420, 3421, 3422, 3423, 3424, 3425, 3426, 3427, 3428, 3429, 3430, 3431, 3432, 3433, 3434, 3435, 3436, 3437, 3438, 3439, 3440, 3441, 3442, 3443, 3444, 3445, 3446, 3447, 3448
*/
};

//##########################################################################3

// Struct to hold fit results
struct FitResult { float rate_signal = -9, rate_bg = -9, ratio = -9; };

// Functions 
void      makePlots();
void      makeHistograms();
void      configure();
FitResult fitdT(TH1D*,bool);

// ROOT objects
TTree*      fTree;
TFile*      fOutFile;
TDirectory* fOutFile_plots;

// Histograms
TH1D* h_blip_charge;
TH1D* h_alpha_charge;
TH1D* h_blip_nhits;
TH1D* h_alpha_nhits;
TH1D* h_clust_timespan;
TH1D* h_nclusts_wire_ave;
TH1D* h_nclusts_inwindow;
TH1D* h_ncands_inwindow;
TH1D* h_cand_precut_dT;
TH1D* h_cand_dT;
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
TH1D* h_clust_prox;
TH1D* h_blip_prox;
TH1D* h_clust_sphere_mult;
TH1D* h_blip_sphere_mult;
TH2D* h_time_vs_dT;
TH1D* h_time_vs_N;
TH1D* h_time_vs_rate_bipo;
TH1D* h_time_vs_rate_bg;
TH1D* h_time_vs_ratio;


//#################################################################################
void makeHistograms()
{
  if( fDoWireDiagnostics )
    h_nclusts_wire_ave  = new TH1D("nclusts_perwire","Average clusters per wire;Average number of clusters per wire",500,0,1); 
  h_clust_timespan    = new TH1D("clust_timespan","Cluster timespan;Ticks",200,0,200);
  h_blip_charge       = new TH1D("blip_charge","3D Blips;Collection Plane Charge [e];Events", 600,0,3e5);
  h_alpha_charge       = new TH1D("alpha_charge","Candidate alphas;Collection Plane Charge [e];Events", 600,0,3e5);
  h_blip_nhits       = new TH1D("blip_nhits","Candidate blips;Collection Plane Hits",20,0,20);
  h_alpha_nhits       = new TH1D("alpha_nhits","Candidate alphas;Collection Plane Hits",20,0,20);
  h_nclusts_inwindow = new TH1D("nclusts_inwindow","Number of clusters in time window following Bi-candidate",20,0,20);
  h_ncands_inwindow = new TH1D("ncands_inwindow","Number of Po candidates in time window following Bi-candidate",10,0,10);
  h_clust_prox  = new TH1D("clust_prox","Cluster separation 'distance' (exclude same+adjacent wire);Cluster separation [cm]",60,0,30);
  h_blip_prox  = new TH1D("blip_prox","3D blip separation;Distance [cm]",60,0,30);
  h_clust_sphere_mult = new TH1D("clust_sphere_mult",Form("Sphere radius %f cm;Cluster multiplicity",fSphereRadius),50,0,50);
  h_blip_sphere_mult = new TH1D("blip_sphere_mult",Form("Sphere radius %f cm;3D blip multiplicity",fSphereRadius),50,0,50);

  float Zmin = -100;  float Zmax = 1100;  int Zbins = 300;
  float Ymin = -150;  float Ymax = 150;   int Ybins = 75;
  float Tmin = 0;     float Tmax = 6400;  int Tbins = 3200/2; 
  float Wmin = -100;  float Wmax = 3500;  int Wbins = 1800; 
  h_zy_blips = new TH2D("zy_blips","3D blips;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
  h_zy_blips_filt = new TH2D("zy_blips_filt","3D blips (quality cuts);Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
  h_zy_bipos = new TH2D("zy_bipos","BiPo candidates;Z [cm]; Y [cm]",Zbins,Zmin,Zmax,Ybins,Ymin,Ymax);
  h_zy_blips->SetOption("colz");
  h_zy_blips_filt->SetOption("colz");
  h_zy_bipos->SetOption("colz");
  h_wt_clusts = new TH2D("wt_clusts","2D clusts;Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
  h_wt_blips = new TH2D("wt_blips","3D blips;Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
  h_wt_blips_filt = new TH2D("wt_blips_filt","3D blips (quality cuts);Collection Plane Wire; Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
  h_wt_bipos = new TH2D("wt_bipos","BiPo candidates;Collection Plane Wire;Ticks",Wbins,Wmin,Wmax,Tbins,Tmin,Tmax);
  h_wt_clusts->SetOption("colz");
  h_wt_blips->SetOption("colz");
  h_wt_blips_filt->SetOption("colz");
  h_wt_bipos->SetOption("colz");
  
  int dTbins = fdT_max / fdT_binSize;
  h_cand_precut_dT= new TH1D("cand_precut_dT","Possible BiPo Candidates (no cuts);Time Difference [#mus];Events",  dTbins,0,fdT_max);
  h_cand_dT       = new TH1D("cand_dT","Selected BiPo Candidates;Time Difference [#mus];Candidates per second", dTbins,0,fdT_max);
  h_clust_dT      = (TH1D*)h_cand_dT->Clone("clust_dT");
  h_clust_dT      ->SetTitle("Same-wire cluster separations");
  
  int timeBins    = 10; float timeMax = 50;
  h_time_vs_dT    = new TH2D("time_vs_dT",";Event time [hr];Time difference [#mus]",timeBins,0,timeMax, dTbins,0,fdT_max);
  h_time_vs_dT    ->SetOption("colz");
  h_time_vs_N     = new TH1D("time_vs_N",";Event time [hr];Number of entries into dT plot",timeBins,0,timeMax);
  h_time_vs_rate_bipo = new TH1D("time_vs_rate_bipo","BiPo component of dT fit;Event time [hr];Rate [sec^{-1}]",timeBins,0,timeMax);
  h_time_vs_rate_bg   = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_rate_BG");
  h_time_vs_rate_bg   ->SetTitle("Background component");
  h_time_vs_ratio     = (TH1D*)h_time_vs_rate_bipo->Clone("time_vs_ratio");
  h_time_vs_ratio     ->SetTitle("Signal-to-background ratio");
  h_time_vs_ratio     ->GetYaxis()->SetTitle("Signal-to-background ratio");

}

//#################################################################################
void configure()
{
  // open the file and set up the TTree
  std::cout<<"Opening file\n";
  TFile* file = new TFile(fInFilePath.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
  
  // make output file to store plots
  fOutFile = new TFile(fOutFilePath.c_str(), "recreate");
  fOutFile_plots  = fOutFile->mkdir("plots");
  
  // initialize all histograms
  makeHistograms();
}


//#################################################################################
void BiPo_macro()
{
  // Configure histograms and TFile
  configure();

  // Keep list of cluster counts per wire/channel
  std::map<int,int> map_wn;
  std::map<int,int> map_chn;
  
  // Record start-time
  std::time_t loopStart = time(0);

  // ========================================================
  // Loop over the events
  int counter = 0;
  int totalEvents = fTree->GetEntries();
  for(int iEvent=0; iEvent < totalEvents; iEvent++){
    
    counter++;
    if( counter > printInterval ) {
      std::cout<<"========== EVENT "<<iEvent<<" / "<<totalEvents<<" ========================\n";
      counter = 1;
    }
    
    // ..... quick-test options ...........
    //int maxEvt    = 1000; if(  iEvent >= maxEvt ) break;
    //int sparsify  = 100; if(  (iEvent % sparsify) != 0 ) continue; 
    //..................................

    // Retrieve event info
    fTree->GetEntry(iEvent);
    numEvents++;
    
    //if( timestamp < minUnixTime || minUnixTime == 0 ) minUnixTime = timestamp;
    //if( timestamp > maxUnixTime || maxUnixTime == 0 ) maxUnixTime = timestamp;
    float eventHour = (timestamp-T0) / 3600.;
    h_time_vs_N->Fill(eventHour);
 
    // -----------------------------------------------------
    // Map of clust IDs per wire on collection plane
    std::vector<int> v_clusts_pl2;
    std::map<int, std::vector<int>> map_wid;
    for(int i=0; i < nclusts; i++){
      int plane = clust_plane[i];
      if( clust_wire[i] > maxWirePlane[plane] ) maxWirePlane[plane] = clust_wire[i];
      if( plane != 2 ) continue;
      v_clusts_pl2.push_back(i);
      map_wid[clust_wire[i]].push_back(i);
      map_wn[clust_wire[i]] += 1;
      h_clust_timespan->Fill(clust_timespan[i]);
      h_wt_clusts->Fill( clust_wire[i], clust_lhit_peakT[i] );
    }
      
    if( fDoWireDiagnostics ) {
      // For each wire, plot cluster separations
      for(int i=0; i<=nWiresColl; i++){
        if( map_wid.find(i) == map_wid.end() ) continue;
        std::vector<float> v_t;
        for(int j=0; j<(int)map_wid[i].size(); j++){
          int id = map_wid[i].at(j);
          v_t.push_back( clust_time[id] * fSamplePeriod);
        }
        std::sort(v_t.begin(), v_t.end());
        for(int j=1; j<(int)v_t.size(); j++){
          h_clust_dT->Fill(v_t.at(j)-v_t.at(j-1));
        }
      }
    }
    

    // Masks to track blip/cluster availability to avoid double-counts
    std::vector<bool> blipAvailable(nblips, true);
    std::vector<bool> clustAvailable(nclusts, true);

    // -----------------------------------------------------
    // Loop over 3D blips...
    for(int i=0; i<nblips; i++){

      // Find the next highest-charge blip that hasn't yet been used
      int iBlip = -9;
      for(int ii=0; ii<nblips; ii++){
        if( !blipAvailable[ii] ) continue;
        if( (iBlip<0) || (blip_charge[ii]>blip_charge[iBlip]) ) iBlip = ii;
      }
      if( iBlip < 0 ) break;
      blipAvailable[iBlip] = false;

      // Find associated cluster on collection
      int   ic    = blip_clustid_pl2[iBlip];

      // Plot ZY position
      h_zy_blips->Fill( blip_z[iBlip], blip_y[iBlip] );
      h_wt_blips->Fill( clust_wire[ic], clust_lhit_peakT[ic] );

      // skip if this cluster was already included in a BiPo candidate 
      if( !clustAvailable[ic] ) continue;
      
      // skip pulse train hits
      if( clust_lhit_gof[ic] < 0  || clust_timespan[ic] > fMaxClusterSpan ) continue;

			// skip if this is a noisy wire
		  if( std::find(noisyWires.begin(), noisyWires.end(), clust_wire[ic] ) != noisyWires.end() ) continue;
    
      // Hit pulse-height cut
      if( clust_lhit_ph[ic] < fMinHitPH ) continue;

      // 3D match cuts for beta
      if( blip_nplanes[iBlip] < fBetaMinPlanes || blip_maxdiff[iBlip] > fBetaMaxDiff ) continue;

      h_zy_blips_filt->Fill( blip_z[iBlip], blip_y[iBlip] );
      h_wt_blips_filt->Fill( clust_wire[ic], clust_lhit_peakT[ic] );
      h_blip_charge->Fill(clust_charge[ic]);
      h_blip_nhits->Fill(clust_nhits[ic]);
      
      // evaluate if in fiducial volume
      bool inFidVol = false;
      if( blip_z[iBlip] > fZmin && blip_z[iBlip] < fZmax &&
          blip_y[iBlip] > fYmin && blip_y[iBlip] < fYmax ) inFidVol=true;
      
      if( !inFidVol ) continue;

      // apply charge/size cuts on beta
      if(   clust_charge[ic] < fBetaCharge_min 
        ||  clust_charge[ic] > fBetaCharge_max 
        ||  clust_nhits[ic] > fBetaHits_max ) continue;
    
      // check for nearby 3D blips
      int nMultSphere3D = 0;
      for(int ii=0; ii<nblips; ii++){
        if( ii == iBlip ) continue;
        TVector3 p1(blip_x[iBlip],blip_y[iBlip],blip_z[iBlip]);
        TVector3 p2(blip_x[ii],blip_y[ii],blip_z[ii]);
        float ds = (p2-p1).Mag();
        h_blip_prox->Fill(ds);
        if(ds<fSphereRadius) nMultSphere3D++; 
      }
      h_blip_sphere_mult->Fill(nMultSphere3D);
    
      // skip if we are near end of wire (account for 400us/800tick trigger offset)
      if( (clust_lhit_peakT[ic]) < minTick ) continue;
      if( (clust_lhit_peakT[ic]) > maxTick ) continue;
     
      // useful structure to save candidate info in
      struct Candidate { int id1, id2; float dT, q1, q2; };

      // -------------------------------------------------------------------- 
      // analyze all other collection plane clusts (2D) in the context of this 
      // Bi-candidate (3D). Also keep track of other clusters in vicinity/sphere.
      std::vector<Candidate> v_cands;
      int nMultSphere2D = 0;
      int nclusts_inwindow = 0;
      for(auto& jc : v_clusts_pl2 ) {
        if( ic == jc ) continue;
        if( !clustAvailable[jc] ) continue;
        if( clust_lhit_gof[jc] < 0 ) continue;
        if( clust_lhit_ph[jc] < fMinHitPH ) continue;
        float dT  = (clust_time[jc]-clust_time[ic])*fSamplePeriod;
        int   ref = clust_wire[ic] + fRandomWireShift;
        if( fRandomWireShift != 0 ) {
          if( ref < 0 ) ref = nWiresColl - abs(ref);
          else if ( ref > nWiresColl ) ref -= nWiresColl;
        }
        float dW  = clust_wire[jc] - ref;
        float q1  = clust_charge[ic];
        float q2  = clust_charge[jc];
        
        // Basic time/wire requirements for BiPo
        bool wireProx = false;
        if( abs(dW) <= fWireRange ) wireProx = true;

        // Check 'distance' in 2D wire-time coordinates, and keep
        // count of clusters that are within sphere but not on the 
        // same or adjacent wires as the Bi-candidate.
        float dx_wt = dW * 0.3;     // 0.3cm wire spacings
        float dy_wt = dT * 0.1041;  // ~0.1041 cm/us drift speed
        float ds_wt = sqrt( pow(dx_wt,2) + pow(dy_wt,2) );
        if( !wireProx ) {
          h_clust_prox->Fill( ds_wt, 1./pow(ds_wt,3) );//weight by 1/r^3
          if( ds_wt < fSphereRadius ) nMultSphere2D++;
        }
       

        // If we are requiring the alha be plane-matched, it must
        // have an associated 3D blip.
        if( fAlphaReq3D && clust_blipid[jc] < 0 ) continue;

        // --- time window cuts ---
        if( dT > 0 && dT < fdT_max && wireProx ) {
          nclusts_inwindow++;
          h_cand_precut_dT->Fill(dT);
          h_alpha_charge->Fill(q2);
          h_alpha_nhits->Fill(clust_nhits[jc]);
          
          // --- alpha charge/nhits cut ---
          if(   q2 > fAlphaCharge_min 
            &&  q2 < fAlphaCharge_max 
            &&  clust_nhits[jc] <= fAlphaHits_max ) {
          
            // Save this as a possible candidate alpha
            Candidate c = { ic, jc, dT, q1, q2 };
            v_cands.push_back(c);
         
          }//alpha-like cut
        
        }//end cut on wire proximity and dT > 0

      }//endloop over secondary clusters
      
      h_clust_sphere_mult->Fill(nMultSphere2D);
      
      // -- skip this candidate if there were a high
      //    number of clusters within proximity
      if( fMinSphereMult2D >= 0 && nMultSphere2D < fMinSphereMult2D ) continue;
      if( fMaxSphereMult2D >= 0 && nMultSphere2D > fMaxSphereMult2D ) continue;
      
      h_nclusts_inwindow->Fill(nclusts_inwindow);
      h_ncands_inwindow->Fill((int)v_cands.size());

      // Make a final candidate selection stored as "finalCand"
      Candidate finalCand;
      bool candFound = false;
     
      // ..if candidates were found (but not too many)...
      if( v_cands.size() >0 && v_cands.size() <=  fMaxCandidates
          && nclusts_inwindow <=  fMaxClustMult 
        ) {
        
        // if just one, take it!
        if( v_cands.size() == 1 ) {
          finalCand = v_cands.at(0);
          candFound = true;
        }

        // otherwise, choose the one with 
        // the smallest dT
        else {
          float min = 9999;
          for(size_t k=0; k<v_cands.size(); k++){
            if( v_cands.at(k).dT < min ) {
              finalCand = v_cands.at(k);
              candFound = true;
            }
          }
        }
      }//end final cand selection

      
      // If a final selection was made, mark the used
      // clusters and blips as unavailable and make a
      // final fiducial cut.
      if( candFound ) {
        
        clustAvailable[finalCand.id1] = false;
        clustAvailable[finalCand.id2] = false;
        int blipid = clust_blipid[finalCand.id2];
        if( blipid >= 0 ) blipAvailable[blipid] = false;
    
        // Require min dT
        if( finalCand.dT < fdT_min ) continue;
        
        // Plot candidate locations
        h_zy_bipos->Fill( blip_z[iBlip], blip_y[iBlip] );
        h_wt_bipos->Fill( clust_wire[finalCand.id1], clust_lhit_peakT[finalCand.id1] );
        
        // Final fiducial cut
        if( !inFidVol ) continue;
       
        numBiPo++;
        if( finalCand.dT >= 150 && finalCand.dT <= 300 ) numBiPo_150_300++;
        h_cand_dT->Fill(finalCand.dT);
        h_time_vs_dT->Fill(eventHour,finalCand.dT);
      }


    }//end loop over 3D blips      

  }//endloop over events

  double loopDuration = ( time(NULL) - loopStart );
  
  // Calculate total live time in seconds (sampling period 0.5us/tick)
  totalLiveTime = float(numEvents) * liveTimePerEvt;
  
  // Scale dT plots so they're 'per second'. 
  float scaleFact = 1./ ( totalLiveTime );
  h_cand_precut_dT->Scale( scaleFact );
  h_cand_dT->Scale( scaleFact );
  h_clust_dT->Scale( scaleFact );

  // Check for noisy wires
  if( fDoWireDiagnostics ) {
    std::cout<<"Checking noisy wires...\n";
    int nNoisy = 0;
    float noiseThresh = 0.08;
    // Populate histogram of average cluster multiplicity per wire
    if( map_wn.size() ) {
      for(int iwire=0; iwire<=nWiresColl; iwire++){
        if( map_wn.find(iwire) != map_wn.end() ) {
          float a = map_wn.at(iwire) / (float)numEvents;
          h_nclusts_wire_ave->Fill(a);
          if( a > noiseThresh ) {
            std::cout<<iwire<<", ";
            nNoisy++;
          } 
        } else { 
          h_nclusts_wire_ave->Fill(0);
        }
      }
    }
    std::cout<<"\n--> Found "<<nNoisy<<" wires with ave clusters/evt > "<<noiseThresh<<"\n";
  }

 
  // Write all histos currently in stack
  fOutFile->Write(); 

  // Make plots
  makePlots();
  
  printf("\n*******************************************\n");
  printf("File                : %s\n",      fFileName.c_str()); 
  printf("Total events        : %i\n",      numEvents); 
  printf("Total live time     : %f sec\n",  totalLiveTime); 
  printf("Ave cands / evt     : %f\n",      h_cand_dT->GetEntries()/(float)numEvents);
  printf("Processing time     : %f sec/evt\n", loopDuration/float(numEvents));
  printf("Excluding %i noisy wires \n",     (int)noisyWires.size());	
  printf("\n*******************************************\n\n");
  
  //std::cout<<"min: "<<minUnixTime<<"\n";
  //std::cout<<"max: "<<maxUnixTime<<"\n";
  //std::cout<<"diff: "<<(maxUnixTime - minUnixTime)/3600.<<" hours\n";
  //std::cout<<"max wire plane 0 "<<maxWirePlane[0]<<"\n";
  //std::cout<<"max wire plane 1 "<<maxWirePlane[1]<<"\n";
  //std::cout<<"max wire plane 2 "<<maxWirePlane[2]<<"\n";

  // Close file
  fOutFile->Close();

}


//#################################################################################
void makePlots()
{
  
  // ============================================
  // Do final fit on dT spectrum
  // ============================================
  fitdT( h_cand_dT, true );
  

  // ============================================
  // Rate vs time plots
  // ============================================
  std::cout<<"\nMaking rate vs time plot...\n";
  TH1D* h_slice;
  for(int i=1; i<=(int)h_time_vs_N->GetXaxis()->GetNbins(); i++){
    h_slice = Make1DSlice( h_time_vs_dT, i, i, Form("dTfit_%i",i) );
    float N = h_time_vs_N->GetBinContent(i);
    h_slice->Scale( 1. / (N*liveTimePerEvt) );
    std::cout<<"time: "<<h_time_vs_N->GetXaxis()->GetBinCenter(i)<<" hrs\n";
    FitResult fr = fitdT(h_slice,false);
    if( fr.rate_signal > 0 ) {
      h_time_vs_rate_bipo ->SetBinContent(i,fr.rate_signal);
      h_time_vs_rate_bg   ->SetBinContent(i,fr.rate_bg);
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
  
  std::string name = "c_time_vs_rate";
  TCanvas* c = new TCanvas(name.c_str(),name.c_str(),600,500);
  float max = std::max( GetHistMax(h1), GetHistMax(h2) );
  h1->GetYaxis()->SetRangeUser(0,max*1.2);
  h1->DrawCopy();
  h2->DrawCopy("same");
  fOutFile_plots->cd();
  c->Write();
  
} 

//################################################################################
// Function that performs the dT fit
//
FitResult fitdT(TH1D* h, bool writeCanvas = false ){

  FitResult out;

  if( h->GetEntries() < 20 ) return out;

  std::cout<<"Fitting dT spectrum "<<h->GetTitle()<<", "<<h->GetEntries()<<"\n";
  std::string label = h->GetName();
  TCanvas* c = new TCanvas(Form("c_fit_%s",label.c_str()),Form("c_fit_%s",label.c_str()),600,500);
  TH1D* hc = (TH1D*)h->Clone();
  
  // Define fit function, initialize parameters
  float histMax = GetHistMax(hc);
  TF1* fit = new TF1("FullFit","[0] + [1]*exp(-x/[2]) + [3]*exp(-x/[4])");
  fit->SetParameter(0, histMax/20 );
  fit->SetParLimits(0, 0, histMax );
  fit->SetParameter(1, histMax/5 );
  fit->SetParLimits(1, 0, histMax );
  fit->SetParameter(2, 15 );
  fit->SetParLimits(2, 0, 50);
  fit->SetParameter(3, histMax/10 );
  fit->SetParLimits(3, 0, histMax );
  fit->FixParameter(4, 164.5 );
  
  // TEMP: fix flat background
  //fit->FixParameter(0, 0.001597);

  // Draw plot and fit
  c->cd();
  hc->Fit(fit,"Q"); 
  hc->DrawCopy();
  fit->Draw("same");
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

  TF1* f_expBG = (TF1*)f_exp->Clone("expBG");
  f_expBG->FixParameter(0, fit->GetParameter(1) );
  f_expBG->FixParameter(1, fit->GetParameter(2) );
  
  TF1* f_bipo = (TF1*)f_exp->Clone("bipo");
  f_bipo->FixParameter(0, fit->GetParameter(3) );
  f_bipo->FixParameter(1, fit->GetParameter(4) );
  
  // Full extrapolated BiPo component; integral of A*exp(-x/B) from 0-infinity is A*B
  float n_bipo  = (f_bipo->GetParameter(0)*f_bipo->GetParameter(1))/fdT_binSize;
  
  // Components within the dT selection window
  float N_total = hc      ->Integral(0,hc->GetXaxis()->GetNbins());
  float N_fit   = fit     ->Integral(fdT_min,fdT_max)/fdT_binSize;
  float N_flat  = f_flat  ->Integral(fdT_min,fdT_max)/fdT_binSize;
  float N_expbg = f_expBG ->Integral(fdT_min,fdT_max)/fdT_binSize;
  float N_bipo  = f_bipo  ->Integral(fdT_min,fdT_max)/fdT_binSize;
  float N_bg    = N_flat + N_expbg;
  float sbratio = N_bipo / N_bg;
  
  printf("------------- dT fit --------------------\n");
  printf("dT min/max          : %i-%i us\n",  int(fdT_min),int(fdT_max));
  printf("Chi2/ndf            : %f\n",        fit->GetChisquare()/fit->GetNDF());
  printf("Total entries       : %f\n",        hc->GetEntries());
  printf("Total rate          : %f per sec\n",N_total);
  printf(" - BiPo component   : %f per sec\n",N_bipo);
  printf(" - ExpBG component  : %f per sec\n",N_expbg);
  printf(" - Flat component   : %f per sec\n",N_flat);
  printf("S/BG ratio          : %f \n",sbratio);
  printf("\n");
  printf("Fiducial-corrected  : %f BiPos per sec\n",n_bipo*fFidCorFactor);
  printf("-----------------------------------------\n");

  out.rate_signal = n_bipo;
  out.rate_bg     = N_bg;
  out.ratio       = sbratio;
  return out;

}


