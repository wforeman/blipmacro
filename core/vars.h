//##################################
// Custom typedefs
//##################################
typedef std::vector<string> vstring;
typedef std::vector<short> vshort;

void configure();
void makeHistograms();
void makePlots();

TTree* fTree;
TFile* fOutFile;

// Set global constants and max array sizes
const int kMaxHits  = 100000;
const int kMaxBlips = 5000;
const int kMaxTrks  = 1000;
const int kMaxG4    = 10000;
const int kMaxEDeps = 10000;
const int kNplanes  = 3;  
    
// --- Event information ---   
int     event;                    // event number
int     run;                      // run number
double  timestamp;                // UNIX time

// --- G4 information ---
float total_depEnergy;          // total deposited energy in AV
float total_numElectrons;       // total electrons reaching anode wires
int   nparticles;               // number of G4 particles
int   isPrimary[kMaxG4];        // is primary particle
int   trackID[kMaxG4];          // G4 track ID
int   pdg[kMaxG4];              // PDG
int   nDaughters[kMaxG4];       // number of daughters
int   mother[kMaxG4];           // mother particle
float E[kMaxG4];                // initial energy (MeV)
float endE[kMaxG4];             // final energy (MeV)
float mass[kMaxG4];             // mass (MeV)
float P[kMaxG4];                // momentum (MeV)
float Px[kMaxG4];               // momentum x (MeV)
float Py[kMaxG4];               // momentum y (MeV)
float Pz[kMaxG4];               // momentum z (MeV)
float startPointx[kMaxG4];      // starting x (cm)
float startPointy[kMaxG4];      // starting y (cm)
float startPointz[kMaxG4];      // starting y (cm)
float endPointx[kMaxG4];        // ending x (cm)
float endPointy[kMaxG4];        // ending y (cm)
float endPointz[kMaxG4];        // ending y (cm)
float startT[kMaxG4];           // starting time (us)
float endT[kMaxG4];             // ending time (us)
float pathlen[kMaxG4];          // path length (cm)
int   numElectrons[kMaxG4];     // electrons reaching anode wires
float depEnergy[kMaxG4];        // energy deposited in AV (MeV)
std::vector<std::string>* process=0;// process name

// --- True energy deposit info (derived) ---
int   nedeps;                   // number of true localized energy depositions
int   edep_g4id[kMaxEDeps];     // leading G4 track ID
int   edep_clustid[kMaxEDeps];  // hitclust ID
int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
float edep_charge[kMaxEDeps];   // total electrons reaching anode wires
float edep_x[kMaxEDeps];        // x (cm)
float edep_y[kMaxEDeps];        // y (cm)
float edep_z[kMaxEDeps];        // z (cm)
float edep_ds[kMaxEDeps];       // extent (cm)

// --- Hit information ---
int	  nhits;                    // number of hits
int	  hit_tpc[kMaxHits];        // tpc number
int	  hit_plane[kMaxHits];      // plane number
int	  hit_wire[kMaxHits];       // wire number
int	  hit_channel[kMaxHits];    // channel ID
float	hit_peakT[kMaxHits];      // raw peak time (tick)
float	hit_time[kMaxHits];       // corrected peak time (tick)
float hit_rms[kMaxHits];        // shape RMS
float	hit_ph[kMaxHits];         // amplitude
float	hit_area[kMaxHits];       // charge (area) in ADC units
float hit_charge[kMaxHits];     // reconstructed number of electrons
int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
int   hit_isreal[kMaxHits];     // is this hit real?
int	  hit_g4id[kMaxHits];       // G4 TrackID of leading particle
float hit_g4frac[kMaxHits];     // fraction of hit energy from leading MCParticle
float hit_g4energy[kMaxHits];   // true energy
float hit_g4charge[kMaxHits];   // true number of electrons (drift-attenuated)
int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
float hit_gof[kMaxHits];        // hit goodness-of-fit

// --- Hit cluster information ---
int   nclusts;
int   clust_tpc[kMaxHits];
int   clust_plane[kMaxHits];
int   clust_wire[kMaxHits];
int   clust_nwires[kMaxHits];
int   clust_nhits[kMaxHits];
float clust_charge[kMaxHits];
float clust_time[kMaxHits];
float clust_startTime[kMaxHits];
float clust_endTime[kMaxHits];
float clust_xpos[kMaxHits];
float clust_g4energy[kMaxHits];
float clust_g4charge[kMaxHits];
int   clust_g4id[kMaxHits];
int   clust_ismatched[kMaxHits];
int   clust_blipid[kMaxHits];
int   clust_edepid[kMaxHits];

// --- Blip information ---
int   nblips;
int   blip_tpc[kMaxBlips];
int   blip_bestplane[kMaxBlips];
int   blip_ncrossings[kMaxBlips];
float blip_x[kMaxBlips];
float blip_y[kMaxBlips];
float blip_z[kMaxBlips];
float blip_prms[kMaxBlips];
float blip_charge[kMaxBlips];
int   blip_edepid[kMaxBlips];


//################################
// AnaTree branches
//################################

void setBranches(TTree *tree){
  
  // for use on saved vectors of strings, for example:
  // tree->SetBranchAddress("processname",&_processname,&br);
  TBranch *br = 0; 
  tree->SetBranchAddress("event",&event);                         
  tree->SetBranchAddress("run",&run);                             
  tree->SetBranchAddress("total_depEnergy",&total_depEnergy);     
  tree->SetBranchAddress("total_numElectrons",&total_numElectrons);
  tree->SetBranchAddress("nparticles",&nparticles);               
  tree->SetBranchAddress("isPrimary",&isPrimary);                  
  tree->SetBranchAddress("trackID",&trackID);                      
  tree->SetBranchAddress("pdg",&pdg);                              
  tree->SetBranchAddress("nDaughters",&nDaughters);                
  tree->SetBranchAddress("mother",&mother);                        
  tree->SetBranchAddress("E",&E);                                  
  tree->SetBranchAddress("endE",&endE);                            
  tree->SetBranchAddress("mass",&mass);                            
  tree->SetBranchAddress("P",&P);                                  
  tree->SetBranchAddress("Px",&Px);                                
  tree->SetBranchAddress("Py",&Py);                                
  tree->SetBranchAddress("Pz",&Pz);                                
  tree->SetBranchAddress("startPointx",&startPointx);              
  tree->SetBranchAddress("startPointy",&startPointy);              
  tree->SetBranchAddress("startPointz",&startPointz);              
  tree->SetBranchAddress("endPointx",&endPointx);                  
  tree->SetBranchAddress("endPointy",&endPointy);                  
  tree->SetBranchAddress("endPointz",&endPointz);                  
  tree->SetBranchAddress("startT",&startT);                        
  tree->SetBranchAddress("endT",&endT);                            
  tree->SetBranchAddress("pathlen",&pathlen);                      
  tree->SetBranchAddress("numElectrons",&numElectrons);            
  tree->SetBranchAddress("depEnergy",&depEnergy);                  
  tree->SetBranchAddress("process",&process,&br);
  tree->SetBranchAddress("nedeps",&nedeps);                       
  tree->SetBranchAddress("edep_g4id",&edep_g4id);                  
  tree->SetBranchAddress("edep_blipid",&edep_blipid);              
  tree->SetBranchAddress("edep_clustid",&edep_clustid);            
  tree->SetBranchAddress("edep_energy",&edep_energy);              
  tree->SetBranchAddress("edep_charge",&edep_charge);              
  tree->SetBranchAddress("edep_x",&edep_x);                        
  tree->SetBranchAddress("edep_y",&edep_y);                        
  tree->SetBranchAddress("edep_z",&edep_z);                        
  tree->SetBranchAddress("edep_ds",&edep_ds);                      
  tree->SetBranchAddress("nhits",&nhits);                         
  tree->SetBranchAddress("hit_tpc",&hit_tpc);                      
  tree->SetBranchAddress("hit_plane",&hit_plane);                  
  tree->SetBranchAddress("hit_wire",&hit_wire);                    
  tree->SetBranchAddress("hit_channel",&hit_channel);              
  tree->SetBranchAddress("hit_peakT",&hit_peakT);                  
  tree->SetBranchAddress("hit_time",&hit_time);                    
  tree->SetBranchAddress("hit_rms",&hit_rms);                      
  tree->SetBranchAddress("hit_ph",&hit_ph);                        
  tree->SetBranchAddress("hit_area",&hit_area);                    
  tree->SetBranchAddress("hit_charge",&hit_charge);                
  tree->SetBranchAddress("hit_isreal",&hit_isreal);                
  tree->SetBranchAddress("hit_trkid",&hit_trkid);                  
  tree->SetBranchAddress("hit_g4id",&hit_g4id);                    
  tree->SetBranchAddress("hit_g4frac",&hit_g4frac);                
  tree->SetBranchAddress("hit_g4energy",&hit_g4energy);            
  tree->SetBranchAddress("hit_g4charge",&hit_g4charge);            
  tree->SetBranchAddress("hit_clustid",&hit_clustid);              
  tree->SetBranchAddress("nclusts",&nclusts);                     
  tree->SetBranchAddress("clust_tpc",&clust_tpc);                  
  tree->SetBranchAddress("clust_plane",&clust_plane);              
  tree->SetBranchAddress("clust_wire",&clust_wire);                
  tree->SetBranchAddress("clust_nwires",&clust_nwires);            
  tree->SetBranchAddress("clust_nhits",&clust_nhits);              
  tree->SetBranchAddress("clust_charge",&clust_charge);            
  tree->SetBranchAddress("clust_time",&clust_time);                
  tree->SetBranchAddress("clust_startTime",&clust_startTime);      
  tree->SetBranchAddress("clust_endTime",&clust_endTime);          
  tree->SetBranchAddress("clust_xpos",&clust_xpos);          
  tree->SetBranchAddress("clust_g4charge",&clust_g4charge);        
  tree->SetBranchAddress("clust_g4energy",&clust_g4energy);        
  tree->SetBranchAddress("clust_ismatched",&clust_ismatched);      
  tree->SetBranchAddress("clust_edepid",&clust_edepid);            
  tree->SetBranchAddress("clust_blipid",&clust_blipid);            
  tree->SetBranchAddress("nblips",&nblips);                       
  tree->SetBranchAddress("blip_tpc",&blip_tpc);                    
  tree->SetBranchAddress("blip_bestplane",&blip_bestplane);        
  tree->SetBranchAddress("blip_ncrossings",&blip_ncrossings);      
  tree->SetBranchAddress("blip_x",&blip_x);                        
  tree->SetBranchAddress("blip_y",&blip_y);                        
  tree->SetBranchAddress("blip_z",&blip_z);                        
  tree->SetBranchAddress("blip_prms",&blip_prms);                  
  tree->SetBranchAddress("blip_charge",&blip_charge);              
  tree->SetBranchAddress("blip_edepid",&blip_edepid);              
}

void printG4Particle(int i) {
  printf("  %5i  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, dL=%7.1f, KE0=%9.3f, Edep=%9.3f, T=%8.1f --> %8.1f, moth=%5i, %12s, ND=%i\n",
  i,
  trackID[i],
  pdg[i],
  startPointx[i],
  startPointy[i],
  startPointz[i],
  pathlen[i],
  E[i]-mass[i],
  depEnergy[i],
  startT[i],
  endT[i],
  mother[i],
  process->at(i).c_str(),
  nDaughters[i]
  );
}
