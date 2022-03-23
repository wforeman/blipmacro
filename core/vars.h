
  // Set global constants

  // Set max array sizes
	const int kMaxHits  = 5000;
	const int kMaxBlips = 1000;
	const int kMaxTrks  = 100;
	const int kMaxG4    = 10000;
	const int kMaxEDeps = 10000;
	const int kNplanes  = 3;  
  const int kMaxShwrs = 1000;    
  

  // --- Event information ---   
  int           event;                    // event number
  int           run;                      // run number
  unsigned int  timestamp;                // unix time of event
  float         lifetime;                 // electron lifetime

  // --- G4 information ---
  float total_depEnergy;          // total deposited energy in AV
  float total_numElectrons;       // total electrons reaching anode wires
  float gamma_depEnergy;          // total gamma-induced energy deposited
  float gamma_numElectrons;       // total electrons from gamma-induced depositions
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
  std::vector<std::string> *process;// process name

// --- True energy deposit info (derived) ---
  int   nedeps;                   // number of true localized energy depositions
  int   edep_tpc[kMaxEDeps];      // TPC
  int   edep_g4id[kMaxEDeps];     // leading G4 track ID
  int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
  int   edep_clustid[kMaxEDeps];  // hitclust ID
  int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
  float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
  float edep_charge[kMaxEDeps];   // total electrons reaching anode wires
  float edep_x[kMaxEDeps];        // x (cm)
  float edep_y[kMaxEDeps];        // y (cm)
  float edep_z[kMaxEDeps];        // z (cm)
  float edep_ds[kMaxEDeps];       // extent (cm)

  // --- Hit information ---
  int   nhits;                    // number of hits
  int   hit_tpc[kMaxHits];        // tpc number
  int   hit_plane[kMaxHits];      // plane number
  int   hit_wire[kMaxHits];       // wire number
  int   hit_channel[kMaxHits];    // channel ID
  float hit_peakT[kMaxHits];      // raw peak time (tick)
  float hit_time[kMaxHits];       // corrected peak time (tick)
  float hit_rms[kMaxHits];        // shape RMS
  float hit_amp[kMaxHits];         // amplitude
  float hit_area[kMaxHits];       // charge (area) in ADC units
  float hit_sumadc[kMaxHits];     // summed ADC
  float hit_charge[kMaxHits];     // reconstructed number of electrons
  int   hit_mult[kMaxHits];       // multiplicity
  int   hit_trkid[kMaxHits];      // is this hit associated with a reco track?
  int   hit_shwrid[kMaxHits];      // is this hit associated with a reco shower?
  int   hit_ismatch[kMaxHits];    // does hit have time match on another plane?
  int   hit_isgood[kMaxHits];     // does hit pass shape-based cuts?
  int   hit_isreal[kMaxHits];     // is this hit real?
  //int   hit_isneartrk[kMaxHits];  // is this hit within a track cylinder cut
  int   hit_g4id[kMaxHits];       // G4 TrackID of leading particle
  float hit_g4frac[kMaxHits];     // fraction of hit energy from leading MCParticle
  float hit_g4energy[kMaxHits];   // true energy
  float hit_g4charge[kMaxHits];   // true number of electrons (drift-attenuated)
  int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
  int   hit_blipid[kMaxHits];     // key of Blip in which hit was included
  float hit_gof[kMaxHits];   // goodness of fit (default -1)

  // --- Track hits ---
  float ave_trkhit_amp;            // average pulse height for hits in tracks

	// --- Track information ---
  int   ntrks;                    // number tracks
  int   trk_id[kMaxTrks];         // trackID
  int   trk_npts[kMaxTrks];       // number 3D trajectory points
  float trk_length[kMaxTrks];     // track length [cm]
  float trk_startx[kMaxTrks];     // starting X coordinate
  float trk_starty[kMaxTrks];     // starting Y coordinate
  float trk_startz[kMaxTrks];     // starting Z coordinate
  float trk_startd[kMaxTrks];     // starting distance to boundary
  float trk_endx[kMaxTrks];       // ending X coordinate
  float trk_endy[kMaxTrks];       // ending Y coordinate
  float trk_endz[kMaxTrks];       // ending Z coordinate
  float trk_endd[kMaxTrks];       // ending distance to boundary
  //int   trk_origin[kMaxTrks];     // 0: unknown, 1: cosmic, 2: neutrino, 3: supernova, 4: singles
  //int   trk_g4id[kMaxTrks];       // G4 TrackID of leading contributing particle
  //int   trk_g4pdg[kMaxTrks];      // G4 PDG
  //float trk_purity[kMaxTrks];     // track hit purity
  //float trk_pitch[kMaxTrks];      // track pitch
  //float trk_ke[kMaxTrks];         // track kinetic energy
  //int   trk_cosmictag[kMaxTrks];  // cosmic tagg
  //float trk_cosmicscore[kMaxTrks];// cosmic score
  //int   trk_pidpdg[kMaxTrks];     // particle PID Pdg
  //float trk_pidchi[kMaxTrks];     // chisquared for PID
  //int   trk_bestplane[kMaxTrks];  // plane w/ most hits

  // --- Shower information ---
  /*
  int   nshwrs;                     // number showers
  int   shwr_id[kMaxShwrs];         // shower ID
  float shwr_dirx[kMaxShwrs];       // direction X
  float shwr_diry[kMaxShwrs];       // direction Y
  float shwr_dirz[kMaxShwrs];       // direction Z
  float shwr_startx[kMaxShwrs];     // starting X coordinate
  float shwr_starty[kMaxShwrs];     // starting X coordinate
  float shwr_startz[kMaxShwrs];     // starting X coordinate
  float shwr_length[kMaxShwrs];     // shower length
  float shwr_openangle[kMaxShwrs];  // opening angle
  */

	// --- Hit cluster information ---
  int   nclusts;
  int   clust_tpc[kMaxHits];
  int   clust_plane[kMaxHits];
  int   clust_wire[kMaxHits];
  int   clust_chan[kMaxHits];
  int   clust_nwires[kMaxHits];
  int   clust_nhits[kMaxHits];
  int   clust_lhit_id[kMaxHits];
  float clust_lhit_amp[kMaxHits];
  float clust_lhit_rms[kMaxHits];
  float clust_lhit_time[kMaxHits];
  float clust_lhit_peakT[kMaxHits];
  float clust_lhit_gof[kMaxHits];
  float clust_charge[kMaxHits];
  float clust_time[kMaxHits];
  float clust_time_w[kMaxHits];
  float clust_time_err[kMaxHits];
  float clust_startTime[kMaxHits];
  float clust_endTime[kMaxHits];
  float clust_timespan[kMaxHits];
  float clust_g4energy[kMaxHits];
  float clust_g4charge[kMaxHits];
  int   clust_g4id[kMaxHits];
  int   clust_ismatch[kMaxHits];
  //int   clust_isneartrk[kMaxHits];
  int   clust_blipid[kMaxHits];
  int   clust_edepid[kMaxHits];

  // --- Blip information ---
  float total_blip_energy;
  int   nblips;
  int   blip_tpc[kMaxBlips];
  int   blip_nplanes[kMaxBlips];
  float blip_x[kMaxBlips];
  float blip_y[kMaxBlips];
  float blip_z[kMaxBlips];
  float blip_maxdiff[kMaxBlips];
  float blip_charge[kMaxBlips];
  float blip_energy[kMaxBlips];
  int   blip_edepid[kMaxBlips];
  int   blip_clustid[kNplanes][kMaxBlips];
  int   blip_clustid_pl2[kMaxBlips];
  float blip_trkdist[kMaxBlips];


//################################
// Print G4 particle
//################################
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

