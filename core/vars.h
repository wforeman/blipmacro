
  // Set max array sizes
	const int kMaxHits  = 100000;
	const int kMaxBlips = 10000;
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
  int   nparticles;               // number of G4 particles
  bool  part_isPrimary[kMaxG4];        // is primary particle
  int   part_trackID[kMaxG4];          // G4 track ID
  int   part_pdg[kMaxG4];              // PDG
  int   part_nDaughters[kMaxG4];       // number of daughters
  int   part_mother[kMaxG4];           // mother particle
  float part_E[kMaxG4];                // initial energy (MeV)
  float part_KE[kMaxG4];               // initial kinetic energy (MeV)
  float part_endE[kMaxG4];             // final energy (MeV)
  float part_endKE[kMaxG4];             // final energy (MeV)
  float part_mass[kMaxG4];             // mass (MeV)
  float part_P[kMaxG4];                // momentum (MeV)
  float part_Px[kMaxG4];               // momentum x (MeV)
  float part_Py[kMaxG4];               // momentum y (MeV)
  float part_Pz[kMaxG4];               // momentum z (MeV)
  float part_startPointx[kMaxG4];      // starting x (cm)
  float part_startPointy[kMaxG4];      // starting y (cm)
  float part_startPointz[kMaxG4];      // starting y (cm)
  float part_endPointx[kMaxG4];        // ending x (cm)
  float part_endPointy[kMaxG4];        // ending y (cm)
  float part_endPointz[kMaxG4];        // ending y (cm)
  float part_startT[kMaxG4];           // starting time (us)
  float part_endT[kMaxG4];             // ending time (us)
  float part_pathlen[kMaxG4];          // path length (cm)
  float part_depEnergy[kMaxG4];        // energy deposited in AV (MeV)
  int   part_depElectrons[kMaxG4];     // electrons deposited
  float part_numElectrons[kMaxG4];     // electrons reaching anode wires
  std::vector<std::string> *part_process;// process name
  //float total_depEnergy;          // total deposited energy in AV
  //int   total_depElectrons;       // total deposited ionization electrons in AV
  //float total_numElectrons;       // total electrons reaching anode wires

  // --- True energy deposit info (derived from SimChannels and SimEnergyDeposits) ---
  int   nedeps;                   // number of true localized energy depositions
  int   edep_tpc[kMaxEDeps];      // TPC
  int   edep_g4id[kMaxEDeps];     // leading G4 track ID
  int   edep_g4index[kMaxEDeps];  // leading G4 track index
  int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
  int   edep_clustid[kMaxEDeps];  // hitclust ID
  int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
  float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
  int   edep_depne[kMaxEDeps];    // total ionization electrons deposited
  float edep_charge[kMaxEDeps];   // total electrons reaching anode wires
  float edep_x[kMaxEDeps];        // x (cm)
  float edep_y[kMaxEDeps];        // y (cm)
  float edep_z[kMaxEDeps];        // z (cm)

  // --- Hit information ---
  int	  nhits;                    // number of hits
  int	  hit_tpc[kMaxHits];        // tpc number
  int	  hit_plane[kMaxHits];      // plane number
  int	  hit_wire[kMaxHits];       // wire number
  int	  hit_channel[kMaxHits];    // channel ID
  float	hit_peakT[kMaxHits];      // raw peak time (tick)
  float	hit_time[kMaxHits];       // corrected peak time (tick)
  float hit_rms[kMaxHits];        // shape RMS
  float	hit_amp[kMaxHits];         // amplitude
  float	hit_area[kMaxHits];       // charge (area) in ADC units
  float hit_sumadc[kMaxHits];     // summed ADC
  float hit_charge[kMaxHits];     // reconstructed number of electrons
  int   hit_mult[kMaxHits];       // multiplicity
  int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
  int   hit_ismatch[kMaxHits];    // does hit have time match on another plane?
  int	  hit_g4id[kMaxHits];       // G4 TrackID of leading particle
  float hit_g4frac[kMaxHits];     // fraction of hit energy from leading MCParticle
  float hit_g4energy[kMaxHits];   // true energy
  float hit_g4charge[kMaxHits];   // true number of electrons (drift-attenuated)
  int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
  int   hit_blipid[kMaxHits];     // key of Blip in which hit was included
  float hit_gof[kMaxHits];        // goodness of fit (default -1)

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

  // --- Hit cluster information ---
  int   nclusts;                      // total clusters made
  int   clust_id[kMaxHits];           // cluster ID (index)
  int   clust_tpc[kMaxHits];          // cluster TPC ID
  int   clust_plane[kMaxHits];        // cluster plane
  int   clust_wire[kMaxHits];         // central-most wire of cluster
  int   clust_startwire[kMaxHits];    // starting wire
  int   clust_endwire[kMaxHits];      // ending wire
  int   clust_nwires[kMaxHits];       // number of wires in this cluster
  int   clust_nhits[kMaxHits];        // number of hits
  float clust_time[kMaxHits];         // charge-weighted time
  float clust_timespan[kMaxHits];     // cluster timespan
  float clust_rms[kMaxHits];          // charge-weighted RMS
  float clust_starttime[kMaxHits];    // cluster start tick
  float clust_endtime[kMaxHits];      // cluster end tick
  float clust_amp[kMaxHits];          // maximum hit amplitude [ADC]
  float clust_charge[kMaxHits];       // cluster charge at anode [e-]
  float clust_g4charge[kMaxHits];     // true cluster charge at anode
  float clust_g4energy[kMaxHits];     // true cluster energy from G4
  int   clust_g4id[kMaxHits];         // true MCParticle ID (index for particle branches)
  int   clust_blipid[kMaxHits];       // blip ID for this nlusteer (if it was made into one)
  int   clust_edepid[kMaxHits];       // true energy dep ID
  bool  clust_ismatch[kMaxHits];      // was this cluster plane-matched?
  //int   clust_lhit_wire[kMaxHits];    // cluster wire (lead hit wire)
  //int   clust_lhit_chan[kMaxHits];    // cluster channel (lead hit wire)
  //int   clust_lhit_id[kMaxHits];      // lead hit ID (index for hit_X[i] branches)
  //float clust_lhit_amp[kMaxHits];     // lead hit peak amplitude [ADC]
  //float clust_lhit_rms[kMaxHits];     // lead hit RMS [ADC]
  //float clust_lhit_time[kMaxHits];    // lead hit time [ticks]
  //float clust_lhit_gof[kMaxHits];     // lead hit goodness-of-fit; pulse train = -1
  //bool  clust_lhit_isfit[kMaxHits];   // is there a valid goodness of fit for lead hit?

  // --- 3D Blip information ---
  int   nblips;                       // number of blips in event
  int   blip_id[kMaxBlips];           // blip ID / index
  int   blip_tpc[kMaxBlips];          // blip TPC
  int   blip_nplanes[kMaxBlips];      // number of planes matched (2 or 3)
  float blip_x[kMaxBlips];            // X position [cm]
  float blip_y[kMaxBlips];            // Y position [cm]
  float blip_z[kMaxBlips];            // Z position [cm]
  float blip_sigmayz[kMaxBlips];      // difference in wire intersection points
  float blip_dx[kMaxBlips];           // dX [cm]
  float blip_dyz[kMaxBlips];
  float blip_sumadc[kMaxBlips];          // integrated ADCs 
  float blip_charge[kMaxBlips];       // blip charge at anode [e-]
  float blip_energy[kMaxBlips];       // blip energy [MeV]
  int   blip_edepid[kMaxBlips];       // true energy dep ID
  float blip_trkdist[kMaxBlips];      // distance to nearest track
  int   blip_trkid[kMaxBlips];        // index of nearest trk
  bool  blip_incylinder[kMaxBlips];   // is blip within a cylinder near a track
  float blip_time[kMaxBlips];                   // blip drift time [ticks]
  int   blip_clustid[kNplanes][kMaxBlips];     // cluster ID per plane

//################################
// Print G4 particle
//################################
void printG4Particle(int i) {
  printf("  %5i  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, dL=%7.1f, KE0=%9.3f, Edep=%9.3f, T=%8.1f --> %8.1f, moth=%5i, %12s, ND=%i\n",
  i,
  part_trackID[i],
  part_pdg[i],
  part_startPointx[i],
  part_startPointy[i],
  part_startPointz[i],
  part_pathlen[i],
  part_KE[i],
  part_depEnergy[i],
  part_startT[i],
  part_endT[i],
  part_mother[i],
  part_process->at(i).c_str(),
  part_nDaughters[i]
  );
}

