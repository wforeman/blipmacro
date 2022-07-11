
  // Struct for input file parameters
  struct infile_t { 
    std::string   fileName;
    std::string   treeName;
    bool          isMC;
    unsigned int  t0;
    unsigned int  t1;
  };
  
  // Struct to hold fit results
  struct FitResult { 
    float p0, p0_err, p1, p1_err;
    float rate_signal     = -9;
    float rate_signal_err = 0;
    float rate_bg         = -9;
    float rate_bg_err     = 0;
    float ratio           = -9;
    float activity        = -9;
    float activity_err        = -9;
  };
  
  // useful structure to save candidate info in
  struct BiPoCandidate { int blipID, id1, id2; float dT, q1, q2; };
