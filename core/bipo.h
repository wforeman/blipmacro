
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
    double p0, p0_err, p1, p1_err;
    double rate_signal     = -9;
    double rate_signal_err = 0;
    double rate_bg         = -9;
    double rate_bg_err     = 0;
    double ratio           = -9;
    double activity        = -9;
    double activity_err        = -9;
  };
  
  // useful structure to save candidate info in
  struct BiPoCandidate { int blipID, id1, id2; float dT, q1, q2; };
