
  // NOTES:
  //
  // Mike Z.'s elog entry 126161: bypass started 7/29/21 16:33 --> 1627594380
  //
  // Run 3 data//set: 
  //  T0 = 1528520000 (9 June 2018, 23:53)
  //  T1 = 1531760000 (16 July 2018, 11:53)
  
  TRandom2* fRand;
  float fRecomb       = 0.584; 
  
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
    double N_signal    = -9;
    double N_signal_err = 0;
    double N_bg         = -9;
    double N_bg_err     = 0;
    double ratio           = -9;
    double activity        = -9;
    double activity_err    = -9;
    double activity_err_stat = -9;
    double activity_err_syst = -9;
  };
  
  // useful structure to save candidate info in
  struct BiPoCandidate { 
    int blipID, id1, id2; 
    float dT;
    float q1, q2;
    float e1, e2;
  };

  //##########################################################################
  // Set default ROOT plot style
  //##########################################################################
  void setRootStyle() {
    
    // set margin sizes
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    // set title offsets (for axis label)
    gStyle->SetTitleXOffset(1.2);
    gStyle->SetTitleYOffset(1.3);

    // use large fonts
    //Int_t font=42; // Helvetica
    Double_t tsize=0.05;
    Double_t lsize=0.05;
    gStyle->SetTitleSize(tsize,"x");
    gStyle->SetLabelSize(lsize,"x");
    gStyle->SetTitleSize(tsize,"y");
    gStyle->SetLabelSize(lsize,"y");
    gStyle->SetTitleSize(tsize,"z");
    gStyle->SetLabelSize(lsize,"z");
    gStyle->SetTextSize(tsize);
  
    // use bold lines and markers
    //gStyle->SetLineWidth(2);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(0.7);
    //gStyle->SetMarkerColor(kAzure-6);
    gStyle->SetMarkerColor(kBlack);
    gStyle->SetLineColor(kBlack);

    // get rid of X error bars 
    //gStyle->SetErrorX(0.001);

    // do not display any of the standard histogram decorations
    //atlasStyle->SetOptTitle(0);
    //atlasStyle->SetOptStat(1111);
    //atlasStyle->SetOptStat(0);
    //atlasStyle->SetOptFit(1111);
    //atlasStyle->SetOptFit(0);

    // put tick marks on top and RHS of plots
    //atlasStyle->SetPadTickX(1);
    //atlasStyle->SetPadTickY(1);
  }

  //##########################################################################
  float charge_to_energy(float q, float recomb){
    return (23.6e-6)*(q/recomb);
  }

  
