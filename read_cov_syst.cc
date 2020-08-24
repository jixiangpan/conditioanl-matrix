void read_cov_syst()
{
  gROOT->ProcessLine(".x ../lhcbStyle_edit.C");
  TString roostr = "";

  // TGaxis *yaxis = new TGaxis();
  // yaxis->SetMaxDigits(3);
  // https://root.cern/doc/v610/classTHistPainter.html
  // gStyle->SetPalette(58);
  
  //////////////////////////////////////////////////////////////////////// global variables

  const int Nnue = 11;
  const int Nnumu = 9;
  const int Ntotal = 20;

  double array_err_syst[Ntotal] = {0};
  double array_err_stat[Ntotal] = {0};
  double array_err_full[Ntotal] = {0};
  double array_val_mean[Ntotal] = {0};
  
  TMatrixD matrix_syst_covariance(Ntotal, Ntotal);
  TMatrixD matrix_stat_covariance(Ntotal, Ntotal);
  TMatrixD matrix_full_covariance(Ntotal, Ntotal);

  double full_covariance_min = 7e4;
  double full_covariance_max = 0.11e-1;
  
  //////////////////////////////////////////////////////////////////////// spectra of nue and numu
  
  TFile *input_spectra = new TFile("DLLEECovar.SBNspec.root", "read");

  roostr = "nu_uBooNE_nue_intrinsic";
  TH1D *nu_uBooNE_nue_intrinsic = (TH1D*)input_spectra->Get(roostr);
  roostr = "nu_uBooNE_nue_cocktail";
  TH1D *nu_uBooNE_nue_cocktail = (TH1D*)input_spectra->Get(roostr);
  roostr = "nu_uBooNE_nue_extbnb";
  TH1D *nu_uBooNE_nue_extbnb = (TH1D*)input_spectra->Get(roostr);

  roostr = "nu_uBooNE_numu_cocktail";
  TH1D *nu_uBooNE_numu_cocktail = (TH1D*)input_spectra->Get(roostr);
  roostr = "nu_uBooNE_numu_extbnb";
  TH1D *nu_uBooNE_numu_extbnb = (TH1D*)input_spectra->Get(roostr);

  TH1D *h1_nue = (TH1D*)nu_uBooNE_nue_intrinsic->Clone("h1_nue");
  h1_nue->Reset();
  h1_nue->Add( nu_uBooNE_nue_intrinsic );
  h1_nue->Add( nu_uBooNE_nue_cocktail );
  h1_nue->Add( nu_uBooNE_nue_extbnb );

  TH1D *h1_numu = (TH1D*)nu_uBooNE_numu_cocktail->Clone("h1_numu");
  h1_numu->Reset();
  h1_numu->Add( nu_uBooNE_numu_cocktail );
  h1_numu->Add( nu_uBooNE_numu_extbnb );

  for(int ibin=1; ibin<=Nnue; ibin++) {
    double value = h1_nue->GetBinContent(ibin);
    h1_nue->SetBinError( ibin, sqrt(value) );
  }

  for(int ibin=1; ibin<=Nnumu; ibin++) {
    double value = h1_numu->GetBinContent(ibin);
    h1_numu->SetBinError( ibin, sqrt(value) );
  }

  //////////////////////////////////////////////////////////////////////// covaraince matrix of flux and cross-section

  TFile *input_matrix = new TFile("DLLEECovar.SBNcollapsedcovar.root", "read");
  TMatrixD *syst_matrix = (TMatrixD*)input_matrix->Get("collapsed_covariance_DLLEECovar");

  for(int i=0; i<Ntotal; i++) {// systematics
    for(int j=0; j<Ntotal; j++) {
      double error = syst_matrix(i,j);
      matrix_syst_covariance(i,j) = error;
    }
  }
  
  for(int i=1; i<=Nnue; i++) {// statistics nue, 0,1,2,3,4,5,6,7,8,9,10
    double value = h1_nue->GetBinContent(i);
    matrix_stat_covariance(i-1,i-1) = value;
  }

  for(int i=1; i<=Nnumu; i++) {// statistics numu
    double value = h1_numu->GetBinContent(i);
    matrix_stat_covariance(10+i,10+i) = value;
  }

  matrix_full_covariance = matrix_syst_covariance + matrix_stat_covariance;// full covarinace matrix

  //////////////////////////////////////////////////////////////////////// checking and ploting

  roostr = "h2_check_syst_matrix";
  TH2D *h2_check_syst_matrix = new TH2D(roostr, roostr, Ntotal, 0.5, Ntotal+0.5, Ntotal, 0.5, Ntotal+0.5);
 
  roostr = "h2_check_stat_matrix";
  TH2D *h2_check_stat_matrix = new TH2D(roostr, roostr, Ntotal, 0.5, Ntotal+0.5, Ntotal, 0.5, Ntotal+0.5);
 
  roostr = "h2_check_full_matrix";
  TH2D *h2_check_full_matrix = new TH2D(roostr, roostr, Ntotal, 0.5, Ntotal+0.5, Ntotal, 0.5, Ntotal+0.5);
 
  roostr = "h2_check_correlation_syst_matrix";
  TH2D *h2_check_correlation_syst_matrix = new TH2D(roostr, roostr, Ntotal, 0.5, Ntotal+0.5, Ntotal, 0.5, Ntotal+0.5);
  
  roostr = "h2_check_correlation_full_matrix";
  TH2D *h2_check_correlation_full_matrix = new TH2D(roostr, roostr, Ntotal, 0.5, Ntotal+0.5, Ntotal, 0.5, Ntotal+0.5);
  
  for(int i=0; i<Ntotal; i++) {
    for(int j=0; j<Ntotal; j++) {
      double value = 0;
      
      value = matrix_syst_covariance(i,j);
      h2_check_syst_matrix->SetBinContent(i+1, j+1, value);

      value = matrix_stat_covariance(i,j);
      h2_check_stat_matrix->SetBinContent(i+1, j+1, (int(value*100+0.5))*1./100. );

      value = matrix_full_covariance(i,j);
      h2_check_full_matrix->SetBinContent(i+1, j+1, value);
    }
  }

  for(int i=0; i<Ntotal; i++) {
    for(int j=0; j<Ntotal; j++) {
      double value    = matrix_syst_covariance(i,j);
      double sigma2_i = matrix_syst_covariance(i,i);
      double sigma2_j = matrix_syst_covariance(j,j);
      h2_check_correlation_syst_matrix->SetBinContent(i+1, j+1, value/sqrt(sigma2_i*sigma2_j) );

      value    = matrix_full_covariance(i,j);
      sigma2_i = matrix_full_covariance(i,i);
      sigma2_j = matrix_full_covariance(j,j);
      h2_check_correlation_full_matrix->SetBinContent(i+1, j+1, value/sqrt(sigma2_i*sigma2_j) );
    }
  }

  //////
  roostr = "canv_h2_check_syst_matrix";
  TCanvas *canv_h2_check_syst_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_syst_matrix->SetLeftMargin(0.12);
  canv_h2_check_syst_matrix->SetRightMargin(0.14);
  canv_h2_check_syst_matrix->SetTopMargin(0.09);
  canv_h2_check_syst_matrix->SetBottomMargin(0.18);
  canv_h2_check_syst_matrix->SetLogz();
  h2_check_syst_matrix->Draw("colz");
  h2_check_syst_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_syst_matrix->GetYaxis()->SetNdivisions(310);
  h2_check_syst_matrix->SetXTitle("Bin index");
  h2_check_syst_matrix->SetYTitle("Bin index");
  h2_check_syst_matrix->GetXaxis()->CenterTitle();
  h2_check_syst_matrix->GetYaxis()->CenterTitle();
  
  //////
  roostr = "canv_h2_check_stat_matrix";
  TCanvas *canv_h2_check_stat_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_stat_matrix->SetLeftMargin(0.12);
  canv_h2_check_stat_matrix->SetRightMargin(0.14);
  canv_h2_check_stat_matrix->SetTopMargin(0.09);
  canv_h2_check_stat_matrix->SetBottomMargin(0.18);
  canv_h2_check_stat_matrix->SetLogz();
  h2_check_stat_matrix->Draw("colz text");
  h2_check_stat_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_stat_matrix->GetYaxis()->SetNdivisions(310);

  //////
  roostr = "canv_h2_check_full_matrix";
  TCanvas *canv_h2_check_full_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_full_matrix->SetLeftMargin(0.12);
  canv_h2_check_full_matrix->SetRightMargin(0.14);
  canv_h2_check_full_matrix->SetTopMargin(0.09);
  canv_h2_check_full_matrix->SetBottomMargin(0.18);
  canv_h2_check_full_matrix->SetLogz();
  h2_check_full_matrix->Draw("colz");
  h2_check_full_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_full_matrix->GetYaxis()->SetNdivisions(310);

  //////
  roostr = "canv_h2_check_correlation_syst_matrix";
  TCanvas *canv_h2_check_correlation_syst_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_correlation_syst_matrix->SetLeftMargin(0.12);
  canv_h2_check_correlation_syst_matrix->SetRightMargin(0.14);
  canv_h2_check_correlation_syst_matrix->SetTopMargin(0.09);
  canv_h2_check_correlation_syst_matrix->SetBottomMargin(0.18);
  // canv_h2_check_correlation_syst_matrix->SetLogz();
  h2_check_correlation_syst_matrix->Draw("colz");
  h2_check_correlation_syst_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_correlation_syst_matrix->GetYaxis()->SetNdivisions(310);
  h2_check_correlation_syst_matrix->GetZaxis()->SetRangeUser(0.01e-1,1);
  // h2_check_correlation_syst_matrix->SetXTitle("Bin index");
  // h2_check_correlation_syst_matrix->SetYTitle("Bin index");
  h2_check_correlation_syst_matrix->GetXaxis()->CenterTitle();
  h2_check_correlation_syst_matrix->GetYaxis()->CenterTitle();
  canv_h2_check_correlation_syst_matrix->SaveAs("canv_h2_check_correlation_syst_matrix.png");
  
  //////
  roostr = "canv_h2_check_correlation_full_matrix";
  TCanvas *canv_h2_check_correlation_full_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_correlation_full_matrix->SetLeftMargin(0.12);
  canv_h2_check_correlation_full_matrix->SetRightMargin(0.14);
  canv_h2_check_correlation_full_matrix->SetTopMargin(0.09);
  canv_h2_check_correlation_full_matrix->SetBottomMargin(0.18);
  // canv_h2_check_correlation_full_matrix->SetLogz();
  h2_check_correlation_full_matrix->Draw("colz");
  h2_check_correlation_full_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_correlation_full_matrix->GetYaxis()->SetNdivisions(310);
  h2_check_correlation_full_matrix->GetZaxis()->SetRangeUser(0.01e-1,1);
  // h2_check_correlation_full_matrix->SetXTitle("Bin index");
  // h2_check_correlation_full_matrix->SetYTitle("Bin index");
  h2_check_correlation_full_matrix->GetXaxis()->CenterTitle();
  h2_check_correlation_full_matrix->GetYaxis()->CenterTitle();
  canv_h2_check_correlation_full_matrix->SaveAs("canv_h2_check_correlation_full_matrix.png");

  //////////
  cout<<endl;
  for(int i=0; i<Ntotal; i++) {
    double err_stat = matrix_stat_covariance(i,i);
    double err_syst = matrix_syst_covariance(i,i);
    double err_full = matrix_full_covariance(i,i);
    double val_mean = err_stat;

    array_err_syst[i] = err_syst;
    array_err_stat[i] = err_stat;
    array_err_full[i] = err_full;
    array_val_mean[i] = val_mean;

    cout<<TString::Format("%2d, %10.3f(relstat %5.3f), %10.3f(relsyst %5.3f), %10.3f(relsum %5.3f), %10.3f(relful %5.3f)",
			  i+1,
			  array_err_stat[i], sqrt(array_err_stat[i])/array_val_mean[i],
			  array_err_syst[i], sqrt(array_err_syst[i])/array_val_mean[i],
			  array_err_stat[i]+array_err_syst[i], sqrt(array_err_syst[i]+array_err_stat[i])/array_val_mean[i],
			  array_err_full[i], sqrt(array_err_full[i])/array_val_mean[i]
			  )<<endl;
  }
  cout<<endl;
  
  
  //////////////////////////////////////////////////////////////////////// checking and ploting

  TH1D *h1_nue_errSyst = (TH1D*)h1_nue->Clone("h1_nue_errSyst");
  for(int ibin=1; ibin<=Nnue; ibin++) {
    double value = h1_nue->GetBinContent(ibin);
    double err_stat = sqrt(value);
    double err_syst = h2_check_syst_matrix->GetBinContent(ibin,ibin);
    err_syst = sqrt(err_syst);
    h1_nue_errSyst->SetBinError(ibin, err_syst);
  }

  TH1D *h1_numu_errSyst = (TH1D*)h1_numu->Clone("h1_numu_errSyst");
  for(int ibin=1; ibin<=Nnumu; ibin++) {
    double value = h1_numu->GetBinContent(ibin);
    double err_stat = sqrt(value);
    double err_syst = h2_check_syst_matrix->GetBinContent(ibin+11,ibin+11);
    err_syst = sqrt(err_syst);
    h1_numu_errSyst->SetBinError(ibin, err_syst);
  }

  //////
  roostr = "canv_spectrum_nue";
  TCanvas *canv_spectrum_nue = new TCanvas(roostr, roostr, 900, 650);
  canv_spectrum_nue->SetRightMargin(0.05);
  canv_spectrum_nue->SetTopMargin(0.09);
  canv_spectrum_nue->SetBottomMargin(0.18);
  canv_spectrum_nue->SetLeftMargin(0.12);
  h1_nue->Draw("hist");
  h1_nue->SetMaximum( h1_nue->GetMaximum(  )*1.5 );
  h1_nue->SetLineColor(kBlack);
  h1_nue->SetMarkerColor(kBlack);
  h1_nue->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_nue->GetXaxis()->CenterTitle();
  h1_nue->GetYaxis()->CenterTitle();
  h1_nue->GetYaxis()->SetNdivisions(508);
  h1_nue->SetYTitle("Entries");
  h1_nue->GetXaxis()->SetLabelSize(0.07);
  h1_nue->GetXaxis()->SetTitleSize(0.07);
  h1_nue->GetYaxis()->SetLabelSize(0.07);
  h1_nue->GetYaxis()->SetTitleSize(0.07);  
  h1_nue->GetXaxis()->SetTitleOffset(1.18);
  h1_nue->GetYaxis()->SetTitleOffset(0.85);

  h1_nue_errSyst->Draw("same E2");
  h1_nue_errSyst->SetFillColor(1);
  h1_nue_errSyst->SetFillStyle(3005);
  h1_nue_errSyst->SetMarkerStyle(1);
  h1_nue->Draw("same hist");
  h1_nue->Draw("same axis");
  
  //////
  roostr = "canv_spectrum_numu";
  TCanvas *canv_spectrum_numu = new TCanvas(roostr, roostr, 900, 650);
  canv_spectrum_numu->SetRightMargin(0.05);
  canv_spectrum_numu->SetTopMargin(0.09);
  canv_spectrum_numu->SetBottomMargin(0.18);
  canv_spectrum_numu->SetLeftMargin(0.14);
  h1_numu->Draw("hist");
  h1_numu->SetMaximum( h1_numu->GetMaximum(  )*1.5 );
  h1_numu->SetLineColor(kBlack);
  h1_numu->SetMarkerColor(kBlack);
  h1_numu->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_numu->GetXaxis()->CenterTitle();
  h1_numu->GetYaxis()->CenterTitle();
  h1_numu->GetXaxis()->SetNdivisions(508);
  h1_numu->GetYaxis()->SetNdivisions(507);
  h1_numu->SetYTitle("Entries");
  h1_numu->GetXaxis()->SetLabelSize(0.07);
  h1_numu->GetXaxis()->SetTitleSize(0.07);
  h1_numu->GetYaxis()->SetLabelSize(0.07);
  h1_numu->GetYaxis()->SetTitleSize(0.07);  
  h1_numu->GetXaxis()->SetTitleOffset(1.18);
  h1_numu->GetYaxis()->SetTitleOffset(1.05);

  h1_numu_errSyst->Draw("same E2");
  h1_numu_errSyst->SetFillColor(1);
  h1_numu_errSyst->SetFillStyle(3005);
  h1_numu_errSyst->SetMarkerStyle(1);
  h1_numu->Draw("same hist");
  h1_numu->Draw("same axis");

  //////////
 
  roostr = "h1_nue_stat";
  TH1D *h1_nue_stat = new TH1D(roostr, roostr, Nnue, 100, 1200);
  
  cout<<endl;
  for(int i=0; i<Nnue; i++) {
    cout<<TString::Format("%2d %5.2f (stat/full)",
			  i+1, array_err_stat[i]/array_err_full[i]
			  )<<endl;

    h1_nue_stat->SetBinContent(i+1, array_err_stat[i]/array_err_full[i]);
  }
  cout<<endl;

  roostr = "canv_h1_nue_stat";
  TCanvas *canv_h1_nue_stat = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_nue_stat->SetRightMargin(0.05);
  canv_h1_nue_stat->SetTopMargin(0.09);
  canv_h1_nue_stat->SetBottomMargin(0.18);
  canv_h1_nue_stat->SetLeftMargin(0.12);
  h1_nue_stat->Draw("hist");
  h1_nue_stat->SetLineColor(kBlack);
  h1_nue_stat->SetMarkerColor(kBlack);
  h1_nue_stat->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_nue_stat->GetXaxis()->CenterTitle();
  h1_nue_stat->GetYaxis()->CenterTitle();
  h1_nue_stat->GetYaxis()->SetNdivisions(508);
  h1_nue_stat->SetYTitle("");
  h1_nue_stat->GetXaxis()->SetLabelSize(0.07);
  h1_nue_stat->GetXaxis()->SetTitleSize(0.07);
  h1_nue_stat->GetYaxis()->SetLabelSize(0.07);
  h1_nue_stat->GetYaxis()->SetTitleSize(0.07);  
  h1_nue_stat->GetXaxis()->SetTitleOffset(1.18);
  h1_nue_stat->GetYaxis()->SetTitleOffset(0.85);
 

  ///////////////////////////////////////////////////////////////////////////////// numu constrains nue
  ////////////////////////////////// case_numu ////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  
  const int case_numu_Nx = 9;
  const int case_numu_Ny = 11;

  double array_case_numu_constraint_full[case_numu_Ny] = {0};
  double array_case_numu_constraint_syst[case_numu_Ny] = {0};
  
  TMatrixD case_numu_convert_matrix( Ntotal, case_numu_Nx+case_numu_Ny );
  TMatrixD case_numu_convert_matrixT( case_numu_Nx+case_numu_Ny, Ntotal );
  
  for(int i=0; i<case_numu_Nx; i++) {// for x
    case_numu_convert_matrix[i+case_numu_Ny][i+0] = 1;// [old-index][new-index]
  }

  for(int i=0; i<case_numu_Ny; i++) {// for y
    case_numu_convert_matrix[i][i+case_numu_Nx] = 1;// [old-index][new-index]
  }

  ////////
  case_numu_convert_matrixT.Transpose( case_numu_convert_matrix );

  TMatrixD case_numu_cov_matrix = case_numu_convert_matrixT * matrix_full_covariance * case_numu_convert_matrix;

  //////
  roostr = "h2_check_case_numu_cov_matrix";
  TH2D *h2_check_case_numu_cov_matrix = new TH2D(roostr, roostr,
					     case_numu_Nx+case_numu_Ny, 0.5, case_numu_Nx+case_numu_Ny+0.5,
					     case_numu_Nx+case_numu_Ny, 0.5, case_numu_Nx+case_numu_Ny+0.5);
  
  for(int i=0; i<case_numu_Nx+case_numu_Ny; i++) {
    for(int j=0; j<case_numu_Nx+case_numu_Ny; j++) {
      double value = case_numu_cov_matrix(i,j);
      h2_check_case_numu_cov_matrix->SetBinContent(i+1, j+1, value );
    }
  }

  roostr = "canv_h2_check_case_numu_cov_matrix";
  TCanvas *canv_h2_check_case_numu_cov_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_case_numu_cov_matrix->SetLeftMargin(0.12);
  canv_h2_check_case_numu_cov_matrix->SetRightMargin(0.14);
  canv_h2_check_case_numu_cov_matrix->SetTopMargin(0.09);
  canv_h2_check_case_numu_cov_matrix->SetBottomMargin(0.18);
  canv_h2_check_case_numu_cov_matrix->SetLogz();
  h2_check_case_numu_cov_matrix->Draw("colz");
  h2_check_case_numu_cov_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_case_numu_cov_matrix->GetYaxis()->SetNdivisions(310);
  
  ////////
  TMatrixD case_numu_cov_xx(case_numu_Nx, case_numu_Nx);
  TMatrixD case_numu_cov_yy(case_numu_Ny, case_numu_Ny);
  TMatrixD case_numu_cov_xy(case_numu_Nx, case_numu_Ny);
  TMatrixD case_numu_cov_yx(case_numu_Ny, case_numu_Nx);
  TMatrixD case_numu_cov_xx_inv(case_numu_Nx, case_numu_Nx);
  TMatrixD case_numu_cov_y_under_x(case_numu_Ny, case_numu_Ny);

  for(int i=0; i<case_numu_Nx; i++) {
    for(int j=0; j<case_numu_Nx; j++) {
      case_numu_cov_xx[i][j] = case_numu_cov_matrix[i][j];
    }
  }

  for(int i=0; i<case_numu_Ny; i++) {
    for(int j=0; j<case_numu_Ny; j++) {
      case_numu_cov_yy[i][j] = case_numu_cov_matrix[i+case_numu_Nx][j+case_numu_Nx];
    }
  }
  
  for(int i=0; i<case_numu_Nx; i++) {
    for(int j=0; j<case_numu_Ny; j++) {
      case_numu_cov_xy[i][j] = case_numu_cov_matrix[i][j+case_numu_Nx];
    }
  }
  
  case_numu_cov_yx.Transpose( case_numu_cov_xy );
  case_numu_cov_xx_inv = case_numu_cov_xx;
  case_numu_cov_xx_inv.Invert();
  case_numu_cov_y_under_x = case_numu_cov_yy - case_numu_cov_yx*case_numu_cov_xx_inv*case_numu_cov_xy;
  
  cout<<endl<<" ---> numu_constraint_full"<<endl;
  for(int ibin=1; ibin<=case_numu_Ny; ibin++) {
    double error = case_numu_cov_y_under_x(ibin-1, ibin-1);
    error = sqrt(error);
    double value = array_val_mean[ibin-1];
    cout<<TString::Format("%2d rel.error %5.3f(full) %5.3f(syst)",
			  ibin,
			  error/value,
			  sqrt( error*error - array_err_stat[ibin-1]  )/value
			  )<<endl;

    array_case_numu_constraint_full[ibin-1] = error*error;
    array_case_numu_constraint_syst[ibin-1] = error*error - array_err_stat[ibin-1];
  }
  cout<<endl;

  //////////////////////////////////////
  roostr = "h1_case_numu_original_full";
  TH1D *h1_case_numu_original_full = (TH1D*)h1_nue->Clone(roostr);
  h1_case_numu_original_full->Reset();
  roostr = "h1_case_numu_constraint_full";
  TH1D *h1_case_numu_constraint_full = (TH1D*)h1_nue->Clone(roostr);
  h1_case_numu_constraint_full->Reset();
  for(int ibin=1; ibin<=case_numu_Ny; ibin++) {
    h1_case_numu_original_full->SetBinContent( ibin, sqrt(array_err_full[ibin-1])/array_val_mean[ibin-1] );
    h1_case_numu_constraint_full->SetBinContent( ibin, sqrt(array_case_numu_constraint_full[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_numu_full";
  TCanvas *canv_case_numu_full = new TCanvas(roostr, roostr, 900, 650);
  canv_case_numu_full->SetRightMargin(0.05);
  canv_case_numu_full->SetTopMargin(0.09);
  canv_case_numu_full->SetBottomMargin(0.18);
  canv_case_numu_full->SetLeftMargin(0.12);

  h1_case_numu_original_full->Draw("hist");
  double case_numu_original_full_max = 0;
  for(int ibin=1; ibin<=case_numu_Ny; ibin++) {
    double content = h1_case_numu_original_full->GetBinContent(ibin);
    if( content>case_numu_original_full_max ) case_numu_original_full_max = content;
  }
  h1_case_numu_original_full->SetMaximum( case_numu_original_full_max*1.2 );
  h1_case_numu_original_full->SetMinimum(0);
  h1_case_numu_original_full->SetLineColor(kRed);
  h1_case_numu_original_full->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_numu_original_full->GetXaxis()->CenterTitle();
  h1_case_numu_original_full->GetYaxis()->CenterTitle();
  h1_case_numu_original_full->GetYaxis()->SetNdivisions(508);
  h1_case_numu_original_full->SetYTitle("Relative uncertainty");
  h1_case_numu_original_full->GetXaxis()->SetLabelSize(0.07);
  h1_case_numu_original_full->GetXaxis()->SetTitleSize(0.07);
  h1_case_numu_original_full->GetYaxis()->SetLabelSize(0.07);
  h1_case_numu_original_full->GetYaxis()->SetTitleSize(0.07);  
  h1_case_numu_original_full->GetXaxis()->SetTitleOffset(1.18);
  h1_case_numu_original_full->GetYaxis()->SetTitleOffset(0.85);

  h1_case_numu_constraint_full->Draw("same");
  h1_case_numu_constraint_full->SetLineColor(kBlue);
  h1_case_numu_constraint_full->Draw("same axis");

  TLegend *lg_case_numu_full = new TLegend(0.154+0.3, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_numu_full->SetBorderSize(1);
  lg_case_numu_full->SetTextFont(42);
  lg_case_numu_full->SetTextSize(0.068);
  lg_case_numu_full->AddEntry(h1_case_numu_original_full,   "no constraint", "l" );
  lg_case_numu_full->AddEntry(h1_case_numu_constraint_full, "with #nu_{#mu}", "l" );
  lg_case_numu_full->Draw();

  canv_case_numu_full->SaveAs("canv_case_numu_full.png");
  
  //////////////////////////////////////
  roostr = "h1_case_numu_original_syst";
  TH1D *h1_case_numu_original_syst = (TH1D*)h1_nue->Clone(roostr);
  h1_case_numu_original_syst->Reset();
  roostr = "h1_case_numu_constraint_syst";
  TH1D *h1_case_numu_constraint_syst = (TH1D*)h1_nue->Clone(roostr);
  h1_case_numu_constraint_syst->Reset();
  for(int ibin=1; ibin<=case_numu_Ny; ibin++) {
    h1_case_numu_original_syst->SetBinContent( ibin, sqrt(array_err_syst[ibin-1])/array_val_mean[ibin-1] );
    h1_case_numu_constraint_syst->SetBinContent( ibin, sqrt(array_case_numu_constraint_syst[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_numu_syst";
  TCanvas *canv_case_numu_syst = new TCanvas(roostr, roostr, 900, 650);
  canv_case_numu_syst->SetRightMargin(0.05);
  canv_case_numu_syst->SetTopMargin(0.09);
  canv_case_numu_syst->SetBottomMargin(0.18);
  canv_case_numu_syst->SetLeftMargin(0.12);

  h1_case_numu_original_syst->Draw("hist");
  double case_numu_original_syst_max = 0;
  for(int ibin=1; ibin<=case_numu_Ny; ibin++) {
    double content = h1_case_numu_original_syst->GetBinContent(ibin);
    if( content>case_numu_original_syst_max ) case_numu_original_syst_max = content;
  }
  h1_case_numu_original_syst->SetMaximum( case_numu_original_syst_max*1.6 );
  h1_case_numu_original_syst->SetMinimum(0);
  h1_case_numu_original_syst->SetLineColor(kRed);
  h1_case_numu_original_syst->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_numu_original_syst->GetXaxis()->CenterTitle();
  h1_case_numu_original_syst->GetYaxis()->CenterTitle();
  h1_case_numu_original_syst->GetYaxis()->SetNdivisions(508);
  h1_case_numu_original_syst->SetYTitle("Relative uncertainty");
  h1_case_numu_original_syst->GetXaxis()->SetLabelSize(0.07);
  h1_case_numu_original_syst->GetXaxis()->SetTitleSize(0.07);
  h1_case_numu_original_syst->GetYaxis()->SetLabelSize(0.07);
  h1_case_numu_original_syst->GetYaxis()->SetTitleSize(0.07);  
  h1_case_numu_original_syst->GetXaxis()->SetTitleOffset(1.18);
  h1_case_numu_original_syst->GetYaxis()->SetTitleOffset(0.85);

  h1_case_numu_constraint_syst->Draw("same");
  h1_case_numu_constraint_syst->SetLineColor(kBlue);
  h1_case_numu_constraint_syst->Draw("same axis");

  TLegend *lg_case_numu_syst = new TLegend(0.154+0.3, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_numu_syst->SetBorderSize(1);
  lg_case_numu_syst->SetTextFont(42);
  lg_case_numu_syst->SetTextSize(0.068);
  lg_case_numu_syst->AddEntry(h1_case_numu_original_syst,   "no constraint", "l" );
  lg_case_numu_syst->AddEntry(h1_case_numu_constraint_syst, "with #nu_{#mu}", "l" );
  lg_case_numu_syst->Draw();

  canv_case_numu_syst->SaveAs("canv_case_numu_syst.png");
  
  ///////////////////////////////////////////////////////////////////////////////// nue low constrained by nue hgh
  ////////////////////////////////// case_nueH ////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  
  const int case_nueH_Nx = 6;
  const int case_nueH_Ny = 5;

  double array_case_nueH_constraint_full[case_nueH_Ny] = {0};
  double array_case_nueH_constraint_syst[case_nueH_Ny] = {0};
  
  TMatrixD case_nueH_convert_matrix( Ntotal, case_nueH_Nx+case_nueH_Ny );
  TMatrixD case_nueH_convert_matrixT( case_nueH_Nx+case_nueH_Ny, Ntotal );
  
  for(int i=0; i<case_nueH_Nx; i++) {// for x
    case_nueH_convert_matrix[i+case_nueH_Ny][i+0] = 1;// [old-index][new-index]
  }

  for(int i=0; i<case_nueH_Ny; i++) {// for y
    case_nueH_convert_matrix[i][i+case_nueH_Nx] = 1;// [old-index][new-index]
  }

  ////////
  case_nueH_convert_matrixT.Transpose( case_nueH_convert_matrix );

  TMatrixD case_nueH_cov_matrix = case_nueH_convert_matrixT * matrix_full_covariance * case_nueH_convert_matrix;

  //////
  roostr = "h2_check_case_nueH_cov_matrix";
  TH2D *h2_check_case_nueH_cov_matrix = new TH2D(roostr, roostr,
					     case_nueH_Nx+case_nueH_Ny, 0.5, case_nueH_Nx+case_nueH_Ny+0.5,
					     case_nueH_Nx+case_nueH_Ny, 0.5, case_nueH_Nx+case_nueH_Ny+0.5);
  
  for(int i=0; i<case_nueH_Nx+case_nueH_Ny; i++) {
    for(int j=0; j<case_nueH_Nx+case_nueH_Ny; j++) {
      double value = case_nueH_cov_matrix(i,j);
      h2_check_case_nueH_cov_matrix->SetBinContent(i+1, j+1, value );
    }
  }

  roostr = "canv_h2_check_case_nueH_cov_matrix";
  TCanvas *canv_h2_check_case_nueH_cov_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_case_nueH_cov_matrix->SetLeftMargin(0.12);
  canv_h2_check_case_nueH_cov_matrix->SetRightMargin(0.14);
  canv_h2_check_case_nueH_cov_matrix->SetTopMargin(0.09);
  canv_h2_check_case_nueH_cov_matrix->SetBottomMargin(0.18);
  canv_h2_check_case_nueH_cov_matrix->SetLogz();
  h2_check_case_nueH_cov_matrix->GetZaxis()->SetRangeUser(1.1e-2, 7e4);
  h2_check_case_nueH_cov_matrix->Draw("colz");
  h2_check_case_nueH_cov_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_case_nueH_cov_matrix->GetYaxis()->SetNdivisions(310);
  
  ////////
  TMatrixD case_nueH_cov_xx(case_nueH_Nx, case_nueH_Nx);
  TMatrixD case_nueH_cov_yy(case_nueH_Ny, case_nueH_Ny);
  TMatrixD case_nueH_cov_xy(case_nueH_Nx, case_nueH_Ny);
  TMatrixD case_nueH_cov_yx(case_nueH_Ny, case_nueH_Nx);
  TMatrixD case_nueH_cov_xx_inv(case_nueH_Nx, case_nueH_Nx);
  TMatrixD case_nueH_cov_y_under_x(case_nueH_Ny, case_nueH_Ny);

  for(int i=0; i<case_nueH_Nx; i++) {
    for(int j=0; j<case_nueH_Nx; j++) {
      case_nueH_cov_xx[i][j] = case_nueH_cov_matrix[i][j];
    }
  }

  for(int i=0; i<case_nueH_Ny; i++) {
    for(int j=0; j<case_nueH_Ny; j++) {
      case_nueH_cov_yy[i][j] = case_nueH_cov_matrix[i+case_nueH_Nx][j+case_nueH_Nx];
    }
  }
  
  for(int i=0; i<case_nueH_Nx; i++) {
    for(int j=0; j<case_nueH_Ny; j++) {
      case_nueH_cov_xy[i][j] = case_nueH_cov_matrix[i][j+case_nueH_Nx];
    }
  }
  
  case_nueH_cov_yx.Transpose( case_nueH_cov_xy );
  case_nueH_cov_xx_inv = case_nueH_cov_xx;
  case_nueH_cov_xx_inv.Invert();
  case_nueH_cov_y_under_x = case_nueH_cov_yy - case_nueH_cov_yx*case_nueH_cov_xx_inv*case_nueH_cov_xy;
  
  cout<<endl<<" ---> nueH_constraint_full"<<endl;
  for(int ibin=1; ibin<=case_nueH_Ny; ibin++) {
    double error = case_nueH_cov_y_under_x(ibin-1, ibin-1);
    error = sqrt(error);
    double value = array_val_mean[ibin-1];
    cout<<TString::Format("%2d rel.error %5.3f(full) %5.3f(syst)",
			  ibin,
			  error/value,
			  sqrt( error*error - array_err_stat[ibin-1]  )/value
			  )<<endl;
    
    array_case_nueH_constraint_full[ibin-1] = error*error;
    array_case_nueH_constraint_syst[ibin-1] = error*error - array_err_stat[ibin-1];
  }
  cout<<endl;

  //////////////////////////////////////
  roostr = "h1_case_nueH_original_full";
  TH1D *h1_case_nueH_original_full = new TH1D(roostr, roostr, case_nueH_Ny, 100, 600);
  roostr = "h1_case_nueH_constraint_full";
  TH1D *h1_case_nueH_constraint_full = new TH1D(roostr, roostr, case_nueH_Ny, 100, 600);
  for(int ibin=1; ibin<=case_nueH_Ny; ibin++) {
    h1_case_nueH_original_full->SetBinContent( ibin, sqrt(array_err_full[ibin-1])/array_val_mean[ibin-1] );
    h1_case_nueH_constraint_full->SetBinContent( ibin, sqrt(array_case_nueH_constraint_full[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_nueH_full";
  TCanvas *canv_case_nueH_full = new TCanvas(roostr, roostr, 900, 650);
  canv_case_nueH_full->SetRightMargin(0.05);
  canv_case_nueH_full->SetTopMargin(0.09);
  canv_case_nueH_full->SetBottomMargin(0.18);
  canv_case_nueH_full->SetLeftMargin(0.12);

  h1_case_nueH_original_full->Draw("hist");
  double case_nueH_original_full_max = 0;
  for(int ibin=1; ibin<=case_nueH_Ny; ibin++) {
    double content = h1_case_nueH_original_full->GetBinContent(ibin);
    if( content>case_nueH_original_full_max ) case_nueH_original_full_max = content;
  }
  h1_case_nueH_original_full->SetMaximum( case_nueH_original_full_max*1.2 );
  h1_case_nueH_original_full->SetMinimum(0);
  h1_case_nueH_original_full->SetLineColor(kRed);
  h1_case_nueH_original_full->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_nueH_original_full->GetXaxis()->CenterTitle();
  h1_case_nueH_original_full->GetYaxis()->CenterTitle();
  h1_case_nueH_original_full->GetYaxis()->SetNdivisions(508);
  h1_case_nueH_original_full->SetYTitle("Relative uncertainty");
  h1_case_nueH_original_full->GetXaxis()->SetLabelSize(0.07);
  h1_case_nueH_original_full->GetXaxis()->SetTitleSize(0.07);
  h1_case_nueH_original_full->GetYaxis()->SetLabelSize(0.07);
  h1_case_nueH_original_full->GetYaxis()->SetTitleSize(0.07);  
  h1_case_nueH_original_full->GetXaxis()->SetTitleOffset(1.18);
  h1_case_nueH_original_full->GetYaxis()->SetTitleOffset(0.85);

  h1_case_nueH_constraint_full->Draw("same");
  h1_case_nueH_constraint_full->SetLineColor(kBlue);
  h1_case_nueH_constraint_full->Draw("same axis");

  TLegend *lg_case_nueH_full = new TLegend(0.154+0.3, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_nueH_full->SetBorderSize(1);
  lg_case_nueH_full->SetTextFont(42);
  lg_case_nueH_full->SetTextSize(0.068);
  lg_case_nueH_full->AddEntry(h1_case_nueH_original_full,   "no constraint", "l" );
  lg_case_nueH_full->AddEntry(h1_case_nueH_constraint_full, "with #nu_{e} high E", "l" );
  lg_case_nueH_full->Draw();
  
  canv_case_nueH_full->SaveAs("canv_case_nueH_full.png");
  
  //////////////////////////////////////
  roostr = "h1_case_nueH_original_syst";
  TH1D *h1_case_nueH_original_syst = new TH1D(roostr, roostr, case_nueH_Ny, 100, 600);
  roostr = "h1_case_nueH_constraint_syst";
  TH1D *h1_case_nueH_constraint_syst = new TH1D(roostr, roostr, case_nueH_Ny, 100, 600);
  for(int ibin=1; ibin<=case_nueH_Ny; ibin++) {
    h1_case_nueH_original_syst->SetBinContent( ibin, sqrt(array_err_syst[ibin-1])/array_val_mean[ibin-1] );
    h1_case_nueH_constraint_syst->SetBinContent( ibin, sqrt(array_case_nueH_constraint_syst[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_nueH_syst";
  TCanvas *canv_case_nueH_syst = new TCanvas(roostr, roostr, 900, 650);
  canv_case_nueH_syst->SetRightMargin(0.05);
  canv_case_nueH_syst->SetTopMargin(0.09);
  canv_case_nueH_syst->SetBottomMargin(0.18);
  canv_case_nueH_syst->SetLeftMargin(0.12);

  h1_case_nueH_original_syst->Draw("hist");
  double case_nueH_original_syst_max = 0;
  for(int ibin=1; ibin<=case_nueH_Ny; ibin++) {
    double content = h1_case_nueH_original_syst->GetBinContent(ibin);
    if( content>case_nueH_original_syst_max ) case_nueH_original_syst_max = content;
  }
  h1_case_nueH_original_syst->SetMaximum( case_nueH_original_syst_max*1.6 );
  h1_case_nueH_original_syst->SetMinimum(0);
  h1_case_nueH_original_syst->SetLineColor(kRed);
  h1_case_nueH_original_syst->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_nueH_original_syst->GetXaxis()->CenterTitle();
  h1_case_nueH_original_syst->GetYaxis()->CenterTitle();
  h1_case_nueH_original_syst->GetYaxis()->SetNdivisions(508);
  h1_case_nueH_original_syst->SetYTitle("Relative uncertainty");
  h1_case_nueH_original_syst->GetXaxis()->SetLabelSize(0.07);
  h1_case_nueH_original_syst->GetXaxis()->SetTitleSize(0.07);
  h1_case_nueH_original_syst->GetYaxis()->SetLabelSize(0.07);
  h1_case_nueH_original_syst->GetYaxis()->SetTitleSize(0.07);  
  h1_case_nueH_original_syst->GetXaxis()->SetTitleOffset(1.18);
  h1_case_nueH_original_syst->GetYaxis()->SetTitleOffset(0.85);

  h1_case_nueH_constraint_syst->Draw("same");
  h1_case_nueH_constraint_syst->SetLineColor(kBlue);
  h1_case_nueH_constraint_syst->Draw("same axis");

  TLegend *lg_case_nueH_syst = new TLegend(0.154+0.3, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_nueH_syst->SetBorderSize(1);
  lg_case_nueH_syst->SetTextFont(42);
  lg_case_nueH_syst->SetTextSize(0.068);
  lg_case_nueH_syst->AddEntry(h1_case_nueH_original_syst,   "no constraint", "l" );
  lg_case_nueH_syst->AddEntry(h1_case_nueH_constraint_syst, "with #nu_{e} high E", "l" );
  lg_case_nueH_syst->Draw();
  
  canv_case_nueH_syst->SaveAs("canv_case_nueH_syst.png");
  
  ///////////////////////////////////////////////////////////////////////////////// nue low constrained by nue hgh, numu low
  ////////////////////////////////// case_nueH_numuL ////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  
  const int case_nueH_numuL_Nx = 6+5;
  const int case_nueH_numuL_Ny = 5;

  double array_case_nueH_numuL_constraint_full[case_nueH_numuL_Ny] = {0};
  double array_case_nueH_numuL_constraint_syst[case_nueH_numuL_Ny] = {0};
  
  TMatrixD case_nueH_numuL_convert_matrix( Ntotal, case_nueH_numuL_Nx+case_nueH_numuL_Ny );
  TMatrixD case_nueH_numuL_convert_matrixT( case_nueH_numuL_Nx+case_nueH_numuL_Ny, Ntotal );
  
  for(int i=0; i<6; i++) {// for x
    case_nueH_numuL_convert_matrix[i+5][i+0] = 1;// [old-index][new-index]
  }

  for(int i=0; i<5; i++) {// for x
    case_nueH_numuL_convert_matrix[i+11][i+6] = 1;// [old-index][new-index]
  }

  for(int i=0; i<case_nueH_numuL_Ny; i++) {// for y
    case_nueH_numuL_convert_matrix[i][i+case_nueH_numuL_Nx] = 1;// [old-index][new-index]
  }

  ////////
  case_nueH_numuL_convert_matrixT.Transpose( case_nueH_numuL_convert_matrix );

  TMatrixD case_nueH_numuL_cov_matrix = case_nueH_numuL_convert_matrixT * matrix_full_covariance * case_nueH_numuL_convert_matrix;

  //////
  roostr = "h2_check_case_nueH_numuL_cov_matrix";
  TH2D *h2_check_case_nueH_numuL_cov_matrix = new TH2D(roostr, roostr,
					     case_nueH_numuL_Nx+case_nueH_numuL_Ny, 0.5, case_nueH_numuL_Nx+case_nueH_numuL_Ny+0.5,
					     case_nueH_numuL_Nx+case_nueH_numuL_Ny, 0.5, case_nueH_numuL_Nx+case_nueH_numuL_Ny+0.5);
  
  for(int i=0; i<case_nueH_numuL_Nx+case_nueH_numuL_Ny; i++) {
    for(int j=0; j<case_nueH_numuL_Nx+case_nueH_numuL_Ny; j++) {
      double value = case_nueH_numuL_cov_matrix(i,j);
      h2_check_case_nueH_numuL_cov_matrix->SetBinContent(i+1, j+1, value );
    }
  }

  roostr = "canv_h2_check_case_nueH_numuL_cov_matrix";
  TCanvas *canv_h2_check_case_nueH_numuL_cov_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_case_nueH_numuL_cov_matrix->SetLeftMargin(0.12);
  canv_h2_check_case_nueH_numuL_cov_matrix->SetRightMargin(0.14);
  canv_h2_check_case_nueH_numuL_cov_matrix->SetTopMargin(0.09);
  canv_h2_check_case_nueH_numuL_cov_matrix->SetBottomMargin(0.18);
  canv_h2_check_case_nueH_numuL_cov_matrix->SetLogz();
  h2_check_case_nueH_numuL_cov_matrix->GetZaxis()->SetRangeUser(1.1e-2, 7e4);
  h2_check_case_nueH_numuL_cov_matrix->Draw("colz");
  h2_check_case_nueH_numuL_cov_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_case_nueH_numuL_cov_matrix->GetYaxis()->SetNdivisions(310);
  
  ////////
  TMatrixD case_nueH_numuL_cov_xx(case_nueH_numuL_Nx, case_nueH_numuL_Nx);
  TMatrixD case_nueH_numuL_cov_yy(case_nueH_numuL_Ny, case_nueH_numuL_Ny);
  TMatrixD case_nueH_numuL_cov_xy(case_nueH_numuL_Nx, case_nueH_numuL_Ny);
  TMatrixD case_nueH_numuL_cov_yx(case_nueH_numuL_Ny, case_nueH_numuL_Nx);
  TMatrixD case_nueH_numuL_cov_xx_inv(case_nueH_numuL_Nx, case_nueH_numuL_Nx);
  TMatrixD case_nueH_numuL_cov_y_under_x(case_nueH_numuL_Ny, case_nueH_numuL_Ny);

  for(int i=0; i<case_nueH_numuL_Nx; i++) {
    for(int j=0; j<case_nueH_numuL_Nx; j++) {
      case_nueH_numuL_cov_xx[i][j] = case_nueH_numuL_cov_matrix[i][j];
    }
  }

  for(int i=0; i<case_nueH_numuL_Ny; i++) {
    for(int j=0; j<case_nueH_numuL_Ny; j++) {
      case_nueH_numuL_cov_yy[i][j] = case_nueH_numuL_cov_matrix[i+case_nueH_numuL_Nx][j+case_nueH_numuL_Nx];
    }
  }
  
  for(int i=0; i<case_nueH_numuL_Nx; i++) {
    for(int j=0; j<case_nueH_numuL_Ny; j++) {
      case_nueH_numuL_cov_xy[i][j] = case_nueH_numuL_cov_matrix[i][j+case_nueH_numuL_Nx];
    }
  }
  
  case_nueH_numuL_cov_yx.Transpose( case_nueH_numuL_cov_xy );
  case_nueH_numuL_cov_xx_inv = case_nueH_numuL_cov_xx;
  case_nueH_numuL_cov_xx_inv.Invert();
  case_nueH_numuL_cov_y_under_x = case_nueH_numuL_cov_yy - case_nueH_numuL_cov_yx*case_nueH_numuL_cov_xx_inv*case_nueH_numuL_cov_xy;
  
  cout<<endl<<" ---> nueH_numuL_constraint_full"<<endl;
  for(int ibin=1; ibin<=case_nueH_numuL_Ny; ibin++) {
    double error = case_nueH_numuL_cov_y_under_x(ibin-1, ibin-1);
    error = sqrt(error);
    double value = array_val_mean[ibin-1];
    cout<<TString::Format("%2d rel.error %5.3f(full) %5.3f(syst)",
			  ibin,
			  error/value,
			  sqrt( error*error - array_err_stat[ibin-1]  )/value
			  )<<endl;
    
    array_case_nueH_numuL_constraint_full[ibin-1] = error*error;
    array_case_nueH_numuL_constraint_syst[ibin-1] = error*error - array_err_stat[ibin-1];
  }
  cout<<endl;

  ////////////////////////////////////////////
  roostr = "h1_case_nueH_numuL_original_full";
  TH1D *h1_case_nueH_numuL_original_full = new TH1D(roostr, roostr, case_nueH_numuL_Ny, 100, 600);
  roostr = "h1_case_nueH_numuL_constraint_full";
  TH1D *h1_case_nueH_numuL_constraint_full = new TH1D(roostr, roostr, case_nueH_numuL_Ny, 100, 600);
  for(int ibin=1; ibin<=case_nueH_numuL_Ny; ibin++) {
    h1_case_nueH_numuL_original_full->SetBinContent( ibin, sqrt(array_err_full[ibin-1])/array_val_mean[ibin-1] );
    h1_case_nueH_numuL_constraint_full->SetBinContent( ibin, sqrt(array_case_nueH_numuL_constraint_full[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_nueH_numuL_full";
  TCanvas *canv_case_nueH_numuL_full = new TCanvas(roostr, roostr, 900, 650);
  canv_case_nueH_numuL_full->SetRightMargin(0.05);
  canv_case_nueH_numuL_full->SetTopMargin(0.09);
  canv_case_nueH_numuL_full->SetBottomMargin(0.18);
  canv_case_nueH_numuL_full->SetLeftMargin(0.12);

  h1_case_nueH_numuL_original_full->Draw("hist");
  double case_nueH_numuL_original_full_max = 0;
  for(int ibin=1; ibin<=case_nueH_numuL_Ny; ibin++) {
    double content = h1_case_nueH_numuL_original_full->GetBinContent(ibin);
    if( content>case_nueH_numuL_original_full_max ) case_nueH_numuL_original_full_max = content;
  }
  h1_case_nueH_numuL_original_full->SetMaximum( case_nueH_numuL_original_full_max*1.2 );
  h1_case_nueH_numuL_original_full->SetMinimum(0);
  h1_case_nueH_numuL_original_full->SetLineColor(kRed);
  h1_case_nueH_numuL_original_full->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_nueH_numuL_original_full->GetXaxis()->CenterTitle();
  h1_case_nueH_numuL_original_full->GetYaxis()->CenterTitle();
  h1_case_nueH_numuL_original_full->GetYaxis()->SetNdivisions(508);
  h1_case_nueH_numuL_original_full->SetYTitle("Relative uncertainty");
  h1_case_nueH_numuL_original_full->GetXaxis()->SetLabelSize(0.07);
  h1_case_nueH_numuL_original_full->GetXaxis()->SetTitleSize(0.07);
  h1_case_nueH_numuL_original_full->GetYaxis()->SetLabelSize(0.07);
  h1_case_nueH_numuL_original_full->GetYaxis()->SetTitleSize(0.07);  
  h1_case_nueH_numuL_original_full->GetXaxis()->SetTitleOffset(1.18);
  h1_case_nueH_numuL_original_full->GetYaxis()->SetTitleOffset(0.85);

  h1_case_nueH_numuL_constraint_full->Draw("same");
  h1_case_nueH_numuL_constraint_full->SetLineColor(kBlue);
  h1_case_nueH_numuL_constraint_full->Draw("same axis");

  TLegend *lg_case_nueH_numuL_full = new TLegend(0.154+0.3-0.15, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_nueH_numuL_full->SetBorderSize(1);
  lg_case_nueH_numuL_full->SetTextFont(42);
  lg_case_nueH_numuL_full->SetTextSize(0.068);
  lg_case_nueH_numuL_full->AddEntry(h1_case_nueH_numuL_original_full,   "no constraint", "l" );
  lg_case_nueH_numuL_full->AddEntry(h1_case_nueH_numuL_constraint_full, "with #nu_{e} high E, #nu_{#mu} low", "l" );
  lg_case_nueH_numuL_full->Draw();
  
  canv_case_nueH_numuL_full->SaveAs("canv_case_nueH_numuL_full.png");
    
  ////////////////////////////////////////////
  roostr = "h1_case_nueH_numuL_original_syst";
  TH1D *h1_case_nueH_numuL_original_syst = new TH1D(roostr, roostr, case_nueH_numuL_Ny, 100, 600);
  roostr = "h1_case_nueH_numuL_constraint_syst";
  TH1D *h1_case_nueH_numuL_constraint_syst = new TH1D(roostr, roostr, case_nueH_numuL_Ny, 100, 600);
  for(int ibin=1; ibin<=case_nueH_numuL_Ny; ibin++) {
    h1_case_nueH_numuL_original_syst->SetBinContent( ibin, sqrt(array_err_syst[ibin-1])/array_val_mean[ibin-1] );
    h1_case_nueH_numuL_constraint_syst->SetBinContent( ibin, sqrt(array_case_nueH_numuL_constraint_syst[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_nueH_numuL_syst";
  TCanvas *canv_case_nueH_numuL_syst = new TCanvas(roostr, roostr, 900, 650);
  canv_case_nueH_numuL_syst->SetRightMargin(0.05);
  canv_case_nueH_numuL_syst->SetTopMargin(0.09);
  canv_case_nueH_numuL_syst->SetBottomMargin(0.18);
  canv_case_nueH_numuL_syst->SetLeftMargin(0.12);

  h1_case_nueH_numuL_original_syst->Draw("hist");
  double case_nueH_numuL_original_syst_max = 0;
  for(int ibin=1; ibin<=case_nueH_numuL_Ny; ibin++) {
    double content = h1_case_nueH_numuL_original_syst->GetBinContent(ibin);
    if( content>case_nueH_numuL_original_syst_max ) case_nueH_numuL_original_syst_max = content;
  }
  h1_case_nueH_numuL_original_syst->SetMaximum( case_nueH_numuL_original_syst_max*1.6 );
  h1_case_nueH_numuL_original_syst->SetMinimum(0);
  h1_case_nueH_numuL_original_syst->SetLineColor(kRed);
  h1_case_nueH_numuL_original_syst->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_nueH_numuL_original_syst->GetXaxis()->CenterTitle();
  h1_case_nueH_numuL_original_syst->GetYaxis()->CenterTitle();
  h1_case_nueH_numuL_original_syst->GetYaxis()->SetNdivisions(508);
  h1_case_nueH_numuL_original_syst->SetYTitle("Relative uncertainty");
  h1_case_nueH_numuL_original_syst->GetXaxis()->SetLabelSize(0.07);
  h1_case_nueH_numuL_original_syst->GetXaxis()->SetTitleSize(0.07);
  h1_case_nueH_numuL_original_syst->GetYaxis()->SetLabelSize(0.07);
  h1_case_nueH_numuL_original_syst->GetYaxis()->SetTitleSize(0.07);  
  h1_case_nueH_numuL_original_syst->GetXaxis()->SetTitleOffset(1.18);
  h1_case_nueH_numuL_original_syst->GetYaxis()->SetTitleOffset(0.85);

  h1_case_nueH_numuL_constraint_syst->Draw("same");
  h1_case_nueH_numuL_constraint_syst->SetLineColor(kBlue);
  h1_case_nueH_numuL_constraint_syst->Draw("same axis");

  TLegend *lg_case_nueH_numuL_syst = new TLegend(0.154+0.3-0.15, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_nueH_numuL_syst->SetBorderSize(1);
  lg_case_nueH_numuL_syst->SetTextFont(42);
  lg_case_nueH_numuL_syst->SetTextSize(0.068);
  lg_case_nueH_numuL_syst->AddEntry(h1_case_nueH_numuL_original_syst,   "no constraint", "l" );
  lg_case_nueH_numuL_syst->AddEntry(h1_case_nueH_numuL_constraint_syst, "with #nu_{e} high E, #nu_{#mu} low", "l" );
  lg_case_nueH_numuL_syst->Draw();
  
  canv_case_nueH_numuL_syst->SaveAs("canv_case_nueH_numuL_syst.png");
    
  ///////////////////////////////////////////////////////////////////////////////// nue low constrained by numu low
  ////////////////////////////////// case_numuL ////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  
  const int case_numuL_Nx = 5;
  const int case_numuL_Ny = 5;

  double array_case_numuL_constraint_full[case_numuL_Ny] = {0};
  double array_case_numuL_constraint_syst[case_numuL_Ny] = {0};
  
  TMatrixD case_numuL_convert_matrix( Ntotal, case_numuL_Nx+case_numuL_Ny );
  TMatrixD case_numuL_convert_matrixT( case_numuL_Nx+case_numuL_Ny, Ntotal );
  
  for(int i=0; i<case_numuL_Nx; i++) {// for x
    case_numuL_convert_matrix[i+11][i+0] = 1;// [old-index][new-index]
  }

  for(int i=0; i<case_numuL_Ny; i++) {// for y
    case_numuL_convert_matrix[i][i+case_numuL_Nx] = 1;// [old-index][new-index]
  }

  ////////
  case_numuL_convert_matrixT.Transpose( case_numuL_convert_matrix );

  TMatrixD case_numuL_cov_matrix = case_numuL_convert_matrixT * matrix_full_covariance * case_numuL_convert_matrix;

  //////
  roostr = "h2_check_case_numuL_cov_matrix";
  TH2D *h2_check_case_numuL_cov_matrix = new TH2D(roostr, roostr,
					     case_numuL_Nx+case_numuL_Ny, 0.5, case_numuL_Nx+case_numuL_Ny+0.5,
					     case_numuL_Nx+case_numuL_Ny, 0.5, case_numuL_Nx+case_numuL_Ny+0.5);
  
  for(int i=0; i<case_numuL_Nx+case_numuL_Ny; i++) {
    for(int j=0; j<case_numuL_Nx+case_numuL_Ny; j++) {
      double value = case_numuL_cov_matrix(i,j);
      h2_check_case_numuL_cov_matrix->SetBinContent(i+1, j+1, value );
    }
  }

  roostr = "canv_h2_check_case_numuL_cov_matrix";
  TCanvas *canv_h2_check_case_numuL_cov_matrix = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_check_case_numuL_cov_matrix->SetLeftMargin(0.12);
  canv_h2_check_case_numuL_cov_matrix->SetRightMargin(0.14);
  canv_h2_check_case_numuL_cov_matrix->SetTopMargin(0.09);
  canv_h2_check_case_numuL_cov_matrix->SetBottomMargin(0.18);
  canv_h2_check_case_numuL_cov_matrix->SetLogz();
  h2_check_case_numuL_cov_matrix->GetZaxis()->SetRangeUser(1.1e-2, 7e4);
  h2_check_case_numuL_cov_matrix->Draw("colz");
  h2_check_case_numuL_cov_matrix->GetXaxis()->SetNdivisions(310);
  h2_check_case_numuL_cov_matrix->GetYaxis()->SetNdivisions(310);
  
  ////////
  TMatrixD case_numuL_cov_xx(case_numuL_Nx, case_numuL_Nx);
  TMatrixD case_numuL_cov_yy(case_numuL_Ny, case_numuL_Ny);
  TMatrixD case_numuL_cov_xy(case_numuL_Nx, case_numuL_Ny);
  TMatrixD case_numuL_cov_yx(case_numuL_Ny, case_numuL_Nx);
  TMatrixD case_numuL_cov_xx_inv(case_numuL_Nx, case_numuL_Nx);
  TMatrixD case_numuL_cov_y_under_x(case_numuL_Ny, case_numuL_Ny);

  for(int i=0; i<case_numuL_Nx; i++) {
    for(int j=0; j<case_numuL_Nx; j++) {
      case_numuL_cov_xx[i][j] = case_numuL_cov_matrix[i][j];
    }
  }

  for(int i=0; i<case_numuL_Ny; i++) {
    for(int j=0; j<case_numuL_Ny; j++) {
      case_numuL_cov_yy[i][j] = case_numuL_cov_matrix[i+case_numuL_Nx][j+case_numuL_Nx];
    }
  }
  
  for(int i=0; i<case_numuL_Nx; i++) {
    for(int j=0; j<case_numuL_Ny; j++) {
      case_numuL_cov_xy[i][j] = case_numuL_cov_matrix[i][j+case_numuL_Nx];
    }
  }
  
  case_numuL_cov_yx.Transpose( case_numuL_cov_xy );
  case_numuL_cov_xx_inv = case_numuL_cov_xx;
  case_numuL_cov_xx_inv.Invert();
  case_numuL_cov_y_under_x = case_numuL_cov_yy - case_numuL_cov_yx*case_numuL_cov_xx_inv*case_numuL_cov_xy;
  
  cout<<endl<<" ---> numuL_constraint_full"<<endl;
  for(int ibin=1; ibin<=case_numuL_Ny; ibin++) {
    double error = case_numuL_cov_y_under_x(ibin-1, ibin-1);
    error = sqrt(error);
    double value = array_val_mean[ibin-1];
    cout<<TString::Format("%2d rel.error %5.3f(full) %5.3f(syst)",
			  ibin,
			  error/value,
			  sqrt( error*error - array_err_stat[ibin-1]  )/value
			  )<<endl;
    
    array_case_numuL_constraint_full[ibin-1] = error*error;
    array_case_numuL_constraint_syst[ibin-1] = error*error - array_err_stat[ibin-1];
  }
  cout<<endl;

  //////////////////////////////////////
  roostr = "h1_case_numuL_original_full";
  TH1D *h1_case_numuL_original_full = new TH1D(roostr, roostr, case_numuL_Ny, 100, 600);
  roostr = "h1_case_numuL_constraint_full";
  TH1D *h1_case_numuL_constraint_full = new TH1D(roostr, roostr, case_numuL_Ny, 100, 600);
  for(int ibin=1; ibin<=case_numuL_Ny; ibin++) {
    h1_case_numuL_original_full->SetBinContent( ibin, sqrt(array_err_full[ibin-1])/array_val_mean[ibin-1] );
    h1_case_numuL_constraint_full->SetBinContent( ibin, sqrt(array_case_numuL_constraint_full[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_numuL_full";
  TCanvas *canv_case_numuL_full = new TCanvas(roostr, roostr, 900, 650);
  canv_case_numuL_full->SetRightMargin(0.05);
  canv_case_numuL_full->SetTopMargin(0.09);
  canv_case_numuL_full->SetBottomMargin(0.18);
  canv_case_numuL_full->SetLeftMargin(0.12);

  h1_case_numuL_original_full->Draw("hist");
  double case_numuL_original_full_max = 0;
  for(int ibin=1; ibin<=case_numuL_Ny; ibin++) {
    double content = h1_case_numuL_original_full->GetBinContent(ibin);
    if( content>case_numuL_original_full_max ) case_numuL_original_full_max = content;
  }
  h1_case_numuL_original_full->SetMaximum( case_numuL_original_full_max*1.2 );
  h1_case_numuL_original_full->SetMinimum(0);
  h1_case_numuL_original_full->SetLineColor(kRed);
  h1_case_numuL_original_full->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_numuL_original_full->GetXaxis()->CenterTitle();
  h1_case_numuL_original_full->GetYaxis()->CenterTitle();
  h1_case_numuL_original_full->GetYaxis()->SetNdivisions(508);
  h1_case_numuL_original_full->SetYTitle("Relative uncertainty");
  h1_case_numuL_original_full->GetXaxis()->SetLabelSize(0.07);
  h1_case_numuL_original_full->GetXaxis()->SetTitleSize(0.07);
  h1_case_numuL_original_full->GetYaxis()->SetLabelSize(0.07);
  h1_case_numuL_original_full->GetYaxis()->SetTitleSize(0.07);  
  h1_case_numuL_original_full->GetXaxis()->SetTitleOffset(1.18);
  h1_case_numuL_original_full->GetYaxis()->SetTitleOffset(0.85);

  h1_case_numuL_constraint_full->Draw("same");
  h1_case_numuL_constraint_full->SetLineColor(kBlue);
  h1_case_numuL_constraint_full->Draw("same axis");

  TLegend *lg_case_numuL_full = new TLegend(0.154+0.3-0.15, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_numuL_full->SetBorderSize(1);
  lg_case_numuL_full->SetTextFont(42);
  lg_case_numuL_full->SetTextSize(0.068);
  lg_case_numuL_full->AddEntry(h1_case_numuL_original_full,   "no constraint", "l" );
  lg_case_numuL_full->AddEntry(h1_case_numuL_constraint_full, "with #nu_{#mu} low E", "l" );
  lg_case_numuL_full->Draw();
  
  canv_case_numuL_full->SaveAs("canv_case_numuL_full.png");
  
  //////////////////////////////////////
  roostr = "h1_case_numuL_original_syst";
  TH1D *h1_case_numuL_original_syst = new TH1D(roostr, roostr, case_numuL_Ny, 100, 600);
  roostr = "h1_case_numuL_constraint_syst";
  TH1D *h1_case_numuL_constraint_syst = new TH1D(roostr, roostr, case_numuL_Ny, 100, 600);
  for(int ibin=1; ibin<=case_numuL_Ny; ibin++) {
    h1_case_numuL_original_syst->SetBinContent( ibin, sqrt(array_err_syst[ibin-1])/array_val_mean[ibin-1] );
    h1_case_numuL_constraint_syst->SetBinContent( ibin, sqrt(array_case_numuL_constraint_syst[ibin-1])/array_val_mean[ibin-1] );
  }

  roostr = "canv_case_numuL_syst";
  TCanvas *canv_case_numuL_syst = new TCanvas(roostr, roostr, 900, 650);
  canv_case_numuL_syst->SetRightMargin(0.05);
  canv_case_numuL_syst->SetTopMargin(0.09);
  canv_case_numuL_syst->SetBottomMargin(0.18);
  canv_case_numuL_syst->SetLeftMargin(0.12);

  h1_case_numuL_original_syst->Draw("hist");
  double case_numuL_original_syst_max = 0;
  for(int ibin=1; ibin<=case_numuL_Ny; ibin++) {
    double content = h1_case_numuL_original_syst->GetBinContent(ibin);
    if( content>case_numuL_original_syst_max ) case_numuL_original_syst_max = content;
  }
  h1_case_numuL_original_syst->SetMaximum( case_numuL_original_syst_max*1.6 );
  h1_case_numuL_original_syst->SetMinimum(0);
  h1_case_numuL_original_syst->SetLineColor(kRed);
  h1_case_numuL_original_syst->SetXTitle("Reconstructed Neutrino Energy [MeV]");
  h1_case_numuL_original_syst->GetXaxis()->CenterTitle();
  h1_case_numuL_original_syst->GetYaxis()->CenterTitle();
  h1_case_numuL_original_syst->GetYaxis()->SetNdivisions(508);
  h1_case_numuL_original_syst->SetYTitle("Relative uncertainty");
  h1_case_numuL_original_syst->GetXaxis()->SetLabelSize(0.07);
  h1_case_numuL_original_syst->GetXaxis()->SetTitleSize(0.07);
  h1_case_numuL_original_syst->GetYaxis()->SetLabelSize(0.07);
  h1_case_numuL_original_syst->GetYaxis()->SetTitleSize(0.07);  
  h1_case_numuL_original_syst->GetXaxis()->SetTitleOffset(1.18);
  h1_case_numuL_original_syst->GetYaxis()->SetTitleOffset(0.85);

  h1_case_numuL_constraint_syst->Draw("same");
  h1_case_numuL_constraint_syst->SetLineColor(kBlue);
  h1_case_numuL_constraint_syst->Draw("same axis");

  TLegend *lg_case_numuL_syst = new TLegend(0.154+0.3-0.15, 0.580+0.075+0.05, 0.6247+0.3, 0.881);
  lg_case_numuL_syst->SetBorderSize(1);
  lg_case_numuL_syst->SetTextFont(42);
  lg_case_numuL_syst->SetTextSize(0.068);
  lg_case_numuL_syst->AddEntry(h1_case_numuL_original_syst,   "no constraint", "l" );
  lg_case_numuL_syst->AddEntry(h1_case_numuL_constraint_syst, "with #nu_{#mu} low E", "l" );
  lg_case_numuL_syst->Draw();
  
  canv_case_numuL_syst->SaveAs("canv_case_numuL_syst.png");
  
  //////////////////////////
  //////////////////////////
  //////////////////////////

  ///////////
  cout<<endl;
  for(int ibin=1; ibin<=5; ibin++) {
    double AA = h1_case_numu_constraint_full->GetBinContent(ibin);// nue constrained by numu
    double BB = h1_case_nueH_constraint_full->GetBinContent(ibin);// nue low constrained by nue hgh
    double CC = h1_case_nueH_numuL_constraint_full->GetBinContent(ibin);// nue low constrained by nue hgh, numu low
    double DD = h1_case_numuL_constraint_full->GetBinContent(ibin);// nue low constrained by numu low

    cout<<TString::Format("%2d %5.3f(0) %5.3f(A) %5.3f(B) %5.3f(C) %5.3f(D)",
			  ibin,
			  sqrt(array_err_full[ibin-1])/array_val_mean[ibin-1],
			  AA, BB, CC, DD
			  )<<endl;
  }
  cout<<endl;

  /////////////////////////////
  roostr = "canv_casefull_full";
  TCanvas *canv_casefull_full = new TCanvas(roostr, roostr, 900, 650);
  canv_casefull_full->SetRightMargin(0.05);
  canv_casefull_full->SetTopMargin(0.09);
  canv_casefull_full->SetBottomMargin(0.18);
  canv_casefull_full->SetLeftMargin(0.12);
  
  h1_case_numuL_original_full->Draw();
  h1_case_numuL_original_full->SetLineColor( kRed );

  h1_case_nueH_constraint_full->Draw("same");
  h1_case_nueH_constraint_full->SetLineColor( kBlack );

  h1_case_numu_constraint_full->Draw("same");
  h1_case_numu_constraint_full->SetLineColor( kOrange );

  h1_case_nueH_numuL_constraint_full->Draw("same");
  h1_case_nueH_numuL_constraint_full->SetLineColor( kCyan );

  h1_case_numuL_constraint_full->Draw("same");
  h1_case_numuL_constraint_full->SetLineColor( kBlue );

  h1_case_numuL_original_full->Draw("same axis");

  TLegend *lg_casefull_full = new TLegend(0.154+0.3-0.15, 0.580+0.075+0.05-0.2, 0.6247+0.3, 0.881);
  lg_casefull_full->SetBorderSize(1);
  lg_casefull_full->SetTextFont(42);
  lg_casefull_full->SetTextSize(0.068);
  lg_casefull_full->AddEntry(h1_case_numuL_original_full,   "no constraint", "l" );
  lg_casefull_full->AddEntry(h1_case_nueH_constraint_full, "with #nu_{e} high E", "l" );
  lg_casefull_full->AddEntry(h1_case_numuL_constraint_full, "with #nu_{#mu} low E", "l" );
  lg_casefull_full->AddEntry(h1_case_nueH_numuL_constraint_full, "with #nu_{e} high, #nu_{#mu} low", "l" );
  lg_casefull_full->AddEntry(h1_case_numu_constraint_full, "with #nu_{#mu}", "l" );  
  lg_casefull_full->Draw();
  
  canv_casefull_full->SaveAs("canv_casefull_full.png");
  
  /////////////////////////////
  roostr = "canv_casesyst_syst";
  TCanvas *canv_casesyst_syst = new TCanvas(roostr, roostr, 900, 650);
  canv_casesyst_syst->SetRightMargin(0.05);
  canv_casesyst_syst->SetTopMargin(0.09);
  canv_casesyst_syst->SetBottomMargin(0.18);
  canv_casesyst_syst->SetLeftMargin(0.12);
  
  h1_case_numuL_original_syst->Draw();
  h1_case_numuL_original_syst->SetLineColor( kRed );
  h1_case_numuL_original_syst->SetMaximum( h1_case_numuL_original_syst->GetMaximum()*1.4 );

  h1_case_nueH_constraint_syst->Draw("same");
  h1_case_nueH_constraint_syst->SetLineColor( kBlack );

  h1_case_numu_constraint_syst->Draw("same");
  h1_case_numu_constraint_syst->SetLineColor( kOrange );

  h1_case_nueH_numuL_constraint_syst->Draw("same");
  h1_case_nueH_numuL_constraint_syst->SetLineColor( kCyan );

  h1_case_numuL_constraint_syst->Draw("same");
  h1_case_numuL_constraint_syst->SetLineColor( kBlue );

  h1_case_numuL_original_syst->Draw("same axis");

  TLegend *lg_casesyst_syst = new TLegend(0.154+0.3-0.15+0.1, 0.580+0.075+0.05-0.2+0.05, 0.6247+0.3, 0.881);
  lg_casesyst_syst->SetBorderSize(1);
  lg_casesyst_syst->SetTextFont(42);
  lg_casesyst_syst->SetTextSize(0.068);
  lg_casesyst_syst->AddEntry(h1_case_numuL_original_syst,   "no constraint", "l" );
  lg_casesyst_syst->AddEntry(h1_case_nueH_constraint_syst, "with #nu_{e} high E", "l" );
  lg_casesyst_syst->AddEntry(h1_case_numuL_constraint_syst, "with #nu_{#mu} low E", "l" );
  lg_casesyst_syst->AddEntry(h1_case_nueH_numuL_constraint_syst, "with #nu_{e} high, #nu_{#mu} low", "l" );
  lg_casesyst_syst->AddEntry(h1_case_numu_constraint_syst, "with #nu_{#mu}", "l" );  
  lg_casesyst_syst->Draw();
  
  canv_casesyst_syst->SaveAs("canv_casesyst_syst.png");
  
  
}
