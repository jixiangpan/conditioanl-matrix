
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// MAIN ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_fitter(int toy_i, int toy_j, int toy_flag)
{

  TString roostr = "";
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
    
  double snWidth = 2;
  
  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);

  gStyle->SetEndErrorSize(4);
  
  ///////////////////////////////////////////////////////// TEST
  ///////////////////////////////////////////////////////// TEST
  ///////////////////////////////////////////////////////// TEST
  /*
  const int N_test_total = 10;
  
  TMatrixD matrix_total_test(N_test_total, N_test_total);
  for(int i=0; i<N_test_total; i++) {
    for(int j=0; j<N_test_total; j++) {
      matrix_total_test(i, j) = i*10 + j;
    }
  }

  roostr = "canv_matrix_total_test";
  TCanvas *canv_matrix_total_test = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_matrix_total_test, 0.15, 0.2,0.1,0.15);
  matrix_total_test.Draw("text");

  /////////// subA
  /////////// sub_matrix = R^T * total_matrix * R
  
  const int N_subA = 5;

  TMatrixD matrix_trans_subA(N_test_total, N_subA);
  TMatrixD T_matrix_trans_subA(N_subA, N_test_total);
  for(int i=0; i<N_subA; i++) {
    matrix_trans_subA(i+0, i+0) = 1;// (old_world_index, new_world_index)
  }  
  T_matrix_trans_subA.Transpose( matrix_trans_subA );

  TMatrixD matrix_subA = T_matrix_trans_subA * matrix_total_test * matrix_trans_subA;

  roostr = "canv_matrix_subA";
  TCanvas *canv_matrix_subA = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_matrix_subA, 0.15, 0.2,0.1,0.15);
  matrix_subA.Draw("text");

  /////////// subB
  /////////// sub_matrix = R^T * total_matrix * R
  
  const int N_subB = 5;

  TMatrixD matrix_trans_subB(N_test_total, N_subB);
  TMatrixD T_matrix_trans_subB(N_subB, N_test_total);
  for(int i=0; i<2; i++) {
    matrix_trans_subB(i+2, i) = 1;// (old_world_index, new_world_index)
  }
  for(int i=2; i<N_subB; i++) {
    matrix_trans_subB(i+4, i) = 1;// (old_world_index, new_world_index)
  } 
  T_matrix_trans_subB.Transpose( matrix_trans_subB );

  TMatrixD matrix_subB = T_matrix_trans_subB * matrix_total_test * matrix_trans_subB;

  roostr = "canv_matrix_subB";
  TCanvas *canv_matrix_subB = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_matrix_subB, 0.15, 0.2,0.1,0.15);
  matrix_subB.Draw("text");

  /////////// subC
  /////////// sub_matrix = R^T * total_matrix * R
  
  const int N_subC = 1;

  TMatrixD matrix_trans_subC(N_test_total, N_subC);
  TMatrixD T_matrix_trans_subC(N_subC, N_test_total);
  //matrix_trans_subC(old_world_index, new_world_index) = 1;
  matrix_trans_subC(0,0) = 1;
  matrix_trans_subC(1,0) = 1;    
  T_matrix_trans_subC.Transpose( matrix_trans_subC );

  TMatrixD matrix_subC = T_matrix_trans_subC * matrix_total_test * matrix_trans_subC;

  roostr = "canv_matrix_subC";
  TCanvas *canv_matrix_subC = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_matrix_subC, 0.15, 0.2,0.1,0.15);
  matrix_subC.Draw("text");

  /////////// subD
  /////////// sub_matrix = R^T * total_matrix * R
  
  const int N_subD = 2;

  TMatrixD matrix_trans_subD(N_test_total, N_subD);
  TMatrixD T_matrix_trans_subD(N_subD, N_test_total);
  //matrix_trans_subD(old_world_index, new_world_index) = 1;
  matrix_trans_subD(0,0) = 1;
  matrix_trans_subD(1,0) = 1;
  matrix_trans_subD(8,1) = 1;
  matrix_trans_subD(9,1) = 1; 
  T_matrix_trans_subD.Transpose( matrix_trans_subD );

  TMatrixD matrix_subD = T_matrix_trans_subD * matrix_total_test * matrix_trans_subD;

  roostr = "canv_matrix_subD";
  TCanvas *canv_matrix_subD = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_matrix_subD, 0.15, 0.2,0.1,0.15);
  matrix_subD.Draw("text");
  */
}
