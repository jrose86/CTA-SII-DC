// josie 31 july 2024
// cleaner analysis code to read new ttree format simulation root files, perform full measurement in one go

double* analyzeRun(TTree* header, TTree* pairtree, TTree* teltree, TCanvas* peaks, int ip);

void analysisClean(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TCanvas* peaks = new TCanvas("peaks","all peaks",1600,1400);
  peaks->Divide(4,4);
  peaks->Print("allpeaks.pdf[");

  TCanvas* junk = new TCanvas;
  junk->cd();

  TGraphErrors* viscurve = new TGraphErrors;     viscurve->SetMarkerStyle(20);      viscurve->SetTitle("visibilities vs baseline;baseline (m);V^{2}");
  int npts(0);

  const int nfiles(5);
  TString files[nfiles] = {
    "gam cas_y2024m7d11h0m0s0.root",
    "gam cas_y2024m7d11h1m0s0.root",
    "gam cas_y2024m7d11h2m0s0.root",
    "gam cas_y2024m7d11h6m0s0.root",
    "gam cas_y2024m7d11h8m0s0.root"
  };
  // -----------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // loop over all data files
  for(int i=0; i<nfiles; i++){

    cout << endl << "starting file " << files[i].Data() << endl << endl;

    // open file and assign trees
    TFile* rootfile = new TFile(files[i], "READONLY");
    TTree* header = new TTree;     rootfile->GetObject("Header",header);
    TTree* pairtree = new TTree;   rootfile->GetObject("Pairs",pairtree);
    TTree* teltree = new TTree;    rootfile->GetObject("Telescopes",teltree);

    const int npairs = pairtree->GetEntries();
    const int ntels = teltree->GetEntries();
    cout << "there are " << npairs << " pairs of telescopes and " << ntels << " individual telescopes" << endl;

    // loop running analysis for all pairs in this file
    for(int ip=0; ip<npairs; ip++){
      double* ptarr;
      ptarr = analyzeRun(header, pairtree, teltree, peaks, ip); 
      cout << ip << "   " << ptarr[0] << " " << ptarr[1] << " " << ptarr[2] << " " << ptarr[3] << endl;
      if(ptarr[3] > 1e-9){ // to exclude the one pt with very small errors
	viscurve->AddPoint(ptarr[0], ptarr[2]);
	viscurve->SetPointError(npts, ptarr[1], ptarr[3]);
        npts++;
      }
    }
    peaks->Print("allpeaks.pdf");
  }
  // -----------------------------------------------------------------------------------------------------------------------------------------------------------

  // draw vis vs baseline and fit size of star
  TCanvas* curvecan = new TCanvas;
  viscurve->Draw("AP");

  TF1* udfit = new TF1("udfit","fabs([0])*pow(2.0*TMath::BesselJ1(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265))/(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265)),2)", 0, 180);
  udfit->SetParName(0, "C_{UD}");                 udfit->SetParameter(0, 1e-3);
  udfit->SetParName(1, "#theta_{UD} (mas)");      udfit->SetParameter(1, 0.5);
  viscurve->Fit(udfit,"EX0");

  peaks->Print("allpeaks.pdf]");
    
} // end of main macro


// ===============================================================================================================================================================

double* analyzeRun(TTree* header, TTree* pairtree, TTree* teltree, TCanvas* peaks, int ip){

  TCanvas* junk = new TCanvas;
  junk->cd();

  // get cf and which telescopes
  TProfile2D* cforig = new TProfile2D;   pairtree->SetBranchAddress("CorrFctnBr",&cforig);
  TString* t1id = new TString;           pairtree->SetBranchAddress("TelName1Br",&t1id); // not necessary now but needed in future
  TString* t2id = new TString;           pairtree->SetBranchAddress("TelName2Br",&t2id);
  int t1index, t2index;                  pairtree->SetBranchAddress("TelIndex1Br",&t1index);    pairtree->SetBranchAddress("TelIndex2Br",&t2index);
  pairtree->GetEvent(ip);
  cout << "this pair is " << t1id->Data() << " and " << t2id->Data() << endl;

  // get matching adcs dists, calculate the avg adc of the distribution to use later
  TH2D* adc1dist = new TH2D;             teltree->SetBranchAddress("AdcDistBr",&adc1dist);      teltree->GetEvent(t1index);
  TH2D* adc2dist = new TH2D;             teltree->SetBranchAddress("AdcDistBr",&adc2dist);      teltree->GetEvent(t2index);
  TH1D* avgadc1 = new TH1D("avgadc1","avg ADC 1",adc1dist->GetYaxis()->GetNbins(),0,adc1dist->GetYaxis()->GetXmax());
  TH1D* avgadc2 = new TH1D("avgadc2","avg ADC 2",adc2dist->GetYaxis()->GetNbins(),0,adc2dist->GetYaxis()->GetXmax());

  for(int iy=1; iy<=adc1dist->GetYaxis()->GetNbins(); iy++){
    TH1D* tempadc1 = new TH1D;   tempadc1 = adc1dist->ProjectionX("",iy,iy);
    double adc1mean = tempadc1->GetMean();
    avgadc1->SetBinContent(iy,adc1mean);
    TH1D* tempadc2 = new TH1D;   tempadc2 = adc2dist->ProjectionX("",iy,iy);
    double adc2mean = tempadc2->GetMean();
    avgadc2->SetBinContent(iy,adc2mean);
  }
  // --------------------------------------------------------------------------------------------------------------------------------------------------------
  // normalize by product of avg ADCs (and subtract 1 to center cf at 0)

  TProfile2D* cfnorm = new TProfile2D("cfnorm","cf normalized",cforig->GetXaxis()->GetNbins(),cforig->GetXaxis()->GetXmin(),cforig->GetXaxis()->GetXmax(),cforig->GetYaxis()->GetNbins(),cforig->GetYaxis()->GetXmin(),cforig->GetYaxis()->GetXmax());
  for(int iy=1; iy<=cforig->GetYaxis()->GetNbins(); iy++){
    for(int ix=1; ix<=cforig->GetXaxis()->GetNbins(); ix++){
      cfnorm->Fill(cforig->GetXaxis()->GetBinCenter(ix), cforig->GetYaxis()->GetBinCenter(iy), (cforig->GetBinContent(ix,iy)/(avgadc1->GetBinContent(iy)*avgadc2->GetBinContent(iy))) - 1.0); 
    }
  }
  // --------------------------------------------------------------------------------------------------------------------------------------------------------
  // opd shift and calc avg baseline

  junk->cd();
  TNtuple* geom = new TNtuple;    pairtree->SetBranchAddress("OpdEtcBr",&geom);     pairtree->GetEvent(ip);
  TGraph* opdtemp;   geom->Draw("OpdNs:timeInRunSec");     opdtemp = new TGraph(geom->GetSelectedRows(), geom->GetV1(), geom->GetV2());
  TGraph* utime;     geom->Draw("uM:timeInRunSec");        utime = new TGraph(geom->GetSelectedRows(), geom->GetV2(), geom->GetV1());
  TGraph* vtime;     geom->Draw("vM:timeInRunSec");        vtime = new TGraph(geom->GetSelectedRows(), geom->GetV2(), geom->GetV1());

  // calc weighted avg baseline
  double avgUsum(0.0), avgVsum(0.0), avgDenom(0.0);
  for(int ib=0; ib<utime->GetN(); ib++){
    avgUsum += avgadc1->GetBinContent(ib)*avgadc2->GetBinContent(ib)*utime->GetPointY(ib);
    avgVsum += avgadc1->GetBinContent(ib)*avgadc2->GetBinContent(ib)*vtime->GetPointY(ib);
    avgDenom += avgadc1->GetBinContent(ib)*avgadc2->GetBinContent(ib);
  }
  double baseline = sqrt(pow((avgUsum/avgDenom), 2.0) + pow((avgVsum/avgDenom), 2.0));

  // opd shift 
  TProfile2D* cfshift = new TProfile2D("cfshift",Form("CF %d;relative time (ns);",ip),cfnorm->GetXaxis()->GetNbins()*4.0,cfnorm->GetXaxis()->GetXmin(),cfnorm->GetXaxis()->GetXmax(),cfnorm->GetYaxis()->GetNbins(),cfnorm->GetYaxis()->GetXmin(),cfnorm->GetYaxis()->GetXmin());
  for(int iy=1; iy<=cfnorm->GetYaxis()->GetNbins(); iy++){
    for(int ix=1; ix<=cfnorm->GetXaxis()->GetNbins(); ix++){
      double xshift = cfnorm->GetXaxis()->GetBinCenter(ix) - (opdtemp->GetPointX(iy) - opdtemp->GetPointX(0));
      cfshift->Fill(xshift, cfnorm->GetYaxis()->GetBinCenter(iy), cfnorm->GetBinContent(ix,iy));
    }
  }
  // --------------------------------------------------------------------------------------------------------------------------------------------------------
  // time avg project peak
  TProfile* cfproj = new TProfile;    cfproj->SetTitle(Form("CF %d;relative time (ns);",ip));
  cfproj = cfshift->ProfileX("",1,-1);
  peaks->cd(ip+1);
  cfproj->GetXaxis()->SetRangeUser(-128,128);
  cfproj->Draw();

  // fit projected peak with gaussian
  TF1* hbtfit = new TF1("hbtfit","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-256,256);
  hbtfit->SetParName(0,"area");       hbtfit->SetParameter(0,0.0);
  hbtfit->SetParName(1,"#tau_{0}");   hbtfit->SetParameter(1,0.0);  hbtfit->SetParLimits(1,-10,10);
  hbtfit->SetParName(2,"#sigma");     hbtfit->SetParameter(2,4.0);  hbtfit->SetParLimits(2,2.0,7.0);
  cfproj->Fit(hbtfit);

  junk->cd();

  // array to return vis and baseline to plot 
  static double ptvals[4];
  ptvals[0] = baseline;
  ptvals[1] = 0; // for now
  ptvals[2] = hbtfit->GetParameter(0);
  ptvals[3] = hbtfit->GetParError(0);

  header->ResetBranchAddresses(); // !!!!! THIS fixes the seg fault
  pairtree->ResetBranchAddresses();
  teltree->ResetBranchAddresses();

  return ptvals;
}
