#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include <utility>
  #include <iostream>

// ######################################################################
// Various tools and algorithms to make life easier

float GetHistMax(TH1D* h) {
  float max = h->GetBinContent(h->GetMaximumBin() );
  for(int i=1; i<=h->GetXaxis()->GetNbins(); i++){
    float v = h->GetBinContent(i)+h->GetBinError(i);
    max = ( v > max )? v : max;
  }
  return max;
}

float GetHistMin(TH1D* h) {
  float min = h->GetBinContent(h->GetMinimumBin() );
  for(int i=1; i<=h->GetXaxis()->GetNbins(); i++){
    float v = h->GetBinContent(i)-h->GetBinError(i);
    min = ( v < min )? v : min;
  }
  return min;
}


std::pair<float,float> GetHistLimits(TH1D* h) {
  float max = h->GetBinContent(h->GetMaximumBin() );
  float min = h->GetBinContent(h->GetMinimumBin() );
  for(int i=1; i<=h->GetXaxis()->GetNbins(); i++){
    float v = h->GetBinContent(i)+h->GetBinError(i);
    max = ( v > max )? v : max;
    min = ( v < min )? v : min;
  }
  std::pair<float,float> out(min,max);
  return out;

}

void AddTextLine(TLatex* t, float x, float y, int lineNum, std::string line){
  t->DrawLatex( x, y-(lineNum-1)*t->GetTextSize(), line.c_str());
}


void ScaleHist(TH1D* h, std::string mode){
  //float factor = 1./h->GetEntries();
  float factor = 1./h->Integral();
  if( mode == "height" ) factor = 1./GetHistMax(h);
  h->Scale(factor);
  h->SetOption("hist");
}

void NormalizeHist(TH1D* h){
  ScaleHist(h,"");
}

void NormalizeHistInRange(TH1D* h, float x1, float x2){
  int nbins = h->GetXaxis()->GetNbins();
  float sum = 0;
  for(int i=1; i<=nbins; i++){
    float x = h->GetXaxis()->GetBinCenter(i);
    if( x > x1 && x < x2 ) sum += h->GetBinContent(i);
  }
  h->Scale(1/sum);
  h->SetOption("hist");
}


TPaveText* MakeTextBox(float x, float y, float textSize, float numLines, float width){
  TPaveText *pt = new TPaveText(x, y-numLines*textSize, x+width, y, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextSize(textSize);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  pt->SetMargin(0);
  return pt;
}

TPaveText* MakeTextBox(float x, float y, float textSize, float numLines ){
  return MakeTextBox(x,y,textSize,numLines,0.4);
}

void UpdateTPaveSize(TPaveText* pt){
  int nlines = pt->GetListOfLines()->GetSize();
  if( !nlines ) return;
  double line_height = pt->GetTextSize();
  pt->SetY1( pt->GetY2() - nlines*line_height);
}

TLegend* MakeLegend(float x1, float y2, float textSize, float numLines, float width){
    TLegend *leg = new TLegend(x1, y2-numLines*textSize, x1+width, y2 );
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(textSize);
    return leg;
}

TLegend* MakeLegend(float x1, float y2, float textSize, float numLines){
  return MakeLegend(x1, y2, textSize, numLines, 0.3);
}

void AddPointToGraph( TGraph* gr, float x, float y){
    gr->SetPoint(gr->GetN(), x, y );
}

TGraphErrors* HistToGraph(TH1D* h) {
  TGraphErrors* gr = new TGraphErrors();
  for(int i=1; i < h->GetXaxis()->GetNbins(); i++) 
    AddPointToGraph(gr, h->GetXaxis()->GetBinCenter(i), h->GetBinContent(i));
  return gr;
}

void FormatAxes(TH1D* h, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  h->GetXaxis()->SetLabelSize(axisLabelSize);
  h->GetYaxis()->SetLabelSize(axisLabelSize);
  h->GetXaxis()->SetTitleSize(axisTitleSize);
  h->GetYaxis()->SetTitleSize(axisTitleSize);
  h->GetXaxis()->SetTitleOffset(xOffset);
  h->GetYaxis()->SetTitleOffset(yOffset);
}

void FormatAxes(TH2D* h, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  h->GetXaxis()->SetLabelSize(axisLabelSize);
  h->GetYaxis()->SetLabelSize(axisLabelSize);
  h->GetXaxis()->SetTitleSize(axisTitleSize);
  h->GetYaxis()->SetTitleSize(axisTitleSize);
  h->GetXaxis()->SetTitleOffset(xOffset);
  h->GetYaxis()->SetTitleOffset(yOffset);
}
void FormatAxes(TGraphErrors* g, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  g->GetXaxis()->SetLabelSize(axisLabelSize);
  g->GetYaxis()->SetLabelSize(axisLabelSize);
  g->GetXaxis()->SetTitleSize(axisTitleSize);
  g->GetYaxis()->SetTitleSize(axisTitleSize);
  g->GetXaxis()->SetTitleOffset(xOffset);
  g->GetYaxis()->SetTitleOffset(yOffset);
}
void FormatAxes(TMultiGraph* g, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  if( axisLabelSize >= 0 ) {
    g->GetXaxis()->SetLabelSize(axisLabelSize);
    g->GetYaxis()->SetLabelSize(axisLabelSize);
  }
  if( axisTitleSize >= 0 ) {
    g->GetXaxis()->SetTitleSize(axisTitleSize);
    g->GetYaxis()->SetTitleSize(axisTitleSize);
  }
  g->GetXaxis()->SetTitleOffset(xOffset);
  g->GetYaxis()->SetTitleOffset(yOffset);
}
void FormatAxes(TGraphAsymmErrors* g, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  g->GetXaxis()->SetLabelSize(axisLabelSize);
  g->GetYaxis()->SetLabelSize(axisLabelSize);
  g->GetXaxis()->SetTitleSize(axisTitleSize);
  g->GetYaxis()->SetTitleSize(axisTitleSize);
  g->GetXaxis()->SetTitleOffset(xOffset);
  g->GetYaxis()->SetTitleOffset(yOffset);
}

void FormatTH1D(TH1D* h, Color_t lcolor, int lstyle, int lwidth){
  h->SetLineColor(lcolor);
  h->SetLineStyle(lstyle);
  h->SetLineWidth(lwidth);
}

void FormatTH1D(TH1D* h, Color_t lcolor, int lstyle, int lwidth, int mstyle, double msize){
  FormatTH1D(h, lcolor, lstyle, lwidth);
  h->SetMarkerColor(lcolor);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(msize);
}


void FormatTGraph(TGraphAsymmErrors* g, Color_t mc, Color_t lc, int ms, int ls, float msize, int lwidth){
  g->SetMarkerColor(mc);
  g->SetLineColor(lc);
  g->SetMarkerStyle(ms);
  g->SetLineStyle(ls);
  g->SetMarkerSize(msize);
  g->SetLineWidth(lwidth);
}
void FormatTGraph(TGraphErrors* g, Color_t mc, Color_t lc, int ms, int ls, float msize, int lwidth){
  g->SetMarkerColor(mc);
  g->SetLineColor(lc);
  g->SetMarkerStyle(ms);
  g->SetLineStyle(ls);
  g->SetMarkerSize(msize);
  g->SetLineWidth(lwidth);
}
void FormatTGraph(TGraph* g, Color_t mc, Color_t lc, int ms, int ls, float msize, int lwidth){
  g->SetMarkerColor(mc);
  g->SetLineColor(lc);
  g->SetMarkerStyle(ms);
  g->SetLineStyle(ls);
  g->SetMarkerSize(msize);
  g->SetLineWidth(lwidth);
}

void FormatTGraph(TGraphAsymmErrors* g, Color_t mc, Color_t lc, int ms, int ls, float msize){
  FormatTGraph(g,mc,lc,ms,ls,msize,1);
}
void FormatTGraph(TGraphErrors* g, Color_t mc, Color_t lc, int ms, int ls, float msize){
  FormatTGraph(g,mc,lc,ms,ls,msize,1);
}

void CopyTGraphFormat(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, bool copyTitles = false){
  g2->SetMarkerColor(g1->GetMarkerColor());
  g2->SetMarkerStyle(g1->GetMarkerStyle());
  g2->SetMarkerSize(g1->GetMarkerSize());
  g2->SetLineColor(g1->GetLineColor());
  g2->SetLineStyle(g1->GetLineStyle());
  g2->SetLineWidth(g1->GetLineWidth());
  FormatAxes(g2, 
    g1->GetXaxis()->GetTitleSize(),
    g1->GetXaxis()->GetLabelSize(), 
    g1->GetXaxis()->GetTitleOffset(),
    g1->GetYaxis()->GetTitleOffset());
  if( copyTitles ) {
    g2->GetXaxis()->SetTitle( g1->GetXaxis()->GetTitle() );
    g2->GetXaxis()->SetTitleOffset( g1->GetXaxis()->GetTitleOffset() );
    g2->GetYaxis()->SetTitle( g1->GetYaxis()->GetTitle() );
    g2->GetYaxis()->SetTitleOffset( g1->GetYaxis()->GetTitleOffset() );
  }
}

void CopyTGraphFormat(TGraphErrors* g1, TGraphErrors* g2, bool copyTitles = false){
  g2->SetMarkerColor(g1->GetMarkerColor());
  g2->SetMarkerStyle(g1->GetMarkerStyle());
  g2->SetMarkerSize(g1->GetMarkerSize());
  g2->SetLineColor(g1->GetLineColor());
  g2->SetLineStyle(g1->GetLineStyle());
  g2->SetLineWidth(g1->GetLineWidth());
  FormatAxes(g2, 
    g1->GetXaxis()->GetTitleSize(),
    g1->GetXaxis()->GetLabelSize(), 
    g1->GetXaxis()->GetTitleOffset(),
    g1->GetYaxis()->GetTitleOffset());
  if( copyTitles ) {
    g2->GetXaxis()->SetTitle( g1->GetXaxis()->GetTitle() );
    g2->GetXaxis()->SetTitleOffset( g1->GetXaxis()->GetTitleOffset() );
    g2->GetYaxis()->SetTitle( g1->GetYaxis()->GetTitle() );
    g2->GetYaxis()->SetTitleOffset( g1->GetYaxis()->GetTitleOffset() );
  }
}

void CopyHistoFormat(TH1D* h1, TH1D* h2, bool copyTitles = false){
  h2->SetFillColor( h1->GetFillColor() );
  h2->SetFillStyle( h1->GetFillStyle() );
  h2->SetMarkerColor(h1->GetMarkerColor());
  h2->SetMarkerStyle(h1->GetMarkerStyle());
  h2->SetMarkerSize(h1->GetMarkerSize());
  h2->SetLineColor(h1->GetLineColor());
  h2->SetLineStyle(h1->GetLineStyle());
  h2->SetLineWidth(h1->GetLineWidth());
  h2->GetXaxis()->SetTitleSize( h1->GetXaxis()->GetTitleOffset() );
  h2->SetTitleSize( h1->GetTitleSize() );
  FormatAxes(h2, 
    h1->GetXaxis()->GetTitleSize(),
    h1->GetXaxis()->GetLabelSize(), 
    h1->GetXaxis()->GetTitleOffset(),
    h1->GetYaxis()->GetTitleOffset());
  if( copyTitles ) {
    h2->SetTitle( h1->GetTitle() );
    h2->GetXaxis()->SetTitle( h1->GetXaxis()->GetTitle() );
    h2->GetXaxis()->SetTitleOffset( h1->GetXaxis()->GetTitleOffset() );
    h2->GetYaxis()->SetTitle( h1->GetYaxis()->GetTitle() );
    h2->GetYaxis()->SetTitleOffset( h1->GetYaxis()->GetTitleOffset() );
  }
}

void CopyHistoFormat(TH2D* h1, TH2D* h2, bool copyTitles = false){
  h2->SetMarkerColor(h1->GetMarkerColor());
  h2->SetMarkerStyle(h1->GetMarkerStyle());
  h2->SetMarkerSize(h1->GetMarkerSize());
  h2->GetXaxis()->SetTitleSize( h1->GetXaxis()->GetTitleOffset() );
  h2->SetTitleSize( h1->GetTitleSize() );
  FormatAxes(h2, 
    h1->GetXaxis()->GetTitleSize(),
    h1->GetXaxis()->GetLabelSize(), 
    h1->GetXaxis()->GetTitleOffset(),
    h1->GetYaxis()->GetTitleOffset());
  if( copyTitles ) {
    h2->SetTitle( h1->GetTitle() );
    h2->GetXaxis()->SetTitle( h1->GetXaxis()->GetTitle() );
    h2->GetXaxis()->SetTitleOffset( h1->GetXaxis()->GetTitleOffset() );
    h2->GetYaxis()->SetTitle( h1->GetYaxis()->GetTitle() );
    h2->GetYaxis()->SetTitleOffset( h1->GetYaxis()->GetTitleOffset() );
  }
}

TH1D* Clone(TH1D* h){
  TH1D* h2 = (TH1D*)h->Clone();
  h2->SetDirectory(0);
  return h2;
}

TH2D* Clone(TH2D* h){
  TH2D* h2 = (TH2D*)h->Clone();
  h2->SetDirectory(0);
  return h2;
}

//#########################################################################################
void DivideHist(TH1D* h1, TH1D* h2){
  int nbins = h1->GetXaxis()->GetNbins();
  if( h2->GetXaxis()->GetNbins() != nbins ) return;
  for(int i=0; i<=nbins; i++){
    double num = h1->GetBinContent(i);
    double denom = h2->GetBinContent(i);
    double newval = 0;
    if( denom != 0 ) newval = num/denom;
    h1->SetBinContent(i,newval);
    h2->SetBinError(i,0);
  }
}

//##########################################################################################
TGraphErrors* MakeGraph(const TH1D* h){
    std::string name  = h->GetName();
    std::string title = h->GetTitle();
    TGraphErrors* gr = new TGraphErrors();
    gr->SetName((name+"_gr").c_str());
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    gr->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
    for(int i=1; i<=h->GetXaxis()->GetNbins(); i++){
      int pt = gr->GetN();
      float x = h->GetBinCenter(i);
      float y = h->GetBinContent(i);
      if( y == 0 ) continue;
      float dx = h->GetBinWidth(i)/sqrt(12);
      float dy = h->GetBinError(i);
      gr->SetPoint(pt, x, y);
      gr->SetPointError(pt, dx, dy);
    }
  return gr;
}

//##########################################################################################
void ShiftBins(TH1D* h, int shift){
  int bins = h->GetXaxis()->GetNbins();
  TH1D* h2 = (TH1D*)h->Clone("h2");
  h->Reset();
  for(int i=1; i<=bins; i++){
    int prevBin = i - shift;
    if( prevBin < 1 || prevBin > bins ) continue;
    h->SetBinContent(i,h2->GetBinContent(prevBin));
    h->SetBinError(i,h2->GetBinError(prevBin));
  }
  return;
}

//#########################################################################################
TH1D* ScaleXAxis(TH1D* h, float factor, std::string new_name = ""){
  int nbins = h->GetXaxis()->GetNbins();
  float ll = h->GetXaxis()->GetXmin() * factor;
  float ul = h->GetXaxis()->GetXmax() * factor;
  std::string title = h->GetTitle();
  std::string name = h->GetName();
  name += "_scaled";
  if( new_name != "" ) {
    name = new_name;
  } 
  TH1D* h2 = new TH1D(name.c_str(),title.c_str(),nbins,ll,ul);
  for(int i=1; i<=nbins; i++){
    h2->SetBinContent(i, h->GetBinContent(i));
    h2->SetBinError(i, h->GetBinError(i));
  }
  return h2;
}

//##########################################################################################
TH1D* Make1DSlice(const TH2D* h2d, int bin1, int bin2, std::string name)
{
  float nbins = h2d->GetYaxis()->GetNbins();
  float xmin = h2d->GetYaxis()->GetXmin();
  float xmax = h2d->GetYaxis()->GetXmax();
  float mean = 0.5*(h2d->GetXaxis()->GetBinCenter(bin1) + h2d->GetXaxis()->GetBinCenter(bin2));
  name = Form("%s_%d_%d",name.c_str(),bin1,bin2);
  
  TH1D* h;
  if(bin1==bin2)  h = new TH1D(name.c_str(),Form("1D slice: bin %d (x= %.2f)",  bin1,mean),nbins,xmin,xmax);
  else            h = new TH1D(name.c_str(),Form("1D slice: bins %d-%d (x= %f)",bin1,bin2,mean),nbins,xmin,xmax);
  
  h->GetXaxis()->SetTitle( h2d->GetYaxis()->GetTitle() );

  for(int j=1; j<=nbins; j++){
    float sum_bc = 0;
    float sum_sqErr = 0;
    for(int i=bin1; i<=bin2; i++){
      sum_bc += h2d->GetBinContent(i,j);
      sum_sqErr += pow(h2d->GetBinError(i,j),2);
    }
    h->SetBinContent(j, sum_bc);
    h->SetBinError(j, sqrt(sum_sqErr) );
  }
  
  return h;

}

TH1D* Make1DSlice(const TH2D* h2d, int bin1, int bin2){
  return Make1DSlice(h2d,bin1,bin2,"h");
}

// ====================================================================================
// Calculates the chi-squared goodness of match between two histograms. It is assumed
// that the 1st histograms (h1) is data and the 2nd (h2) is MC. 
void CalcChi2(TH1D* h1, TH1D* h2, double& chi2Sum, int& ndf, float x1=-9e9, float x2=9e9){
  int nbins = h1->GetNbinsX();
  if( nbins != h2->GetNbinsX() ) { std::cout<<"Bins don't match!!!\n"; return;}
  chi2Sum = 0.;
  ndf     = 0;
  for(int i=1; i<=nbins; i++){
    float x_l = h1->GetXaxis()->GetBinLowEdge(i);
    float x_u = h1->GetXaxis()->GetBinLowEdge(i+1);
    float x   = h1->GetXaxis()->GetBinCenter(i);
    if( x_l < x1 || x_u > x2 ) continue;
    float b1 = float(h1->GetBinContent(i));
    float b2 = float(h2->GetBinContent(i));
    float e1 = h1->GetBinError(i);
    float e2 = h2->GetBinError(i);
    float sig = sqrt( pow(e1,2) + pow(e2,2) );
    if( sig == 0. ) continue;
    //std::cout<<"bin "<<i<<"  val1= "<<b1<<" +/- "<<e1<<", val2= "<<b2<<" +/- "<<e2<<"\n";
    chi2Sum += pow( (b1-b2)/sig, 2 );
    ndf++;
  }
  std::cout<<"ChiSquare comparison: "<<nbins<<" bins total, "<<ndf<<" bins used, sum_chi2 = "<<chi2Sum<<"\n";
  return;
}

double CalcChi2(TH1D* h1, TH1D* h2, float x1=-9e9, float x2=9e9){
  double chi2 = 0;
  int ndf = 0;
  CalcChi2(h1,h2,chi2,ndf,x1,x2);
  return chi2/ndf;
}

TH1D* GetHistogram(std::string fName, std::string hName, bool normalize){
  TFile* file = new TFile(fName.c_str());
  TH1D* h = (TH1D*)file->Get(hName.c_str());
  if(normalize) NormalizeHist(h);
  return h;
  delete file;
}





