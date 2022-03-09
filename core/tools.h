#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLine.h"

// ######################################################################
// Various tools and algorithms to make life easier

float GetHistMax(TH1D* h) {
  return h->GetBinContent(h->GetMaximumBin() );
}

void AddTextLine(TLatex* t, float x, float y, int lineNum, std::string line){
  t->DrawLatex( x, y-(lineNum-1)*t->GetTextSize(), line.c_str());
}

void ScaleHist(TH1D* h, std::string mode){
  float factor = 1./h->GetEntries();
  if( mode == "height" ) factor = 1./GetHistMax(h);
  h->Scale(factor);
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

void FormatTH1D(TH1D* h, Color_t lc, int ls, int lwidth){
  h->SetLineColor(lc);
  h->SetLineStyle(ls);
  h->SetLineWidth(lwidth);
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

//##########################################################################################
TH1D* Make1DSlice(const TH2D* h2d, int bin1, int bin2, std::string name)
{
  float nbins = h2d->GetYaxis()->GetNbins();
  float xmin = h2d->GetYaxis()->GetXmin();
  float xmax = h2d->GetYaxis()->GetXmax();
  float mean = (h2d->GetXaxis()->GetBinCenter(bin1) + h2d->GetXaxis()->GetBinCenter(bin2))/2.;
  name = Form("%s_%d_%d",name.c_str(),bin1,bin2);
  TH1D* h = new TH1D(name.c_str(),Form("1D slice: bins %d-%d (x= %f)",bin1,bin2,mean),nbins,xmin,xmax);
  for(int i=bin1; i<= bin2; i++){
    for(int j=1; j<nbins; j++){
      int bincontent = h2d->GetBinContent(i,j);
      for(int jj=0; jj<bincontent; jj++) h->Fill(h2d->GetYaxis()->GetBinCenter(j)); //h2d->GetBinContent(i,j));
    }
  }

  return h;

}

TH1D* Make1DSlice(const TH2D* h2d, int bin1, int bin2){
  return Make1DSlice(h2d,bin1,bin2,"h");
}



