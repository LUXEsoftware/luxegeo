
//
/// \brief Implementation of the PrimarySpectra class
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <functional>

#include "Randomize.hh"

#include "PrimarySpectra.h"


PrimarySpectra::PrimarySpectra() : vx(0), vy(0), vfc(0), fverbose(0)
{}


PrimarySpectra::PrimarySpectra(const std::string fname) : vx(0), vy(0), vfc(0), fverbose(0)
{ 
  LoadData(fname);  
}


int PrimarySpectra::LoadData(const std::string fname)
{
  std::fstream    fin;
  int nlread, fioerr;

  fin.open(fname.c_str(), std::ios::in);
  if (!fin.is_open()) {
    throw std::runtime_error(std::string("Error open file ") + fname);
  }
  vx.clear();
  vy.clear();
  vfc.clear();
  nlread = fioerr = 0;
  while (!fin.eof()) {
    std::string datstr;
    std::getline(fin, datstr);
    if (fin.fail()) {
      if (!fin.eof()) { ++fioerr; std::cout << "Error reading data from the file " << fname << "!" << std::endl; }
      break;    
    }
    double xx, yy;
    std::stringstream sstr(datstr);
    sstr >> xx >> yy;
    if (sstr.fail()) {
      ++fioerr; 
      std::cout << "Error exracting data from the file " << fname << "!" << std::endl;
      break;    
    }

    if (fverbose) {
      std::cout << xx << " " << yy << std::endl;
    }
    vx.push_back(xx);
    vy.push_back(yy);
    ++nlread;
  }
  fin.close();
  if (fioerr) throw std::runtime_error(std::string("Error reding data from the file ") + fname);
  
  FindCumulative();
  
  if (fverbose) {
    std::cout << nlread << " rows were successfully read from the file " << fname << std::endl;
  }
  return nlread;
}


void PrimarySpectra::FindCumulative()
{
  vfc.resize(vx.size());
  double fc = 0.0;
  vfc[0] = fc;
  for (size_t i = 1; i < vx.size(); ++i) {
    fc += 0.5 * (vy[i]+vy[i-1]) * (vx[i]-vx[i-1]);
    vfc[i] = fc;
  }
  
  if (fc > 0.0) std::for_each(vfc.begin(), vfc.end(), [=](double &x) {x /= fc;});
}


double PrimarySpectra::GetRandom()
{
//   if (!vfc.size()) throw std::runtime_error("Can not generate a random number. The PDF was not loaded!");
//   double x = gRandom->Uniform(1.0);
  double x = G4RandFlat::shoot();
  
  auto itr = std::find_if( vfc.begin(), vfc.end(), std::bind(std::less_equal<double>(), x, std::placeholders::_1) );
  if (itr == vfc.begin())  return vx.front();
  if (itr == vfc.end()  )  return vx.back();

  int ii = itr - vfc.begin();
  double res = (vx[ii]-vx[ii-1]) / (vfc[ii]-vfc[ii-1]) * (x - vfc[ii-1]) + vx[ii-1];
  return res * fscale;
}

  

// int PrimarySpectraTest(const std::string fname)
// {
//   PrimarySpectra fsptr;
//   fsptr.LoadData(fname);
// 
//   std::vector<double>   vx, vy, vfc;
//   fsptr.GetPDF(vx, vy);
//   fsptr.GetCumulative(vx, vfc);
// 
//   TH1D *hh = new TH1D("pdf", "pdf", 100, vx.front(), vx.back());
//   for (int jj = 0; jj < 1000000; ++jj) {
//     double rr = fsptr.GetRandom();  
// //   std::cout << rr << std::endl; 
//     hh->Fill(rr);  
//   }
//   
// 
//   TGraph *gd = new TGraph(vx.size(), &vx[0], &vy[0]);
//   TGraph *gc = new TGraph(vx.size(), &vx[0], &vfc[0]);
//   
//   TCanvas   *c4 = new TCanvas("RndGenTest", "RndGenTest", 1000, 500);
//   c4->Divide(3,1);
//   c4->cd(1);
//   gd->Draw("LA");
//   
//   c4->cd(2);
//   gc->Draw("LA");
//   
//   c4->cd(3);
//   hh->Draw();
//   
//   return 0;
// }



