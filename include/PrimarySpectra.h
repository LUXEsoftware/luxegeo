//
//
/// \brief Definition of the PrimarySpectra class
//

#ifndef PRIMARYSPECTRA_H
#define PRIMARYSPECTRA_H

#include <vector>

class PrimarySpectra
{
public:    
  PrimarySpectra();  
  PrimarySpectra(const std::string fname);
  ~PrimarySpectra(){};  
  int LoadData(const std::string fname);
  double GetRandom();
  void SetScale(const double c = 1.0) { fscale = c; }
  void SetVerbose(const int v = 1) { fverbose = v; }
  void GetPDF(std::vector<double>  &x, std::vector<double>  &y) const { x = vx; y = vy; }
  void GetCumulative(std::vector<double>  &x, std::vector<double>  &y) const { x = vx; y = vfc; }
  
protected:
  void FindCumulative();   
    
  
private:  
  std::vector<double>   vx, vy, vfc;
  double                fscale;
  int fverbose;
};

#endif


