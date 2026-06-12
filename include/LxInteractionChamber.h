//
/// \brief LxInteractionChamber class (DD4hep version)
//

#ifndef LxInteractionChamber_dd4hep_h
#define LxInteractionChamber_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxInteractionChamber
{
public:
  LxInteractionChamber() {};
  virtual ~LxInteractionChamber() {};

  void Construct(dd4hep::Detector&   description,
                 dd4hep::DetElement& sdet,
                 xml_h&              e);

protected:
  void ConstructMirrors(dd4hep::Detector& description,
                        dd4hep::Volume&   logicTAUICContainer);

  dd4hep::Assembly ConstructMirrorAssembly(dd4hep::Detector&  description,
                                           double r, double d, double h,
                                           const std::string& mname);
};

#endif
