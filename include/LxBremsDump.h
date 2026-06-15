//
/// \brief LxBremsDump class (DD4hep version)
//

#ifndef LxBremsDump_dd4hep_h
#define LxBremsDump_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxBremsDump
{
public:
  LxBremsDump() {};
  virtual ~LxBremsDump() {};

  void Construct(dd4hep::Detector&   description,
                 dd4hep::DetElement& sdet,
                 xml_h&              e);

private:
  void ConstructBeamDump(dd4hep::Detector& description, dd4hep::Volume& motherVol);
  void ConstructShielding(dd4hep::Detector& description, dd4hep::Volume& motherVol);
};

#endif
