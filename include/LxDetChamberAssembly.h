//
/// \brief LxDetChamberAssembly class (DD4hep version)
//

#ifndef LxDetChamberAssembly_dd4hep_h
#define LxDetChamberAssembly_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxDetChamberAssembly
{
public:
  LxDetChamberAssembly(const std::string& typeName);
  virtual ~LxDetChamberAssembly() {};

  void Construct(dd4hep::Detector&   description,
                 dd4hep::DetElement& sdet,
                 xml_h&              e);

protected:
  dd4hep::Assembly ConstructVacuumChamber(dd4hep::Detector& description);
  dd4hep::Assembly ConstructFrameAssembly(dd4hep::Detector& description, double vcLength);

protected:
  std::string fChamberType;
};

#endif
