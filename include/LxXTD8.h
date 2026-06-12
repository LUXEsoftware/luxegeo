#ifndef LxXTD8_h
#define LxXTD8_h 1

#include "DD4hep/DetFactoryHelper.h"

class LxXTD8
{
  public:
    LxXTD8() {};
    virtual ~LxXTD8() {};
    virtual void Construct(dd4hep::Detector& description,
                           dd4hep::DetElement& sdet,
                           xml_h& e);

  protected:
    void CreateMaterial();
};

#endif
