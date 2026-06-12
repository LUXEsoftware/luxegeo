//
/// \brief FlashMagnetAssembly class (DD4hep version)
//

#ifndef FlashMagnetAssembly_dd4hep_h
#define FlashMagnetAssembly_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"
#include "LxMagnetAssembly.h"

class FlashMagnetAssembly : public LxMagnetAssembly
{
public:
  FlashMagnetAssembly() : LxMagnetAssembly("FlashMagnet") {};
  virtual ~FlashMagnetAssembly() {};

  // Returns the shared assembly; builds once, caches in TGeoManager.
  virtual dd4hep::Assembly GetAssembly(dd4hep::Detector& description) override;

  // Returns a fresh field volume with the given unique name.
  // Static so it can be called without an instance.
  static dd4hep::Volume GetFieldVolumeStatic(dd4hep::Detector&  description,
                                       const std::string& lvname);

  // Non-static virtual override — delegates to static version.
  virtual dd4hep::Volume GetFieldVolume(dd4hep::Detector&  description,
                                        const std::string& lvname) override
  { return FlashMagnetAssembly::GetFieldVolumeStatic(description, lvname); }

  virtual void ConstructSupport(dd4hep::Detector&       description,
                                dd4hep::Volume&         motherVol,
                                const dd4hep::Position& pos,
                                const std::string&      mname,
                                bool                    rotate = false) override;

protected:
  dd4hep::Assembly BuildMagnet(dd4hep::Detector& description);
};

#endif
