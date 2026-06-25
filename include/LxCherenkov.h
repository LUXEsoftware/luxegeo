//
/// \brief LxCherenkov class (DD4hep version)
//
/// Architecture: single plugin function DECLARE_DETELEMENT(LxCherenkov,...)
/// is called once per <detector type="LxCherenkov"> XML entry.
/// On the first call the class instance is built and stored in description
/// via addExtension<LxCherenkov>. Subsequent calls retrieve the stored
/// instance and place the already-built CerenkovMotherLogical at the new
/// position. Each call produces one DetElement with a unique detector ID.
//

#ifndef LxCherenkov_dd4hep_h
#define LxCherenkov_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxCherenkov
{
public:
  LxCherenkov() : fBuilt(false) {};
  virtual ~LxCherenkov() {};

  typedef std::map<std::string, dd4hep::Transform3D> TransformMapT;

  // Build internal volumes and straw template — called once on first plugin invocation.
  void Build(dd4hep::Detector& description, dd4hep::SensitiveDetector& sd);

  // Place CerenkovMotherLogical into motherVol at the given transform.
  // Creates and returns a DetElement for this placement with straw DetElements attached.
  dd4hep::DetElement Place(dd4hep::Detector&   description,
                           dd4hep::Volume&     motherVol,
                           const std::string&  detName,
                           int                 detID
                           /*, const dd4hep::Transform3D& trf*/);

  // Half-dimensions of CerenkovMotherLogical — useful for computing positions externally
  double halfX() const { return fHalfX; }
  double halfZ() const { return fHalfZ; }
  double totalBoxHeight() const { return fTotalBoxHeight; }

protected:
  void PrintCherenkovTranslations();

private:
  bool                 fBuilt;
  dd4hep::Volume       fCerenkovMotherVol;
  dd4hep::DetElement   fStrawTemplate;   // cloned per straw per placement
  int                  fNStraws;
  double               fHalfX;
  double               fHalfZ;
  double               fTotalBoxHeight;

  TransformMapT        fDetTransformMap;
};

#endif
