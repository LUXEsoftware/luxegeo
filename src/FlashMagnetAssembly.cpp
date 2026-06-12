//
/// \brief Implementation of the FlashMagnetAssembly class (DD4hep version)
//

#include <cmath>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TGeoManager.h"

#include "LxAux.h"
#include "FlashMagnetAssembly.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
Assembly FlashMagnetAssembly::GetAssembly(dd4hep::Detector& description)
{
  const std::string asmName = fModelName + "_assembly";

  // Check TGeoManager first — if already built, return cached assembly
  TGeoVolume* existing = description.manager().GetVolume(asmName.c_str());
  if (existing) {
    printout(INFO, "FlashMagnetAssembly::GetAssembly",
             "Reusing cached assembly '%s' from TGeoManager.", asmName.c_str());
    return *dynamic_cast<dd4hep::Assembly*>(existing);
  }

  // Not yet built — construct and register via BuildMagnet
  printout(INFO, "FlashMagnetAssembly::GetAssembly",
           "Building assembly '%s'.", asmName.c_str());
  return BuildMagnet(description);
}


// ---------------------------------------------------------------------------
dd4hep::Volume FlashMagnetAssembly::GetFieldVolumeStatic(dd4hep::Detector&  description,
                                                   const std::string& lvname)
{
  Material envMaterial     = description.material(description.constant<std::string>("EnvironmentMaterial"));
  double FlashMFieldX      = description.constant<double>("FlashMFieldX");
  double FlashMagneteffY   = description.constant<double>("FlashMagneteffY");
  double FlashMFieldLength = description.constant<double>("FlashMFieldLength");

  Box    solidFlashMFieldVol(FlashMFieldX/2.0, FlashMagneteffY/2.0, FlashMFieldLength/2.0);
  Volume logicFlashMFieldVol(lvname, solidFlashMFieldVol, envMaterial);
  return logicFlashMFieldVol;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly FlashMagnetAssembly::BuildMagnet(dd4hep::Detector& description)
{
  Material magMaterial  = description.material(description.constant<std::string>("MagnetMaterial"));
  Material wireMaterial = description.material(description.constant<std::string>("TypMBWireMaterial"));

  double FlashMagnetX        = description.constant<double>("FlashMagnetX");
  double FlashMagnetY        = description.constant<double>("FlashMagnetY");
  double FlashMagnetZ        = description.constant<double>("FlashMagnetZ");
  double FlashMMiddleCutX    = description.constant<double>("FlashMMiddleCutX");
  double FlashMMiddleCutY    = description.constant<double>("FlashMMiddleCutY");
  double FlashMagnetCutD     = description.constant<double>("FlashMagnetCutD");
  double FlashMFieldGapY     = description.constant<double>("FlashMFieldGapY");
  double FlashMagnetCoreX    = description.constant<double>("FlashMagnetCoreX");
  double FlashMagnetCoreZ    = description.constant<double>("FlashMagnetCoreZ");
  double FlashMagnetCoilH    = description.constant<double>("FlashMagnetCoilH");
  double FlashMagnetCoilD    = description.constant<double>("FlashMagnetCoilD");
  double FlashMagnetCoilGapH = description.constant<double>("FlashMagnetCoilGapH");

  // -----------------------------------------------------------------------
  // Magnet half yoke
  // -----------------------------------------------------------------------
  double hmagy = FlashMagnetY / 2.0;
  double acutx = 2.0 * FlashMagnetCutD * std::sqrt(2.0);
  double acuty = 2.0 * FlashMagnetCutD / std::sqrt(2.0);

  Box solidFlashMagnetHalf1(FlashMagnetX/2.0, hmagy/2.0, FlashMagnetZ/2.0);
  Box solidFlashMMiddleCut (FlashMMiddleCutX/2.0, FlashMMiddleCutY/2.0, FlashMagnetZ);
  Box solidFlashMACut      (acutx/2.0, acuty/2.0, FlashMagnetZ);

  SubtractionSolid solidFlashMagnetHalf2("solidFlashMagnetHalf2",
      solidFlashMagnetHalf1, solidFlashMMiddleCut,
      Position(0.0, -hmagy/2.0, 0.0));
  SubtractionSolid solidFlashMagnetHalf3("solidFlashMagnetHalf3",
      solidFlashMagnetHalf2, solidFlashMACut,
      Transform3D(RotationZYX(-M_PI/4.0, 0.0, 0.0),
                  Position(FlashMagnetX/2.0, hmagy/2.0, 0.0)));
  SubtractionSolid solidFlashMagnetHalf("solidFlashMagnetHalf",
      solidFlashMagnetHalf3, solidFlashMACut,
      Transform3D(RotationZYX(M_PI/4.0, 0.0, 0.0),
                  Position(-FlashMagnetX/2.0, hmagy/2.0, 0.0)));
  Volume logicFlashMagnetHalf("logicFlashMagnetHalf", solidFlashMagnetHalf, magMaterial);

  // -----------------------------------------------------------------------
  // Core
  // -----------------------------------------------------------------------
  double corey     = 0.5 * (FlashMMiddleCutY - FlashMFieldGapY);
  double corecutd  = 20.0 * dd4hep::mm;
  double acorecutx = 2.0 * corecutd * std::sqrt(2.0);
  double acorecuty = 2.0 * corecutd / std::sqrt(2.0);

  Box solidFlashMagnetCore1  (FlashMagnetCoreX/2.0, corey/2.0, FlashMagnetCoreZ/2.0);
  Box solidFlashMagnetCoreCut(acorecutx/2.0, acorecuty/2.0, FlashMagnetCoreZ);

  SubtractionSolid solidFlashMagnetCore2("solidFlashMagnetCore2",
      solidFlashMagnetCore1, solidFlashMagnetCoreCut,
      Transform3D(RotationZYX(M_PI/4.0, 0.0, 0.0),
                  Position(FlashMagnetCoreX/2.0, -corey/2.0, 0.0)));
  SubtractionSolid solidFlashMagnetCore("solidFlashMagnetCore",
      solidFlashMagnetCore2, solidFlashMagnetCoreCut,
      Transform3D(RotationZYX(-M_PI/4.0, 0.0, 0.0),
                  Position(-FlashMagnetCoreX/2.0, -corey/2.0, 0.0)));
  Volume logicFlashMagnetCore("logicFlashMagnetCore", solidFlashMagnetCore, magMaterial);

  // -----------------------------------------------------------------------
  // Coil
  // -----------------------------------------------------------------------
  double lcoilz = FlashMagnetCoreZ + 2.0 * FlashMagnetCoilGapH;
  double hcoilx = FlashMagnetCoreX + 2.0 * FlashMagnetCoilGapH;

  Box    solidFlashMagnetCoilZ(FlashMagnetCoilH/2.0, FlashMagnetCoilD/2.0, lcoilz/2.0);
  Volume logicFlashMagnetCoilZ("logicFlashMagnetCoilZ", solidFlashMagnetCoilZ, wireMaterial);

  Box    solidFlashMagnetCoilX(hcoilx/2.0, FlashMagnetCoilD/2.0, FlashMagnetCoilH/2.0);
  Volume logicFlashMagnetCoilX("logicFlashMagnetCoilX", solidFlashMagnetCoilX, wireMaterial);

  Tube   solidFlashMCoilR(0.0, FlashMagnetCoilH, FlashMagnetCoilD/2.0, 0.0, 0.5*M_PI);
  Volume logicFlashMCoilR("logicFlashMCoilR", solidFlashMCoilR, wireMaterial);

  // -----------------------------------------------------------------------
  // Half-magnet assembly
  // -----------------------------------------------------------------------
  dd4hep::Assembly magnetHalfAssembly("magnetHalfAssembly");

  magnetHalfAssembly.placeVolume(logicFlashMagnetHalf,
                                 Position(0.0, hmagy/2.0, 0.0));
  magnetHalfAssembly.placeVolume(logicFlashMagnetCore,
                                 Position(0.0, (FlashMMiddleCutY - corey)/2.0, 0.0));

  double coillx = 0.5*(FlashMagnetCoilH + FlashMagnetCoreX) + FlashMagnetCoilGapH;
  double coilly = 0.5*(FlashMMiddleCutY  - FlashMagnetCoilD);
  double coilhz = 0.5*(FlashMagnetZ + FlashMagnetCoilH) + FlashMagnetCoilGapH;

  magnetHalfAssembly.placeVolume(logicFlashMagnetCoilZ, Position( coillx, coilly, 0.0));
  magnetHalfAssembly.placeVolume(logicFlashMagnetCoilZ, Position(-coillx, coilly, 0.0));
  magnetHalfAssembly.placeVolume(logicFlashMagnetCoilX, Position(0.0, coilly,  coilhz));
  magnetHalfAssembly.placeVolume(logicFlashMagnetCoilX, Position(0.0, coilly, -coilhz));

  double coilrx = coillx - FlashMagnetCoilH/2.0;
  double coilrz = coilhz - FlashMagnetCoilH/2.0;
  for (int ii = 0; ii < 4; ++ii) {
    Position crpos(std::pow(-1.0,  ii    & 1) * coilrx,
                   coilly,
                   std::pow(-1.0, (ii>>1)& 1) * coilrz);
    double yangle = std::pow(-1.0, ii) * ((ii & 1) + ((ii >> 1) & 1)) * M_PI/2.0;
    // Note: verified empirically — correct form is RotationZYX(-yangle, 0, pi/2)
    RotationZYX crrot(-yangle, 0.0, M_PI/2.0);
    magnetHalfAssembly.placeVolume(logicFlashMCoilR, Transform3D(crrot, crpos));
  }

  // -----------------------------------------------------------------------
  // Full assembly: two halves, second rotated pi around Z
  // Assembly is named fModelName+"_assembly" so TGeoManager can find it
  // -----------------------------------------------------------------------
  dd4hep::Assembly fMagnetAssembly(fModelName + "_assembly");

  LxAux::AddAssemblyVolumes(fMagnetAssembly, magnetHalfAssembly,
                            Position(0.0, 0.0, 0.0));
  LxAux::AddAssemblyVolumes(fMagnetAssembly, magnetHalfAssembly,
                            Position(0.0, 0.0, 0.0),
                            RotationZYX(M_PI, 0.0, 0.0));

  return fMagnetAssembly;
}


// ---------------------------------------------------------------------------
void FlashMagnetAssembly::ConstructSupport(dd4hep::Detector&       description,
                                           dd4hep::Volume&         motherVol,
                                           const dd4hep::Position& pos,
                                           const std::string&      mname,
                                           bool                    rotate)
{
  double FlashMagnetX    = description.constant<double>("FlashMagnetX");
  double FlashMagnetY    = description.constant<double>("FlashMagnetY");
  double FlashMagnetZ    = description.constant<double>("FlashMagnetZ");
  double FlashMagnetCutD = description.constant<double>("FlashMagnetCutD");
  double FloorSurfaceYpos= description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest     = (description.constant<int>("OverlapTest") != 0);

  double ypestal = 0.5 * dd4hep::m;
  double ylevel  = FlashMagnetY/2.0 - pos.Y();
  double tblx    = FlashMagnetX - 2.0*FlashMagnetCutD;
  double tblz    = FlashMagnetZ - 2.0*FlashMagnetCutD;
  if (rotate) {
    ylevel = FlashMagnetX/2.0 - pos.Y();
    tblx   = FlashMagnetY - 2.0*FlashMagnetCutD;
  }
  double tblhight = -ylevel - ypestal - FloorSurfaceYpos;

  dd4hep::Assembly tablesupport = LxAux::BuildTable(description, mname + "Support",
                                            tblx, tblhight, tblz, 3, 0.0);
  motherVol.placeVolume(tablesupport, Position(pos.X(), -ylevel, pos.Z()));

  Volume pedestal = LxAux::BuildPedestal(description, "Profiler", tblx, ypestal, tblz);
  PlacedVolume pvPed = motherVol.placeVolume(pedestal,
      Position(pos.X(), 0.5*ypestal + FloorSurfaceYpos, pos.Z()));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();
}
