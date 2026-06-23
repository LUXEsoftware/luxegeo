//
/// \brief DD4hep plugins for placing target chambers with specific targets.
///
/// Two detector types:
///   type="LxBremsTarget"  — bremsstrahlung target chamber
///   type="LxGammaTarget"  — gamma target chamber (wire or foil)
///
/// Architecture mirrors LxMagnet_geo: the chamber body is shared and cached
/// in TGeoManager; the target volume is fresh per placement.
///
/// Example XML:
///   <detector id="11" name="BremsTargetChamber" type="LxBremsTarget">
///     <target zpos="BTargetZpos" targetmaterial="BTargetMaterial"
///             targetz="BTargetZ"/>
///   </detector>
///
///   <detector id="12" name="GammaTargetChamber" type="LxGammaTarget">
///     <target zpos="GTargetZpos" targetmaterial="GTargetMaterial"
///             targetz="GTargetZ" targettype="foil"/>
///   </detector>
//

#include <stdexcept>
#include <string>
#include <cmath>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TGeoManager.h"

#include "LxTargetChamber.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
// Plugin: LxBremsTarget
// ---------------------------------------------------------------------------
static Ref_t create_LxBremsTarget(dd4hep::Detector& description,
                                   xml_h             e,
                                   dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  Volume       fLogicWorld = description.worldVolume();
  Assembly     envelope(detName + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", detID);
  sdet.setPlacement(envPV);

  bool OverlapTest = (description.constant<int>("OverlapTest") != 0);

  // Read target-specific parameters from <target .../> child element
  xml_comp_t x_target = x_det.child(xml_tag_t("target"));
  double      zpos           = x_target.attr<double>(xml_tag_t("zpos"));
  std::string targetMatName  = x_target.attr<std::string>(xml_tag_t("targetmaterial"));
  double      targetZ        = x_target.attr<double>(xml_tag_t("targetz"));

  double TargetChamberXPos = description.constant<double>("TargetChamberXPos");
  double ChamberTargetVolX = description.constant<double>("ChamberTargetVolX");
  double ChamberTargetVolY = description.constant<double>("ChamberTargetVolY");

  Material targetMaterial = description.material(
      description.constant<std::string>(targetMatName));

  // -----------------------------------------------------------------------
  // Chamber body (shared, cached)
  // -----------------------------------------------------------------------
  LxTargetChamber tc;
  Volume chamberVol = tc.GetChamberVolume(description);

  PlacedVolume pvChamber = envelope.placeVolume(chamberVol,
      Position(TargetChamberXPos, 0.0, zpos));
  pvChamber.addPhysVolID("BremsTargetChamber", 0);
  if (OverlapTest) pvChamber.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Target volume (fresh per placement)
  // -----------------------------------------------------------------------
  Volume targetContVol = LxTargetChamber::GetTargetVolume(description, "BTarget");

  // Brems target: simple foil (box)
  Box    solidBremsTarget(ChamberTargetVolX/2.0, ChamberTargetVolY/2.0, targetZ/2.0);
  Volume logicBremsTarget("logicBremsTarget", solidBremsTarget, targetMaterial);

  PlacedVolume pvTarget = targetContVol.placeVolume(logicBremsTarget, Position(0.0, 0.0, 0.0));
  pvTarget.addPhysVolID("BremsTarget", 0);
  if (OverlapTest) pvTarget.ptr()->CheckOverlaps();

  // Place target container into world at same zpos (no X offset — container
  // is centred on beam axis, target position inside is defined by ChamberTargetVolX offset)
  PlacedVolume pvTargetCont = envelope.placeVolume(targetContVol,
      Position(0.0, 0.0, zpos));
  pvTargetCont.addPhysVolID("BremsTargetContainer", 0);
  if (OverlapTest) pvTargetCont.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Support
  // -----------------------------------------------------------------------
  LxTargetChamber::ConstructSupport(description, envelope, "BremsTargetChamber", zpos);

  return sdet;
}


// ---------------------------------------------------------------------------
// Plugin: LxGammaTarget
// ---------------------------------------------------------------------------
static Ref_t create_LxGammaTarget(dd4hep::Detector& description,
                                   xml_h             e,
                                   dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  Volume       fLogicWorld = description.worldVolume();
  Assembly     envelope(detName + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", detID);
  sdet.setPlacement(envPV);

  bool OverlapTest = (description.constant<int>("OverlapTest") != 0);

  xml_comp_t  x_target      = x_det.child(xml_tag_t("target"));
  double      zpos           = x_target.attr<double>(xml_tag_t("zpos"));
  std::string targetMatName  = x_target.attr<std::string>(xml_tag_t("targetmaterial"));
  double      targetZ        = x_target.attr<double>(xml_tag_t("targetz"));
  std::string targetType     = x_target.hasAttr(xml_tag_t("targettype"))
                               ? x_target.attr<std::string>(xml_tag_t("targettype")) : "foil";

  double TargetChamberXPos = description.constant<double>("TargetChamberXPos");
  double ChamberTargetVolX = description.constant<double>("ChamberTargetVolX");
  double ChamberTargetVolY = description.constant<double>("ChamberTargetVolY");

  Material targetMaterial = description.material(
      description.constant<std::string>(targetMatName));

  // -----------------------------------------------------------------------
  // Chamber body (shared, cached — same volume as bremsstrahlung chamber)
  // -----------------------------------------------------------------------
  LxTargetChamber tc;
  Volume chamberVol = tc.GetChamberVolume(description);

  // copy number 1 — bremsstrahlung chamber uses copy 0
  PlacedVolume pvChamber = envelope.placeVolume(chamberVol,
      Position(TargetChamberXPos, 0.0, zpos));
  pvChamber.addPhysVolID("GammaTargetChamber", 1);
  if (OverlapTest) pvChamber.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Target volume (fresh per placement)
  // -----------------------------------------------------------------------
  Volume targetContVol = LxTargetChamber::GetTargetVolume(description, "GammaTarget");

  if (targetType == "wire") {
    // Wire target: Tubs rotated around X by pi/2 so wire axis lies along Y
    // G4RotationMatrix(0, pi/2, 0) is ZXZ Euler → rotation around X by pi/2
    // → RotationZYX(0, 0, pi/2)
    Tube   solidWire(0.0, targetZ/2.0, ChamberTargetVolY, 0.0, 2.0*M_PI);
    Volume logicWire("logicGWire", solidWire, targetMaterial);
    PlacedVolume pvWire = targetContVol.placeVolume(logicWire,
        Transform3D(RotationZYX(0.0, 0.0, M_PI/2.0), Position(0.0, 0.0, 0.0)));
    pvWire.addPhysVolID("GammaTarget", 0);
    if (OverlapTest) pvWire.ptr()->CheckOverlaps();
  } else {
    // Foil target: box
    Box    solidFoil(ChamberTargetVolX/2.0, ChamberTargetVolY/2.0, targetZ/2.0);
    Volume logicFoil("logicGFoil", solidFoil, targetMaterial);
    PlacedVolume pvFoil = targetContVol.placeVolume(logicFoil, Position(0.0, 0.0, 0.0));
    pvFoil.addPhysVolID("GammaTarget", 0);
    if (OverlapTest) pvFoil.ptr()->CheckOverlaps();
  }

  PlacedVolume pvTargetCont = envelope.placeVolume(targetContVol,
      Position(0.0, 0.0, zpos));
  pvTargetCont.addPhysVolID("GammaTargetContainer", 0);
  if (OverlapTest) pvTargetCont.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Support
  // -----------------------------------------------------------------------
  LxTargetChamber::ConstructSupport(description, envelope, "GammaTargetChamber", zpos);

  return sdet;
}


DECLARE_DETELEMENT(LxBremsTarget, create_LxBremsTarget)
DECLARE_DETELEMENT(LxGammaTarget, create_LxGammaTarget)
