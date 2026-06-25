//
/// \brief Implementation of the LxCherenkov class (DD4hep version)
//

#include <cmath>
#include <string>
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxCherenkov.h"

using namespace dd4hep;

static const std::string kCherenkovVolName = "CerenkovMotherLogical";


// ---------------------------------------------------------------------------
void LxCherenkov::Build(dd4hep::Detector& description, dd4hep::SensitiveDetector& sd)
{
  if (fBuilt) return;
  fBuilt = true;

  // -----------------------------------------------------------------------
  // Materials
  // -----------------------------------------------------------------------
  Material envMaterial          = description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material cerenkovMedium       = description.material(description.constant<std::string>("CerenkovMedium"));
  Material cerenkovMetal        = description.material(description.constant<std::string>("CerenkovMetal"));
  Material graphiteMaterial     = description.material(description.constant<std::string>("CerenkovStrawGraphiteLayerMaterial"));
  Material alMaterial           = description.material(description.constant<std::string>("CerenkovStrawAlLayerMaterial"));
  Material kaptonMaterial       = description.material(description.constant<std::string>("CerenkovStrawKaptonLayerMaterial"));
  Material polyurethaneMaterial = description.material(description.constant<std::string>("CerenkovStrawPolyurethaneLayerMaterial"));
  Material shieldingMaterial    = description.material(description.constant<std::string>("CerenkovShieldingPlateMaterial"));
  Material electronicsMaterial  = description.material(description.constant<std::string>("CerenkovElectronicsBoardMaterial"));

  // -----------------------------------------------------------------------
  // Constants
  // -----------------------------------------------------------------------
  double csh  = description.constant<double>("CerenkovStrawHeight");
  double sir  = description.constant<double>("CerenkovStrawInnerRadius");
  double sglt = description.constant<double>("CerenkovStrawGraphiteLayerThickness");
  double salt = description.constant<double>("CerenkovStrawAlLayerThickness");
  double sklt = description.constant<double>("CerenkovStrawKaptonLayerThickness");
  double splt = description.constant<double>("CerenkovStrawPolyurethaneLayerThickness");
  double cspt = description.constant<double>("CerenkovShieldingPlateThickness");
  double cbt  = description.constant<double>("CerenkovBoxThickness");
  double cspb = description.constant<double>("CerenkovStrawPlateBuffer");
  double cst  = description.constant<double>("CerenkovSupportThickness");
  double cebt = description.constant<double>("CerenkovElectronicsBoardThickness");
  double csxf = description.constant<double>("CerenkovStrawXFrequency");
  double cszf = description.constant<double>("CerenkovStrawZFrequency");
  double cslo = description.constant<double>("CerenkovStrawLayerOffset");
  double CerenkovTotalBoxHeight      = description.constant<double>("CerenkovTotalBoxHeight");
  double CerenkovBeamWindowY         = description.constant<double>("CerenkovBeamWindowY");
  double CerenkovBeamWindowThickness = description.constant<double>("CerenkovBeamWindowThickness");
  double CerenkovSupportPosYtoStage  = description.constant<double>("CerenkovSupportPosYtoStage");
  double ComptonElectronBeamtoStageY = description.constant<double>("ComptonElectronBeamtoStageY");
  bool   OverlapTest = (description.constant<int>("OverlapTest") != 0);

  int cc  = (int)description.constant<double>("CerenkovChannels");
  int ccl = (int)description.constant<double>("CerenkovChannelLayers");

  // Store for external use
  fNStraws       = cc;
  fHalfX         = (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. + sir + 2*(sglt+salt+sklt+splt) + cspb;
  fHalfZ         = (ccl-1)*cszf/2.       + sir + 2*(sglt+salt+sklt+splt) + cspb;
  fTotalBoxHeight= CerenkovTotalBoxHeight;

  // -----------------------------------------------------------------------
  // Straw substructure
  // -----------------------------------------------------------------------
  double sor = sir + 2*(sglt+salt+sklt+splt);  // straw outer radius

  Tube solidStrawMother(0.0, sor,                              csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawInner (0.0, sir,                              csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawGrIn  (sir,           sir+sglt,               csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawAlIn  (sir+sglt,      sir+sglt+salt,          csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawKapIn (sir+sglt+salt, sir+sglt+salt+sklt,     csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawPoly  (sir+sglt+salt+sklt,
                        sir+sglt+salt+sklt+2*splt,             csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawKapOut(sir+sglt+salt+sklt+2*splt,
                        sir+sglt+salt+2*(sklt+splt),           csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawAlOut (sir+sglt+salt+2*(sklt+splt),
                        sir+sglt+2*(salt+sklt+splt),           csh/2., 0.0, 2.0*M_PI);
  Tube solidStrawGrOut (sir+sglt+2*(salt+sklt+splt), sor,      csh/2., 0.0, 2.0*M_PI);

  Volume CerenkovStrawMotherLogical("CerenkovStrawMotherLogical", solidStrawMother, cerenkovMedium);
  Volume CerenkovStrawInnerLogical ("CerenkovStrawInnerLogical",  solidStrawInner,  cerenkovMedium);
  Volume CerenkovStrawGrInLogical  ("CerenkovStrawGraphiteInnerLayerLogical",  solidStrawGrIn,   graphiteMaterial);
  Volume CerenkovStrawAlInLogical  ("CerenkovStrawAlInnerLayerLogical",        solidStrawAlIn,   alMaterial);
  Volume CerenkovStrawKapInLogical ("CerenkovStrawKaptonInnerLayerLogical",    solidStrawKapIn,  kaptonMaterial);
  Volume CerenkovStrawPolyLogical  ("CerenkovStrawPolyurethaneLayerLogical",   solidStrawPoly,   polyurethaneMaterial);
  Volume CerenkovStrawKapOutLogical("CerenkovStrawKaptonOuterLayerLogical",    solidStrawKapOut, kaptonMaterial);
  Volume CerenkovStrawAlOutLogical ("CerenkovStrawAlOuterLayerLogical",        solidStrawAlOut,  alMaterial);
  Volume CerenkovStrawGrOutLogical ("CerenkovStrawGraphiteOuterLayerLogical",  solidStrawGrOut,  graphiteMaterial);

  // Mark inner gas volume as sensitive
  CerenkovStrawInnerLogical.setSensitiveDetector(sd);

  // Place inner gas (sensitive) into straw mother — once, outside any loop
  PlacedVolume pvInner = CerenkovStrawMotherLogical.placeVolume(
      CerenkovStrawInnerLogical, Position(0., 0., 0.));
  // No addPhysVolID — straw number assigned at straw mother placement level

  // Straw template DetElement — cloned per straw per Cherenkov placement
  fStrawTemplate = DetElement("strawInner", 0);
  fStrawTemplate.setPlacement(pvInner);

  // Layer shells
  CerenkovStrawMotherLogical.placeVolume(CerenkovStrawGrInLogical,   Position(0.,0.,0.));
  CerenkovStrawMotherLogical.placeVolume(CerenkovStrawAlInLogical,   Position(0.,0.,0.));
  CerenkovStrawMotherLogical.placeVolume(CerenkovStrawKapInLogical,  Position(0.,0.,0.));
  CerenkovStrawMotherLogical.placeVolume(CerenkovStrawPolyLogical,   Position(0.,0.,0.));
  CerenkovStrawMotherLogical.placeVolume(CerenkovStrawKapOutLogical, Position(0.,0.,0.));
  CerenkovStrawMotherLogical.placeVolume(CerenkovStrawAlOutLogical,  Position(0.,0.,0.));
  CerenkovStrawMotherLogical.placeVolume(CerenkovStrawGrOutLogical,  Position(0.,0.,0.));

  // -----------------------------------------------------------------------
  // Housing volumes
  // -----------------------------------------------------------------------
  Box solidCerenkovMother(fHalfX + cbt + cspt,
                          CerenkovTotalBoxHeight/2.,
                          fHalfZ + cbt + cspt);
  Box solidCerenkovBoxOuter(fHalfX + cbt, CerenkovTotalBoxHeight/2., fHalfZ + cbt);
  Box solidCerenkovBoxInner(fHalfX,       CerenkovTotalBoxHeight/2. - cbt, fHalfZ);
  Box solidCerenkovBeamWindow(fHalfX, CerenkovBeamWindowY/2.,
                              fHalfZ + cbt - CerenkovBeamWindowThickness);

  SubtractionSolid solidCerenkovBox1("CerenkovBoxSolid",
      solidCerenkovBoxOuter, solidCerenkovBoxInner, Position(0.,0.,0.));
  SubtractionSolid solidCerenkovBox("CerenkovBoxSolid2",
      solidCerenkovBox1, solidCerenkovBeamWindow,
      Position(0., ComptonElectronBeamtoStageY - CerenkovTotalBoxHeight/2., 0.));

  Box solidSupportPlate(fHalfX, cst/2.,  fHalfZ);
  Box solidElecBoard   (fHalfX, cebt/2., fHalfZ);

  Box solidShieldUpper(fHalfX + cbt + cspt,
                       (CerenkovTotalBoxHeight - ComptonElectronBeamtoStageY - CerenkovBeamWindowY/2.)/2.,
                       cspt/2.);
  Box solidShieldLower(fHalfX + cbt + cspt,
                       (ComptonElectronBeamtoStageY - CerenkovBeamWindowY/2.)/2.,
                       cspt/2.);
  Box solidShieldSide (cspt/2., CerenkovTotalBoxHeight/2., fHalfZ + cbt + cspt);

  Volume CerenkovBoxLogical     ("CerenkovBoxLogical",              solidCerenkovBox,      cerenkovMetal);
  Volume CerenkovBoxInnerLogical("CerenkovBoxInnerLogical",         solidCerenkovBoxInner, cerenkovMedium);
  Volume CerenkovSupportLogical ("CerenkovSupportPlateLogical",     solidSupportPlate,     cerenkovMetal);
  Volume CerenkovElecLogical    ("CerenkovElectronicsBoardLogical", solidElecBoard,        electronicsMaterial);
  Volume ShieldUpperLogical     ("CerenkovShieldingPlateUpperLogical", solidShieldUpper,   shieldingMaterial);
  Volume ShieldLowerLogical     ("CerenkovShieldingPlateLowerLogical", solidShieldLower,   shieldingMaterial);
  Volume ShieldSideLogical      ("CerenkovShieldingPlateSideLogical",  solidShieldSide,    shieldingMaterial);

  // -----------------------------------------------------------------------
  // Straw placement loop into CerenkovBoxInnerLogical
  // G4RotationMatrix(G4ThreeVector(-1,0,0), pi/2) → RotationZYX(0, 0, -pi/2)
  // rotates tube axis Z → -Y so straw stands vertically along Y
  // -----------------------------------------------------------------------
  RotationZYX strawRot(0.0, 0.0, -M_PI/2.0);
  double straw_y = -CerenkovTotalBoxHeight/2. + csh/2. + cst + CerenkovSupportPosYtoStage;
  int strawIdx = 0;

  for (int i = 0; i < ccl; ++i) {
    double straw_z = -(ccl-1)*cszf/2. + i*cszf;
    for (int j = 0; j < cc/ccl; ++j) {
      double straw_x = -(cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. - i*cslo + j*csxf;
      PlacedVolume pvStraw = CerenkovBoxInnerLogical.placeVolume(
          CerenkovStrawMotherLogical,
          Transform3D(strawRot, Position(straw_x, straw_y, straw_z)));
      pvStraw.addPhysVolID("straw", strawIdx++);
      if (OverlapTest) pvStraw.ptr()->CheckOverlaps();
    }
  }

  // Support and electronics plates
  CerenkovBoxInnerLogical.placeVolume(CerenkovSupportLogical,
      Position(0., straw_y - csh/2. - cst/2., 0.)).addPhysVolID("CerenkovSupportPlate", 0);
  CerenkovBoxInnerLogical.placeVolume(CerenkovSupportLogical,
      Position(0., straw_y + csh/2. + cst/2., 0.)).addPhysVolID("CerenkovSupportPlate", 1);
  CerenkovBoxInnerLogical.placeVolume(CerenkovElecLogical,
      Position(0., straw_y + csh/2. + cst + cebt/2., 0.)).addPhysVolID("CerenkovElecBoard", 0);

  // -----------------------------------------------------------------------
  // Assemble CerenkovMotherLogical
  // -----------------------------------------------------------------------
  fCerenkovMotherVol = Volume(kCherenkovVolName, solidCerenkovMother, envMaterial);

  fCerenkovMotherVol.placeVolume(CerenkovBoxLogical,
      Position(0.,0.,0.)); //.addPhysVolID("CerenkovBox", 0);
  fCerenkovMotherVol.placeVolume(CerenkovBoxInnerLogical,
      Position(0.,0.,0.)); //.addPhysVolID("CerenkovBoxInner", 0);

  // Shielding plates
  double shieldZ = fHalfZ + cbt + cspt/2.;
  double shieldUpperY = CerenkovTotalBoxHeight/2.
                        - (CerenkovTotalBoxHeight - ComptonElectronBeamtoStageY - CerenkovBeamWindowY/2.)/2.;
  double shieldLowerY = -CerenkovTotalBoxHeight/2.
                        + (ComptonElectronBeamtoStageY - CerenkovBeamWindowY/2.)/2.;

  fCerenkovMotherVol.placeVolume(ShieldUpperLogical, Position(0., shieldUpperY,  shieldZ)).addPhysVolID("ShieldUpper", 0);
  fCerenkovMotherVol.placeVolume(ShieldUpperLogical, Position(0., shieldUpperY, -shieldZ)).addPhysVolID("ShieldUpper", 1);
  fCerenkovMotherVol.placeVolume(ShieldLowerLogical, Position(0., shieldLowerY,  shieldZ)).addPhysVolID("ShieldLower", 0);
  fCerenkovMotherVol.placeVolume(ShieldLowerLogical, Position(0., shieldLowerY, -shieldZ)).addPhysVolID("ShieldLower", 1);

  double shieldSideX = fHalfX + cbt + cspt/2.;
  fCerenkovMotherVol.placeVolume(ShieldSideLogical, Position( shieldSideX, 0., 0.)).addPhysVolID("ShieldSide", 0);
  fCerenkovMotherVol.placeVolume(ShieldSideLogical, Position(-shieldSideX, 0., 0.)).addPhysVolID("ShieldSide", 1);

  printout(INFO, "LxCherenkov::Build",
           "Built CerenkovMotherLogical: halfX=%.1f halfZ=%.1f height=%.1f mm, %d straws.",
           fHalfX/dd4hep::mm, fHalfZ/dd4hep::mm, CerenkovTotalBoxHeight/dd4hep::mm, fNStraws);


// Create position map, easy for translation from G4
// Positioning constants fo HICSECherenkov
  double IPMagFieldY             = description.constant<double>("IPMagFieldY");
  double IPMagnetZpos            = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength       = description.constant<double>("FlashMFieldLength");
  double ComptonElBackshift      = description.constant<double>("ComptonElBackshift");
  double ScintFrameThickness     = description.constant<double>("ScintFrameThickness");
  double HICSDetBottomSupportZ   = description.constant<double>("HICSDetBottomSupportZ");
  double HICSCherenkovXPosMag2T  = description.constant<double>("HICSCherenkovXPosMag2T");

  double cerx = ((cc/ccl - 1)*csxf/2. + (ccl-1)*cslo +  sir/2. + (sglt+salt+sklt+splt) + cspb + cbt + cspt);
  double cerz =  (ccl-1)*cszf/2. +  sir/2. + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt;
  double comptonelzpos = IPMagnetZpos + FlashMFieldLength/2.0  + ComptonElBackshift;
  double shift_z = comptonelzpos - ComptonElectronBeamtoStageY - 3.0*ScintFrameThickness/2.0 + HICSDetBottomSupportZ;
  double cermxpos = 0.0;

  if (std::abs(IPMagFieldY + 10000.0*gauss) < 1.0e-10) {
    cermxpos = description.constant<double>("HICSCherenkovXPosMag1T");
  } else if (std::abs(IPMagFieldY + 20000.0*gauss) < 1.0e-10) {
    cermxpos = description.constant<double>("HICSCherenkovXPosMag2T");
  } else {
    printout(FATAL, "LxCherenkov", "Position for the %f is not supported\n", IPMagFieldY);
  }

  double hicsxpos = -(cerx+cermxpos);
  double hicsypos = -ComptonElectronBeamtoStageY + CerenkovTotalBoxHeight/2.0;
  double hicszpos = shift_z - cerz;
  fDetTransformMap.emplace("HICSElectronCerenkov", Transform3D(Position(hicsxpos, hicsypos, hicszpos)));

// Positioning constants for GammaSpectrometerCerenkov
  double ComptonLysoZpos      = description.constant<double>("ComptonLysoZpos");
  double ComptonLysoZ         = description.constant<double>("ComptonLysoZ");
  double ComptonLysoXpos      = description.constant<double>("ComptonLysoXpos");
  double ComptonLysoX         = description.constant<double>("ComptonLysoX");

  double lysozpos = ComptonLysoZpos + ComptonLysoZ/2.0 + 100.0*mm;

  double comptxpos = ComptonLysoXpos - ComptonLysoX/2.0 + cerx;
  double comptypos = -ComptonElectronBeamtoStageY + CerenkovTotalBoxHeight/2.0;
  double comptzpos = lysozpos + cerz + 150.0*mm;
  fDetTransformMap.emplace("GammaSpectrometerCerenkov_L", Transform3D(Position(comptxpos, comptypos, comptzpos)));
  fDetTransformMap.emplace("GammaSpectrometerCerenkov_R", Transform3D(Position(-comptxpos, comptypos, comptzpos)));

//Positioning constants for BremCerenkov
  double CerenkovXpos         = description.constant<double>("CerenkovXpos");
  double CerenkovZpos         = description.constant<double>("CerenkovZpos");
  double CerenkovZrot         = description.constant<double>("CerenkovAngle");
  double bremscxpos = ComptonElectronBeamtoStageY - CerenkovTotalBoxHeight/2.0;
  fDetTransformMap.emplace("BremCerenkov",
                           Transform3D(RotationZYX(-0.5*M_PI, CerenkovZrot, 0.0), Position(bremscxpos, CerenkovXpos, CerenkovZpos)));
}


// ---------------------------------------------------------------------------
dd4hep::DetElement LxCherenkov::Place(dd4hep::Detector&          description,
                                       dd4hep::Volume&            motherVol,
                                       const std::string&         detName,
                                       int                        detID
                                       /*, const dd4hep::Transform3D& trf*/)
{
  if (!fBuilt)
    dd4hep::except("LxCherenkov::Place", "Build() must be called before Place()");

  bool OverlapTest = (description.constant<int>("OverlapTest") != 0);

  auto itr = fDetTransformMap.find(detName);
  if (itr == fDetTransformMap.end()) {
    PrintCherenkovTranslations();
    printout(FATAL, "LxCherenkov", "Position for the Cherenkov detector with name %s is not defined\n", detName.c_str());
    throw std::runtime_error("Cannot place detector");
  }

  Transform3D &trf = itr->second;

  // Place CerenkovMotherLogical into mother volume
  PlacedVolume pv = motherVol.placeVolume(fCerenkovMotherVol, trf);
  pv.addPhysVolID("detector", detID);
  if (OverlapTest) pv.ptr()->CheckOverlaps();

  // DetElement for this placement
  DetElement det(detName, detID);
  det.setPlacement(pv);

  // Clone straw DetElement for each straw in this placement
  for (int s = 0; s < fNStraws; ++s) {
    DetElement strawDE = fStrawTemplate.clone("straw_" + std::to_string(s), s);
    det.add(strawDE);
  }

  return det;
}


void LxCherenkov::PrintCherenkovTranslations()
{
  std::stringstream sstr("Supported detector names:", std::ios_base::ate | std::ios_base::in | std::ios_base::out);
  std::for_each(fDetTransformMap.begin(), fDetTransformMap.end(), [&](const auto pp){sstr << "  " << pp.first;});
  printout(INFO, "LxCherenkov", sstr);
}



// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// Called once per <detector type="LxCherenkov"> XML entry.
// First call builds and registers the LxCherenkov instance.
// All calls place the detector and return a DetElement.
// ---------------------------------------------------------------------------
static Ref_t create_LxCherenkov(dd4hep::Detector& description,
                                 xml_h             e,
                                 dd4hep::SensitiveDetector sd)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  sd.setType("calorimeter");
  // Retrieve or create the shared LxCherenkov instance
  LxCherenkov* cer = nullptr;
  try {
    cer = description.extension<LxCherenkov>();
    printout(INFO, "LxCherenkov", "Reusing existing LxCherenkov instance for '%s'.", detName.c_str());
  } catch (...) {
    printout(INFO, "LxCherenkov", "Creating new LxCherenkov instance for '%s'.", detName.c_str());
    cer = new LxCherenkov();
    cer->Build(description, sd);
    description.addExtension<LxCherenkov>(cer);
  }

  // Each call gets its own envelope assembly and DetElement
  Volume       fLogicWorld = description.worldVolume();
  Assembly     envelope(detName + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.attr<int>(xml_tag_t("system_id")));

  dd4hep::DetElement sdet(detName, detID);
  sdet.setPlacement(envPV);

  // Read placement transform from XML child element <placement .../>
//   xml_comp_t x_pl = x_det.child(xml_tag_t("placement"));
//   double xpos = x_pl.attr<double>(xml_tag_t("xpos"));
//   double ypos = x_pl.attr<double>(xml_tag_t("ypos"));
//   double zpos = x_pl.attr<double>(xml_tag_t("zpos"));
//   double rotZ = x_pl.hasAttr(xml_tag_t("rotZ")) ? x_pl.attr<double>(xml_tag_t("rotZ")) : 0.0;
//   double rotY = x_pl.hasAttr(xml_tag_t("rotY")) ? x_pl.attr<double>(xml_tag_t("rotY")) : 0.0;

//   Transform3D trf(RotationZYX(rotZ, rotY, 0.0), Position(xpos, ypos, zpos));

  // Place and get DetElement for this copy
  dd4hep::DetElement cerDE = cer->Place(description, envelope, detName, detID /*, trf*/);
  sdet.add(cerDE);

  return sdet;
}

DECLARE_DETELEMENT(LxCherenkov, create_LxCherenkov)
