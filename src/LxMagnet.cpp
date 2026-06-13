//
/// \brief Generic DD4hep plugin for placing any LxMagnetAssembly model.
///
/// Usage in XML:
///   <detector id="N" name="DumpMagnet" type="LxMagnet">
///     <magnet model="FlashMagnet"
///             xpos="0.0*mm" ypos="0.0*mm" zpos="1000.0*mm"
///             rotate="true"
///             fieldvolname="logicDumpMagnetField"
///             support="true"/>
///   </detector>
///
/// The plugin:
///   1. Reads placement parameters and magnet model name from XML.
///   2. Checks TGeoManager for an existing assembly of that model.
///   3. If absent: instantiates the appropriate model class and builds it.
///      If present: wraps the cached TGeoVolume in a dd4hep::Assembly.
///   4. Places the assembly and a fresh field volume at the given position.
///   5. Optionally places support structure.
//

#include <stdexcept>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TGeoManager.h"

#include "LxMagnetAssembly.h"
#include "FlashMagnetAssembly.h"
// #include "TypMBMagnetAssembly_dd4hep.hh"   // include as more models are added

using namespace dd4hep;

void AddDumpMagnetFieldGeometry(dd4hep::Detector& description, dd4hep::Volume fieldVol);
void AddIPMagnetFieldGeometry(dd4hep::Detector& description, dd4hep::Volume fieldVol);
void AddGammaMagnetFieldGeometry(dd4hep::Detector& description, dd4hep::Volume fieldVol);


static Ref_t create_LxMagnet(dd4hep::Detector& description,
                              xml_h             e,
                              dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  // Envelope assembly — satisfies DetElement placement requirement
  Volume       fLogicWorld = description.worldVolume();
  Assembly     envelope(detName + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", detID);
  sdet.setPlacement(envPV);

  bool OverlapTest = (description.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Read placement parameters from <magnet ...> child element
  // -----------------------------------------------------------------------
  xml_comp_t x_magnet = x_det.child(xml_tag_t("magnet"));

  std::string model        = x_magnet.attr<std::string>(xml_tag_t("model"));
  std::string fieldvolname = x_magnet.attr<std::string>(xml_tag_t("fieldvolname"));
  bool        rotate       = x_magnet.hasAttr(xml_tag_t("rotate"))
                             ? x_magnet.attr<bool>(xml_tag_t("rotate")) : false;
  bool        support      = x_magnet.hasAttr(xml_tag_t("support"))
                             ? x_magnet.attr<bool>(xml_tag_t("support")) : false;

  double xpos = x_magnet.attr<double>(xml_tag_t("xpos"));
  double ypos = x_magnet.attr<double>(xml_tag_t("ypos"));
  double zpos = x_magnet.attr<double>(xml_tag_t("zpos"));
  Position pos(xpos, ypos, zpos);

  // Placement rotation: the dump magnet is placed rotated pi/2 around Z.
  // G4RotationMatrix(G4ThreeVector(0,0,1), pi/2) → RotationZYX(pi/2, 0, 0)
  // Generalise: always apply pi/2 around Z when rotate==true, identity otherwise.
  RotationZYX rot = rotate ? RotationZYX(M_PI/2.0, 0.0, 0.0)
                           : RotationZYX(0.0, 0.0, 0.0);

  // -----------------------------------------------------------------------
  // Retrieve or build the magnet assembly
  // -----------------------------------------------------------------------
  const std::string asmName = model + "_assembly";
  TGeoVolume*       existing = description.manager().GetVolume(asmName.c_str());

  LxMagnetAssembly* mag = nullptr;
  Assembly          magAssembly("");

  if (!existing) {
    // First placement of this model — instantiate and build
    if      (model == "FlashMagnet") mag = new FlashMagnetAssembly();
    // else if (model == "TypMBMagnet") mag = new TypMBMagnetAssembly();
    else {
      throw std::runtime_error("LxMagnet_geo: unknown magnet model '" + model + "'");
    }
    magAssembly = mag->GetAssembly(description);
    printout(INFO, "LxMagnet_geo",
             "Built new assembly '%s' for detector '%s'.",
             asmName.c_str(), detName.c_str());
  } else {
    // Subsequent placement — reuse cached assembly from TGeoManager
    printout(INFO, "LxMagnet_geo",
             "Reusing cached assembly '%s' for detector '%s'.",
             asmName.c_str(), detName.c_str());
    magAssembly = Assembly(Volume(existing));
  }

  // -----------------------------------------------------------------------
  // Place magnet assembly
  // -----------------------------------------------------------------------
  PlacedVolume pvMag = envelope.placeVolume(magAssembly, Transform3D(rot, pos));
  pvMag.addPhysVolID(model, 0);
  if (OverlapTest) pvMag.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Place field volume — always fresh, unique name per placement
  // Static method: no mag instance needed when reusing cached assembly
  // -----------------------------------------------------------------------
  Volume fieldVol;
  if (model == "FlashMagnet") {
    fieldVol = FlashMagnetAssembly::GetFieldVolumeStatic(description, fieldvolname);
    if (detName == "DumpMagnet") AddDumpMagnetFieldGeometry(description, fieldVol);
    if (detName == "IPMagnet") AddIPMagnetFieldGeometry(description, fieldVol);
    if (detName == "GMagnet") AddGammaMagnetFieldGeometry(description, fieldVol);
  }
  // else if (model == "TypMBMagnet") {
  //   fieldVol = TypMBMagnetAssembly::GetFieldVolume(description, fieldvolname);
  // }
  else {
    throw std::runtime_error("LxMagnet_geo: no GetFieldVolume for model '" + model + "'");
  }

  // Field volume uses opposite Z rotation sign (matches original ConstructDumpMagnet)
  // That was claude's vision, not clear why. Fixed to same Z angle sign as for the magnet
  RotationZYX fieldRot = rotate ? RotationZYX(M_PI/2.0, 0.0, 0.0)
                                : RotationZYX(0.0, 0.0, 0.0);
  PlacedVolume pvField = envelope.placeVolume(fieldVol, Transform3D(fieldRot, pos));
  pvField.addPhysVolID(fieldvolname, 0);
  if (OverlapTest) pvField.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Support structure (optional)
  // -----------------------------------------------------------------------
  if (support) {
    if (mag) {
      // Have instance — call virtual method directly
      mag->ConstructSupport(description, envelope, pos, detName, rotate);
    } else {
      // Reusing cached assembly — need a temporary instance for support
      // Support geometry does not get cached so instantiation is cheap
      if (model == "FlashMagnet") {
        FlashMagnetAssembly tmp;
        tmp.ConstructSupport(description, envelope, pos, detName, rotate);
      }
      // else if (model == "TypMBMagnet") { ... }
    }
  }

  // Clean up temporary instance if we created one
  delete mag;

  return sdet;
}

DECLARE_DETELEMENT(LxMagnet, create_LxMagnet)




// ---------------------------------------------------------------------------
void AddDumpMagnetFieldGeometry(dd4hep::Detector& description, dd4hep::Volume fieldVol)
{
  double BPipeR         = description.constant<double>("BPipeR");
  double BPipeThickness = description.constant<double>("BPipeThickness");
  double DumpMagnetYpos = description.constant<double>("DumpMagnetYpos");

  Material beamPipeMaterial = description.material(description.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial   = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));

  bool OverlapTest = (description.constant<int>("OverlapTest") != 0);

  // Field volume half-length along Z (needed for beam pipe sizing)
  // The field volume is a Box — extract hz from the FlashMFieldLength constant
  double fhlength = description.constant<double>("FlashMFieldLength") / 2.0;

  // -----------------------------------------------------------------------
  // Beam pipe inside field volume (rectangular)
  // -----------------------------------------------------------------------
  double ipmbphx = 2.0*(BPipeR + BPipeThickness);

  Box solidDMBPipeContainer(ipmbphx, BPipeR, fhlength);
  Volume logicDMBPipeContainer("logicDMBPipeContainer", solidDMBPipeContainer, vacuumMaterial);

  Box solidDMBPipeOuter(ipmbphx, BPipeR, fhlength);
  Box solidDMBPipeInner(2.0*BPipeR + BPipeThickness, (BPipeR - BPipeThickness),
                        fhlength - BPipeThickness);
  SubtractionSolid solidDMBPipeXY("solidDMBPipeXY", solidDMBPipeOuter, solidDMBPipeInner);

  // Round hole for incoming pipe (-Z face)
  Tube solidSHoleIn(0.0, BPipeR - BPipeThickness, BPipeThickness, 0.0, 2.0*M_PI);
  SubtractionSolid solidDMBPipe1("solidDMBPipe1", solidDMBPipeXY, solidSHoleIn,
      // G4RotationMatrix(Z, 0): identity rotation
      Position(BPipeR + 2.0*BPipeThickness, 0.0, -fhlength + BPipeThickness/2.0));

  // Rectangular hole for outgoing pipe (+Z face)
  Box solidoutcut(2.0*BPipeR + BPipeThickness, (BPipeR - BPipeThickness), BPipeThickness);
  SubtractionSolid solidDMBPipe("solidDMBPipe", solidDMBPipe1, solidoutcut,
      Position(0.0, 0.0, fhlength - BPipeThickness/2.0));

  Volume logicDMBPipe("logicDMBPipe", solidDMBPipe, beamPipeMaterial);

  PlacedVolume pvPipe = logicDMBPipeContainer.placeVolume(logicDMBPipe, Position(0.0, 0.0, 0.0));
  pvPipe.addPhysVolID("DMBPipe", 0);
  if (OverlapTest) pvPipe.ptr()->CheckOverlaps();

  PlacedVolume pvCont = fieldVol.placeVolume(logicDMBPipeContainer,
      Position(-(BPipeR + 2.0*BPipeThickness) - DumpMagnetYpos, 0.0, 0.0));
  pvCont.addPhysVolID("DMBPipeContainer", 0);
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();
}



// ---------------------------------------------------------------------------
void AddIPMagnetFieldGeometry(dd4hep::Detector& description, dd4hep::Volume fieldVol)
{
  Material beamPipeMaterial = description.material(
      description.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial = description.material(
      description.constant<std::string>("BeamPipeVacuumMaterial"));

  double BPipeThickness       = description.constant<double>("BPipeThickness");
  double IPMAgnetBeamPipeXGap = description.constant<double>("IPMAgnetBeamPipeXGap");
  double fhlength = description.constant<double>("FlashMFieldLength") / 2.0;
  double fhwidth  = description.constant<double>("FlashMFieldX")      / 2.0;
  double fhhight  = description.constant<double>("FlashMagneteffY")   / 2.0;
  bool   OverlapTest = (description.constant<int>("OverlapTest") != 0);

  double ipmbphx = fhwidth - IPMAgnetBeamPipeXGap;

  // Container
  Box    solidIPMBPipeContainer(ipmbphx, fhhight, fhlength);
  Volume logicIPMBPipeContainer("logicIPMBPipeContainer", solidIPMBPipeContainer, vacuumMaterial);

  // Pipe wall
  Box solidIPMBPipeOuter(ipmbphx, fhhight, fhlength);
  Box solidIPMBPipeInner(ipmbphx - BPipeThickness,
                         fhhight - BPipeThickness,
                         fhlength - BPipeThickness);
  SubtractionSolid solidIPMBPipeXY("solidIPMBPipeXY", solidIPMBPipeOuter, solidIPMBPipeInner);

  // Rectangular inlet hole on -Z face
  Box solidHoleIn(ipmbphx - BPipeThickness, fhhight - BPipeThickness, BPipeThickness);
  SubtractionSolid solidIPMBPipe1("solidIPMBPipe1", solidIPMBPipeXY, solidHoleIn,
      Position(0.0, 0.0, -fhlength + BPipeThickness/2.0));

  // Rectangular outlet hole on +Z face
  Box solidOutCut(ipmbphx - BPipeThickness, fhhight - BPipeThickness, BPipeThickness);
  SubtractionSolid solidIPMBPipe("solidIPMBPipe", solidIPMBPipe1, solidOutCut,
      Position(0.0, 0.0, fhlength - BPipeThickness/2.0));

  Volume logicIPMBPipe("logicIPMBPipe", solidIPMBPipe, beamPipeMaterial);

  PlacedVolume pvPipe = logicIPMBPipeContainer.placeVolume(logicIPMBPipe,
                                                           Position(0.0, 0.0, 0.0));
  pvPipe.addPhysVolID("IPMBPipe", 0);
  if (OverlapTest) pvPipe.ptr()->CheckOverlaps();

  PlacedVolume pvCont = fieldVol.placeVolume(logicIPMBPipeContainer,
                                             Position(0.0, 0.0, 0.0));
  pvCont.addPhysVolID("IPMBPipeContainer", 0);
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();
}



// ---------------------------------------------------------------------------
void AddGammaMagnetFieldGeometry(dd4hep::Detector& description, dd4hep::Volume fieldVol)
{
  using namespace dd4hep;

  Material beamPipeMaterial  = description.material(description.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial    = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));
  Material bpipeWindowMaterial = description.material(description.constant<std::string>("BeamPipeWindowMaterial"));

  double gmpipex    = description.constant<double>("QBeamPipeContainerX");
  double gmpiptby   = description.constant<double>("GChamberTBWallThickness");
  double gmpipsidex = description.constant<double>("GammaBPipeWindowThickness");
  double BPipeRLG   = description.constant<double>("BPipeRLG");
  double BPipeThickness = description.constant<double>("BPipeThickness");
  double GammaMagnetBeamPipeXGap = description.constant<double>("GammaMagnetBeamPipeXGap");
  bool   OverlapTest = (description.constant<int>("OverlapTest") != 0);

  // Field volume half-dimensions — read from constants, same source as GetFieldVolume
  double gfieldz = description.constant<double>("FlashMFieldLength") / 2.0;
  double gfieldy = description.constant<double>("FlashMagneteffY")   / 2.0;
  double gfieldx = description.constant<double>("FlashMFieldX")      / 2.0;

  double swy = gfieldy - gmpiptby;  // height of the pipe / chamber interior

  // -----------------------------------------------------------------------
  // Narrow beam pipe container (beam axis side)
  // -----------------------------------------------------------------------
  Box    solidGammaMagnetPipeContainer(gmpipex/2.0, gfieldy, gfieldz);
  Volume logicGammaMagnetPipeContainer("logicGammaMagnetPipeContainer",
                                        solidGammaMagnetPipeContainer, vacuumMaterial);

  // Top / bottom walls
  Box    solidGamMagPipeTopBot(gmpipex/2.0, gmpiptby/2.0, gfieldz);
  Volume logicGamMagPipeTopBot("logicGamMagPipeTopBot", solidGamMagPipeTopBot, beamPipeMaterial);

  // Side walls (commented out in original — preserved)
  Box    solidGamMagPipeSide(gmpipsidex/2.0, swy, gfieldz);
  Volume logicGamMagPipeSide("logicGamMagPipeSide", solidGamMagPipeSide, beamPipeMaterial);

  // Front window (commented out in original — preserved)
  Box    solidGamMagPipeFrontWind(gmpipex/2.0, swy, gmpipsidex/2.0);
  Volume logicGamMagPipeFrontWind("logicGamMagPipeFrontWind", solidGamMagPipeFrontWind, bpipeWindowMaterial);

  // Rear wall with round beam pipe cutout
  Box  solidGamMagPipeRearBox(gmpipex/2.0, swy, gmpiptby/2.0);
  double incutr = BPipeRLG - BPipeThickness;
  if (incutr >= swy) incutr = 0.96 * swy;
  Tube solidGamMagPipeRearBeamPipeCut(0.0, incutr, gmpiptby, 0.0, 2.0*M_PI);
  SubtractionSolid solidGamMagPipeRear("solidGamMagPipeRear",
                                        solidGamMagPipeRearBox, solidGamMagPipeRearBeamPipeCut);
  Volume logicGamMagPipeRear("logicGamMagPipeRear", solidGamMagPipeRear, beamPipeMaterial);

  // Placements into narrow container
  PlacedVolume pvTop = logicGammaMagnetPipeContainer.placeVolume(logicGamMagPipeTopBot,
      Position(0.0, gfieldy - gmpiptby/2.0, 0.0));
  pvTop.addPhysVolID("GMagnetFieldPipeTop", 0);
  if (OverlapTest) pvTop.ptr()->CheckOverlaps();

  PlacedVolume pvBot = logicGammaMagnetPipeContainer.placeVolume(logicGamMagPipeTopBot,
      Position(0.0, -(gfieldy - gmpiptby/2.0), 0.0));
  pvBot.addPhysVolID("GMagnetFieldPipeBottom", 1);
  if (OverlapTest) pvBot.ptr()->CheckOverlaps();

//   logicGammaMagnetPipeContainer.placeVolume(logicGamMagPipeSide,
//       Position((gmpipex - gmpipsidex)/2.0, 0.0, 0.0)).addPhysVolID("GMagnetFieldPipeSide", 0);
//   logicGammaMagnetPipeContainer.placeVolume(logicGamMagPipeSide,
//       Position(-(gmpipex - gmpipsidex)/2.0, 0.0, 0.0)).addPhysVolID("GMagnetFieldPipeSide", 1);
//   logicGammaMagnetPipeContainer.placeVolume(logicGamMagPipeFrontWind,
//       Position(0.0, 0.0, gfieldz - gmpipsidex/2.0)).addPhysVolID("GamMagPipeFrontWind", 1);

  PlacedVolume pvRear = logicGammaMagnetPipeContainer.placeVolume(logicGamMagPipeRear,
      Position(0.0, 0.0, -(gfieldz - gmpiptby/2.0)));
  pvRear.addPhysVolID("GamMagPipeRear", 1);
  if (OverlapTest) pvRear.ptr()->CheckOverlaps();

  // Place narrow container into field volume
  PlacedVolume pvNarrow = fieldVol.placeVolume(logicGammaMagnetPipeContainer,
      Position(0.0, 0.0, 0.0));
  pvNarrow.addPhysVolID("GMagnetFieldBeamPipe", 0);
  if (OverlapTest) pvNarrow.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Wide beam pipe container (electron/positron side)
  // -----------------------------------------------------------------------
  double wpipecontainerx = gfieldx - gmpipex/2.0 - GammaMagnetBeamPipeXGap;

  Box    solidGammaMagnetWidePipeContainer(wpipecontainerx/2.0, gfieldy, gfieldz);
  Volume logicGammaMagnetWidePipeContainer("logicGammaMagnetWidePipeContainer",
                                            solidGammaMagnetWidePipeContainer, vacuumMaterial);

  // Top / bottom walls
  Box    solidGamMagWidePipeTopBot(wpipecontainerx/2.0, gmpiptby/2.0, gfieldz);
  Volume logicGamMagWidePipeTopBot("logicGamMagWidePipeTopBot", solidGamMagWidePipeTopBot, beamPipeMaterial);

  // Side wall
  Box    solidGamMagWidePipeSide(gmpipsidex/2.0, gfieldy - gmpiptby, gfieldz);
  Volume logicGamMagWidePipeSide("logicGamMagWidePipeSide", solidGamMagWidePipeSide, beamPipeMaterial);

  // Front window (commented out in original — preserved)
  Box    solidGamMagWidePipeFront(wpipecontainerx/2.0 - gmpipsidex/2.0,
                                  gfieldy - gmpiptby, gmpipsidex/2.0);
  Volume logicGamMagWidePipeFront("logicGamMagWidePipeFront", solidGamMagWidePipeFront, bpipeWindowMaterial);

  // Rear wall
  Box    solidGamMagWidePipeRear(wpipecontainerx/2.0 - gmpipsidex/2.0,
                                 gfieldy - gmpiptby, gmpiptby/2.0);
  Volume logicGamMagWidePipeRear("logicGamMagWidePipeRear", solidGamMagWidePipeRear, beamPipeMaterial);

  // Placements into wide container
//   double GMagnetZpos = description.constant<double>("GMagnetZpos");
//   double GMagnetCutX = description.constant<double>("GMagnetCutX");
//   logicGammaMagnetWidePipeContainer.placeVolume(logicGamMagPipeSide,
//       Position((gmpipex - gmpipsidex)/2.0, 0.0, GMagnetZpos)).addPhysVolID("GMagnetWidePipeSide", 0);
//   logicGammaMagnetWidePipeContainer.placeVolume(logicGamMagPipeSide,
//       Position(-(GMagnetCutX)/2.0, 0.0, GMagnetZpos)).addPhysVolID("GMagnetWidePipeSide", 1);

  PlacedVolume pvWTop = logicGammaMagnetWidePipeContainer.placeVolume(logicGamMagWidePipeTopBot,
      Position(0.0, gfieldy - gmpiptby/2.0, 0.0));
  pvWTop.addPhysVolID("GMagnetFieldWidePipeTop", 0);
  if (OverlapTest) pvWTop.ptr()->CheckOverlaps();

  PlacedVolume pvWBot = logicGammaMagnetWidePipeContainer.placeVolume(logicGamMagWidePipeTopBot,
      Position(0.0, -(gfieldy - gmpiptby/2.0), 0.0));
  pvWBot.addPhysVolID("GMagnetFieldWidePipeBottom", 1);
  if (OverlapTest) pvWBot.ptr()->CheckOverlaps();

  PlacedVolume pvWSide = logicGammaMagnetWidePipeContainer.placeVolume(logicGamMagWidePipeSide,
      Position((wpipecontainerx - gmpipsidex)/2.0, 0.0, 0.0));
  pvWSide.addPhysVolID("GMagnetFieldWidePipeSide", 0);
  if (OverlapTest) pvWSide.ptr()->CheckOverlaps();

//   logicGammaMagnetWidePipeContainer.placeVolume(logicGamMagWidePipeFront,
//       Position(-gmpipsidex/2.0, 0.0, gfieldz - gmpipsidex/2.0)).addPhysVolID("GamMagWidePipeFront", 0);

  PlacedVolume pvWRear = logicGammaMagnetWidePipeContainer.placeVolume(logicGamMagWidePipeRear,
      Position(-gmpipsidex/2.0, 0.0, -(gfieldz - gmpiptby/2.0)));
  pvWRear.addPhysVolID("GamMagWidePipeRear", 0);
  if (OverlapTest) pvWRear.ptr()->CheckOverlaps();

  // Place wide container into field volume — left side
  PlacedVolume pvWideL = fieldVol.placeVolume(logicGammaMagnetWidePipeContainer,
      Position((gmpipex + wpipecontainerx)/2.0, 0.0, 0.0));
  pvWideL.addPhysVolID("GMagnetWideBeamPipeL", 0);
  if (OverlapTest) pvWideL.ptr()->CheckOverlaps();

  // Place wide container into field volume — right side, rotated pi around Z
  // G4RotationMatrix(G4ThreeVector(0,0,1), pi): axis-angle around Z → RotationZYX(pi, 0, 0)
  PlacedVolume pvWideR = fieldVol.placeVolume(logicGammaMagnetWidePipeContainer,
      Transform3D(RotationZYX(M_PI, 0.0, 0.0),
                  Position(-(gmpipex + wpipecontainerx)/2.0, 0.0, 0.0)));
  pvWideR.addPhysVolID("GMagnetWideBeamPipeR", 1);
  if (OverlapTest) pvWideR.ptr()->CheckOverlaps();
}

