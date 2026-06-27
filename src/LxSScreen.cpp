//
/// \brief Implementation of the LxSScreen class (DD4hep version)
//

#include <cmath>
#include <string>
#include <algorithm>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxAux.h"
#include "LxSScreen.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
void LxSScreen::Build(dd4hep::Detector& description, dd4hep::SensitiveDetector& sd)
{
  if (fBuilt) return;
  fBuilt = true;

  BuildScintArm(description, sd);
  BuildCamera(description);

  // Register placement methods in dispatch map
  fPlaceMap["LxHICSScint"]  = [this](dd4hep::Detector& d, dd4hep::Volume& v, int id)
                               { return PlaceHICSScint(d, v, id); };
  fPlaceMap["LxBremsScint"] = [this](dd4hep::Detector& d, dd4hep::Volume& v, int id)
                               { return PlaceBremsScint(d, v, id); };
}


// ---------------------------------------------------------------------------
void LxSScreen::BuildScintArm(dd4hep::Detector& description, dd4hep::SensitiveDetector& sd)
{
  Material env_mat    = description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material frame_mat  = description.material(description.constant<std::string>("ScintFrameMaterial"));
  Material base_mat   = description.material(description.constant<std::string>("ScintBaseMaterial"));
  Material scint_mat  = description.material(description.constant<std::string>("ScintPhosphorMaterial"));
  Material finish_mat = description.material(description.constant<std::string>("ScintFinishMaterial"));

  double sft      = description.constant<double>("ScintFrameThickness");
  double scx      = description.constant<double>("ScintX");
  double scy      = description.constant<double>("ScintY");
  double scz1     = description.constant<double>("ScintBaseZ");
  double scz2     = description.constant<double>("ScintPhosphorZ");
  double scz3     = description.constant<double>("ScintFinishZ");
  double beam2stage = description.constant<double>("ComptonElectronBeamtoStageY");
  double sca      = description.constant<double>("ScintAngle");
  double bsfwo    = description.constant<double>("BremScintFrameWidthOffset");
  double ScintFrameBeamLoopX = description.constant<double>("ScintFrameBeamLoopX");
  bool   OverlapTest = (description.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Mother volumes
  // -----------------------------------------------------------------------
  Box solidScintArmMother(scx/2. + sft + ScintFrameBeamLoopX, beam2stage, beam2stage + 3*sft/2.);
  Box solidBremScintArmMother(scy/2. + 2*sft + bsfwo,
                              std::cos(sca)*scx/2. + std::cos(sca)*sft + std::sin(sca)*sft/2.,
                              std::sin(sca)*scx/2. + std::sin(sca)*sft + std::cos(sca)*sft/2.);

  fScintArmMotherVol     = Volume("scintArmMotherLogical",    solidScintArmMother,    env_mat);
  fBremScintArmMotherVol = Volume("BremScintArmMotherLogical",solidBremScintArmMother,env_mat);

  // -----------------------------------------------------------------------
  // Scint arm inner structure (shared between both mothers)
  // -----------------------------------------------------------------------
  Box solidScintArm(scx/2., scy/2., (scz1+scz2+scz3)/2.);
  fScintArmVol = Volume("scintArmLogical", solidScintArm, env_mat);

  Box    solidScintBase    (scx/2., scy/2., scz1/2.);
  Box    solidScintPhosphor(scx/2., scy/2., scz2/2.);
  Box    solidScintFinish  (scx/2., scy/2., scz3/2.);

  Volume scintBaseLogical    ("scintBaseLogical",    solidScintBase,     base_mat);
  Volume scintPhosphorLogical("scintPhosphorLogical",solidScintPhosphor, scint_mat);
  Volume scintFinishLogical  ("scintFinishLogical",  solidScintFinish,   finish_mat);

  // Mark phosphor as sensitive
  scintPhosphorLogical.setSensitiveDetector(sd);

  fScintArmVol.placeVolume(scintBaseLogical,
      Position(0., 0., -(scz2+scz3)/2.)).addPhysVolID("scintBase", 0);
  PlacedVolume pvPhosphor = fScintArmVol.placeVolume(scintPhosphorLogical,
      Position(0., 0., (scz1-scz3)/2.));
  // No addPhysVolID — detector ID comes from placement of scintArmMotherVol
  if (OverlapTest) pvPhosphor.ptr()->CheckOverlaps();
  fScintArmVol.placeVolume(scintFinishLogical,
      Position(0., 0., (scz1+scz2)/2.)).addPhysVolID("scintFinish", 0);

  // Phosphor DetElement template — cloned per placement
  fPhosphorTemplate = DetElement("phosphor", 0);
  fPhosphorTemplate.setPlacement(pvPhosphor);

  // -----------------------------------------------------------------------
  // Frame legs and bars for scintArmMotherLogical
  // -----------------------------------------------------------------------
  Box solidScintFrameLeg(sft/2., scy/2. + sft, sft/2.);
  Box solidScintFrameBar(scx/2. + ScintFrameBeamLoopX/2., sft/2., sft/2.);

  Volume scintFrameLegLogical("scintFrameLegLogical", solidScintFrameLeg, frame_mat);
  Volume scintFrameBarLogical("scintFrameBarLogical", solidScintFrameBar, frame_mat);

  // Place scintArm into scintArmMother
  fScintArmMotherVol.placeVolume(fScintArmVol,
      Position(0., 0., 0.));
//cl       Position(0., 0., 0.)).addPhysVolID("scintArm", 0);

  // Frame legs
  fScintArmMotherVol.placeVolume(scintFrameLegLogical,
      Position(-(scx+sft)/2., 0., 0.)).addPhysVolID("scintFrameLeg", 0);
  fScintArmMotherVol.placeVolume(scintFrameLegLogical,
      Position((scx+sft)/2.+ScintFrameBeamLoopX, 0., 0.)).addPhysVolID("scintFrameLeg", 1);

  // Frame bars
  fScintArmMotherVol.placeVolume(scintFrameBarLogical,
      Position(ScintFrameBeamLoopX/2.,  (scy+sft)/2., 0.)).addPhysVolID("scintFrameBar", 0);
  fScintArmMotherVol.placeVolume(scintFrameBarLogical,
      Position(ScintFrameBeamLoopX/2., -(scy+sft)/2., 0.)).addPhysVolID("scintFrameBar", 1);
  fScintArmMotherVol.placeVolume(scintFrameBarLogical,
      Position(ScintFrameBeamLoopX/2., -(beam2stage)+sft/2.,  beam2stage-sft/2.)).addPhysVolID("scintFrameBar", 2);
  fScintArmMotherVol.placeVolume(scintFrameBarLogical,
      Position(ScintFrameBeamLoopX/2., -(beam2stage)+sft/2., -beam2stage+sft/2.)).addPhysVolID("scintFrameBar", 3);

  // Support legs
  double dsupportl = std::sqrt(2.)*(beam2stage - sft/2.);
  Box solidSupportLegB(sft/2., dsupportl/2., sft/2.);
  Box solidSupportLegA(sft, sft, sft);
  Box solidSupportLegC(sft/2., sft/2., beam2stage/2.);

  // G4RotationMatrix rotateX(-45deg) → RotationZYX(0, 0, -pi/4)
  RotationZYX legCutRot(0., 0., -M_PI/4.);
  SubtractionSolid solidSupportLegD("scintSupportLegBoxD", solidSupportLegB, solidSupportLegA,
      Transform3D(legCutRot, Position(0.,  0.5*(dsupportl-sft*(1.-std::sqrt(2.))), sft/std::sqrt(2.))));
  SubtractionSolid solidSupportLegE("scintSupportLegBoxE", solidSupportLegD, solidSupportLegA,
      Transform3D(legCutRot, Position(0., -0.5*(dsupportl-sft*(1.-std::sqrt(2.))), sft/std::sqrt(2.))));

  Volume scintSupportLegLogical ("scintSupportLegLogical",  solidSupportLegC, frame_mat);
  Volume scintSupportLegLogicalE("scintSupportLegLogicalE", solidSupportLegE, frame_mat);

  // G4RotationMatrix rotateY(180deg) → RotationZYX(0, pi, 0)
  RotationZYX flipY(0., M_PI, 0.);
  RotationZYX legRot125(0., 0., 3.0*M_PI/4.0);

  double legposx  = (scx+sft)/2.;
  double legposy  = -(beam2stage)+sft/2.;
  double legposz  = beam2stage/2.;
  double legeposy = 0.5*(sft + (dsupportl-sft)/std::sqrt(2.));
  double legeposz = 0.5*(beam2stage - (dsupportl+sft)/std::sqrt(2.));

  std::vector<double> legposxv {legposx+ScintFrameBeamLoopX, -legposx, legposx+ScintFrameBeamLoopX, -legposx};
  std::vector<double> legposzv {legposz, legposz, -legposz, -legposz};
  std::vector<bool>   legflipv {false,   false,   true,     true};
  std::vector<double> legeposzv{legeposz, legeposz, -legeposz, -legeposz};

  for (int ii = 0; ii < legposxv.size(); ++ii) {
    Position lpos(legposxv[ii], legposy, legposzv[ii]);
    if (legflipv[ii])
      fScintArmMotherVol.placeVolume(scintSupportLegLogical,
          Transform3D(flipY, lpos)).addPhysVolID("scintSupportLeg", ii);
    else
      fScintArmMotherVol.placeVolume(scintSupportLegLogical,
          lpos).addPhysVolID("scintSupportLeg", ii);

    // Diagonal support leg E
    Position lepos(0., legeposy, legeposzv[ii]);
    RotationZYX ler = legflipv[ii] ? flipY * legRot125 : legRot125;
    fScintArmMotherVol.placeVolume(scintSupportLegLogicalE,
        Transform3D(ler, lpos + lepos)).addPhysVolID("scintSupportLegE", ii);
  }

  // -----------------------------------------------------------------------
  // Frame for BremScintArmMotherLogical
  // -----------------------------------------------------------------------
  Box solidBremScintFrameLeg(sft/2., scy/2.+sft+bsfwo, sft/2.);
  Box solidBremScintFrameBar(scx/2., sft/2., sft/2.);
  Box solidBremScintFrameSupportLeg1(sft/2.,
      std::cos(sca)*scx/2.+std::cos(sca)*sft+std::sin(sca)*sft/2., sft/2.);
  Box solidBremScintFrameSupportLeg2(sft/2., sft/2.,
      std::sin(sca)*scx/2.+std::sin(sca)*sft+std::cos(sca)*sft/2.-sft/2.);

  Volume BremScintFrameLegLogical         ("BremScintFrameLegLogical",          solidBremScintFrameLeg,          frame_mat);
  Volume BremScintFrameBarLogical         ("BremScintFrameBarLogical",          solidBremScintFrameBar,          frame_mat);
  Volume BremScintFrameSupportLeg1Logical ("BremScintFrameSupportLeg1Logical",  solidBremScintFrameSupportLeg1,  frame_mat);
  Volume BremScintFrameSupportLeg2Logical ("BremScintFrameSupportLeg2Logical",  solidBremScintFrameSupportLeg2,  frame_mat);

  // scintRot: rotateZ(-90deg) then rotateY(ScintAngle)
  // → RotationZYX(-pi/2, ScintAngle, 0)
//cl   RotationZYX scintRot(-M_PI/2., sca, 0.);
  RotationZYX scintRot(-M_PI/2., 0.0, sca);

  // brembasebarRot: rotateZ(-90deg) → RotationZYX(-pi/2, 0, 0)
  RotationZYX brembasebarRot(-M_PI/2., 0., 0.);

  // Place scintArm into BremScintArmMother
  fBremScintArmMotherVol.placeVolume(fScintArmVol,
      Transform3D(scintRot, Position(0.,0.,0.)));
//cl       Transform3D(scintRot, Position(0.,0.,0.))).addPhysVolID("BremScintArm", 1);

  // Brem frame legs
  fBremScintArmMotherVol.placeVolume(BremScintFrameLegLogical,
       Transform3D(scintRot, Position(0., -std::cos(sca)*(scx+sft)/2., -std::sin(sca)*(scx+sft)/2.))).addPhysVolID("BremScintFrameLeg", 0);
  fBremScintArmMotherVol.placeVolume(BremScintFrameLegLogical,
      Transform3D(scintRot, Position(0.,  std::cos(sca)*(scx+sft)/2.,  std::sin(sca)*(scx+sft)/2.))).addPhysVolID("BremScintFrameLeg", 1);

  // Brem frame bars
  fBremScintArmMotherVol.placeVolume(BremScintFrameBarLogical,
      Transform3D(scintRot, Position( (scy+sft)/2., 0., 0.))).addPhysVolID("BremScintFrameBar", 0);
  fBremScintArmMotherVol.placeVolume(BremScintFrameBarLogical,
      Transform3D(scintRot, Position(-(scy+sft)/2., 0., 0.))).addPhysVolID("BremScintFrameBar", 1);

  // Support legs 1 (vertical)
  double sl1z = std::sin(sca)*scx/2.+std::sin(sca)*sft+std::cos(sca)*sft/2.-sft/2.;
  fBremScintArmMotherVol.placeVolume(BremScintFrameSupportLeg1Logical,
      Position( bsfwo+(scy+3*sft)/2., 0., sl1z)).addPhysVolID("BremScintFrameSupportLeg1", 0);
  fBremScintArmMotherVol.placeVolume(BremScintFrameSupportLeg1Logical,
      Position(-bsfwo-(scy+3*sft)/2., 0., sl1z)).addPhysVolID("BremScintFrameSupportLeg1", 1);

  // Base bar
  double bbary = -(std::cos(sca)*scx/2.+std::cos(sca)*sft+std::sin(sca)*sft/2.-sft/2.);
  fBremScintArmMotherVol.placeVolume(BremScintFrameLegLogical,
      Transform3D(brembasebarRot, Position(0., bbary, sl1z))).addPhysVolID("BremScintFrameBaseBar", 0);

  // Support legs 2 (horizontal)
  fBremScintArmMotherVol.placeVolume(BremScintFrameSupportLeg2Logical,
      Position( bsfwo+(scy+3*sft)/2., bbary, -sft/2.)).addPhysVolID("BremScintFrameSupportLeg2", 0);
  fBremScintArmMotherVol.placeVolume(BremScintFrameSupportLeg2Logical,
      Position(-bsfwo-(scy+3*sft)/2., bbary, -sft/2.)).addPhysVolID("BremScintFrameSupportLeg2", 1);
}


// ---------------------------------------------------------------------------
void LxSScreen::BuildCamera(dd4hep::Detector& description)
{
  Material env_mat       = description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material cameraMat     = description.material(description.constant<std::string>("ScintCameraMaterial"));
  Material apertMat      = description.material(description.constant<std::string>("ScintCameraApertureMaterial"));
  Material shieldMat     = description.material(description.constant<std::string>("ScintCameraShieldMaterial"));

  double ScintCameraX              = description.constant<double>("ScintCameraX");
  double ScintCameraY              = description.constant<double>("ScintCameraY");
  double ScintCameraZ              = description.constant<double>("ScintCameraZ");
  double ScintCameraApertureZ      = description.constant<double>("ScintCameraApertureZ");
  double ScintCameraApertureInner  = description.constant<double>("ScintCameraApertureInner");
  double ScintCameraApertureOuter  = description.constant<double>("ScintCameraApertureOuter");
  double ScintCameraShieldThickness= description.constant<double>("ScintCameraShieldThickness");
  bool   OverlapTest               = (description.constant<int>("OverlapTest") != 0);

  // Mother box
  Box solidCameraMotherBox(ScintCameraX/2. + ScintCameraShieldThickness,
                           ScintCameraY/2. + ScintCameraShieldThickness,
                           (ScintCameraZ + ScintCameraApertureZ + ScintCameraShieldThickness)/2.);
  fCameraMotherVol = Volume("ScintCameraMotherBoxLogical", solidCameraMotherBox, env_mat);

  // Main camera box
  Box    solidCameraBox(ScintCameraX/2., ScintCameraY/2., ScintCameraZ/2.);
  Volume cameraBoxLogical("ScintCameraBoxLogical", solidCameraBox, cameraMat);
  fCameraMotherVol.placeVolume(cameraBoxLogical,
      Position(0., 0., -ScintCameraShieldThickness/2. + ScintCameraApertureZ/2.)).addPhysVolID("ScintCameraBox", 69);

  // Aperture outer (metallic rim)
  Tube   solidApertureOuter(ScintCameraApertureInner/2., ScintCameraApertureOuter/2., ScintCameraApertureZ/2., 0., 2.*M_PI);
  Volume apertureOuterLogical("ScintCameraApertureOuterLogical", solidApertureOuter, cameraMat);
  fCameraMotherVol.placeVolume(apertureOuterLogical,
      Position(0., 0., -ScintCameraShieldThickness/2. - ScintCameraZ/2.)).addPhysVolID("ScintCameraApertureOuter", 0);

  // Aperture inner (vacuum/air tube)
  Tube   solidApertureInner(0., ScintCameraApertureInner/2., ScintCameraApertureZ/2., 0., 2.*M_PI);
  Volume apertureInnerLogical("ScintCameraApertureInnerLogical", solidApertureInner, apertMat);
  fCameraMotherVol.placeVolume(apertureInnerLogical,
      Position(0., 0., -ScintCameraShieldThickness/2. - ScintCameraZ/2.)).addPhysVolID("ScintCameraApertureInner", 0);

  // Shielding — sides, top/bottom, back
  Box solidShieldSide(ScintCameraShieldThickness/2.,
                      ScintCameraY/2. + ScintCameraShieldThickness, ScintCameraZ/2.);
  Volume shieldSideLogical("ScintCameraShieldSideLogical", solidShieldSide, shieldMat);
  double sideX = (ScintCameraX + ScintCameraShieldThickness)/2.;
  double sideZ = -ScintCameraShieldThickness/2. + ScintCameraApertureZ/2.;
  fCameraMotherVol.placeVolume(shieldSideLogical, Position( sideX, 0., sideZ)).addPhysVolID("ScintCameraShieldSide", 0);
  fCameraMotherVol.placeVolume(shieldSideLogical, Position(-sideX, 0., sideZ)).addPhysVolID("ScintCameraShieldSide", 1);

  Box solidShieldTop(ScintCameraX/2., ScintCameraShieldThickness/2., ScintCameraZ/2.);
  Volume shieldTopLogical("ScintCameraShieldTopLogical", solidShieldTop, shieldMat);
  double topY = ScintCameraY/2. + ScintCameraShieldThickness/2.;
  fCameraMotherVol.placeVolume(shieldTopLogical, Position(0.,  topY, sideZ)).addPhysVolID("ScintCameraShieldTop", 0);
  fCameraMotherVol.placeVolume(shieldTopLogical, Position(0., -topY, sideZ)).addPhysVolID("ScintCameraShieldTop", 1);

  Box solidShieldBack(ScintCameraX/2. + ScintCameraShieldThickness,
                      ScintCameraY/2. + ScintCameraShieldThickness, ScintCameraShieldThickness/2.);
  Volume shieldBackLogical("ScintCameraShieldBackLogical", solidShieldBack, shieldMat);
  fCameraMotherVol.placeVolume(shieldBackLogical,
      Position(0., 0., ScintCameraZ/2. + ScintCameraApertureZ/2.)).addPhysVolID("ScintCameraShieldBack", 0);
}


// ---------------------------------------------------------------------------
dd4hep::DetElement LxSScreen::PlaceHICSScint(dd4hep::Detector& description,
                                              dd4hep::Volume&   motherVol,
                                              int               detID)
{
  double sft      = description.constant<double>("ScintFrameThickness");
  double scz1     = description.constant<double>("ScintBaseZ");
  double scz2     = description.constant<double>("ScintPhosphorZ");
  double scz3     = description.constant<double>("ScintFinishZ");
  double scx      = description.constant<double>("ScintX");
  double IPMagnetZpos     = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength= description.constant<double>("FlashMFieldLength");
  double FlashMagnetCoilZ = description.constant<double>("FlashMagnetCoilZ");
  double ComptonElBackshift = description.constant<double>("ComptonElBackshift");
  double OPPPDetZtoMagnet = description.constant<double>("OPPPDetZtoMagnet");
  double VacChambertoOPPPDetZGap = description.constant<double>("VacChambertoOPPPDetZGap");
  double ScintCameraAngle = description.constant<double>("ScintCameraAngle");
  double ScintCameraY     = description.constant<double>("ScintCameraY");
  double ScintCameraZ     = description.constant<double>("ScintCameraZ");
  double ScintCameraApertureZ = description.constant<double>("ScintCameraApertureZ");
  double ScintCameraShieldThickness = description.constant<double>("ScintCameraShieldThickness");
  double ScintCameraApertureInner   = description.constant<double>("ScintCameraApertureInner");
  double ScintCameraYpos  = description.constant<double>("ScintCameraYpos");
  double ScintCameraPlatformZ = description.constant<double>("ScintCameraPlatformZ");
  double ScintCameraSupportThickness = description.constant<double>("ScintCameraSupportThickness");
  double HICSScintilatorXPosMag1T = description.constant<double>("HICSScintilatorXPosMag1T");
  double HICSScintilatorXPosMag2T = description.constant<double>("HICSScintilatorXPosMag2T");
  double IPMagFieldY      = description.constant<double>("IPMagFieldY");
  double ScintXpos        = description.constant<double>("ScintXpos");
  int    StickScintScreentoBeamWindow = description.constant<int>("StickScintScreentoBeamWindow");
  bool   OverlapTest      = (description.constant<int>("OverlapTest") != 0);

  // Field-dependent X position
  double scintmxpos = 0.0;
  double mag1T = -10000.0 *gauss;
  double mag2T = -20000.0 *gauss;
  if (std::abs(IPMagFieldY - mag1T) < 1.0e-10)
    scintmxpos = HICSScintilatorXPosMag1T;
  else if (std::abs(IPMagFieldY - mag2T) < 1.0e-10)
    scintmxpos = HICSScintilatorXPosMag2T;
  else
    dd4hep::except("LxSScreen::PlaceHICSScint",
                   "Field %.1f is not supported!", IPMagFieldY);

  // Z position
  double comptonelzpos = IPMagnetZpos + FlashMFieldLength/2. + ComptonElBackshift;
  double dumpMagnetZ   = FlashMFieldLength;
  double gaptomag      = 0.5*(FlashMagnetCoilZ - FlashMFieldLength);
  double bpipel        = OPPPDetZtoMagnet - VacChambertoOPPPDetZGap - gaptomag;
  double cescintz;
  if (StickScintScreentoBeamWindow) {
    cescintz = IPMagnetZpos + dumpMagnetZ/2. + bpipel/2. + gaptomag
               + (OPPPDetZtoMagnet - VacChambertoOPPPDetZGap - gaptomag)/2.
               + sft/2.;
  } else {
    cescintz = comptonelzpos + (scz1+scz2+scz3)/2.;
  }

  // Place scintArmMotherLogical
  PlacedVolume pvScint = motherVol.placeVolume(fScintArmMotherVol,
      Position(-(scx/2. + scintmxpos), 0., cescintz));
  pvScint.addPhysVolID("detector", detID);
  if (OverlapTest) pvScint.ptr()->CheckOverlaps();

  DetElement det("LxHICSScint", detID);
  det.setPlacement(pvScint);

  // Clone phosphor DetElement
  DetElement phosphorDE = fPhosphorTemplate.clone("phosphor", 0);
  det.add(phosphorDE);

  // Camera rotation: rotateX(ScintCameraAngle) → RotationZYX(0, 0, ScintCameraAngle)
  RotationZYX cameraRot(0., 0., -ScintCameraAngle);

  // HICS cameras (2 placements at 1/4 and 3/4 of screen)
  double camZ = cescintz + 1.*dd4hep::m
                + (ScintCameraShieldThickness - ScintCameraApertureZ + ScintCameraZ)/2.;
  motherVol.placeVolume(fCameraMotherVol,
      Transform3D(cameraRot, Position(-(scintmxpos + 3*scx/4.), ScintCameraYpos, camZ)))
      .addPhysVolID("HICSScintCamera", 2);
  motherVol.placeVolume(fCameraMotherVol,
      Transform3D(cameraRot, Position(-(scintmxpos +   scx/4.), ScintCameraYpos, camZ)))
      .addPhysVolID("HICSScintCamera", 3);

  // Camera support leg
  Material supportMat = description.material(description.constant<std::string>("ScintCameraSupportMaterial"));
  double CeilingSurfaceYpos = description.constant<double>("CeilingSurfaceYpos");

  double camlegdy = CeilingSurfaceYpos - ScintCameraYpos + ScintCameraY/2. + ScintCameraShieldThickness;
  Tube   solidCamSupportLeg(0., ScintCameraSupportThickness/2., camlegdy/2., 0., 2.*M_PI);

  // Cut platform from support leg
  // cameraRot applied to the box cut
  Box solidCamPlatformCut(scx/2. + 5.*dd4hep::cm, ScintCameraPlatformZ/2., ScintCameraSupportThickness/2.);
  double cutZ = (camlegdy + ScintCameraSupportThickness)/2.
                - (1.-std::cos(ScintCameraAngle))*(ScintCameraY/2. + ScintCameraShieldThickness + ScintCameraSupportThickness/2.)
                - std::sin(ScintCameraAngle)*(ScintCameraPlatformZ/2. - ScintCameraZ - ScintCameraShieldThickness
                  + (ScintCameraShieldThickness + ScintCameraApertureZ + ScintCameraZ)/2.);
  SubtractionSolid solidCamSupport("ScintCameraSupport", solidCamSupportLeg, solidCamPlatformCut,
      Transform3D(cameraRot, Position(0., 0., cutZ)));
  Volume camSupportLogical("HICSScintCameraScaffoldLegLogical", solidCamSupport, supportMat);

  // cameraSupportRot: rotateZ(pi), rotateY(pi), rotateX(pi/2)
  // → RotationZYX(pi, pi, pi/2) — applies rx=pi/2, ry=pi, rz=pi
//cld   RotationZYX camSupportRot(M_PI, M_PI, -M_PI/2.);
  RotationZYX camSupportRot(0.0, 0.0, M_PI/2.0);
  double legX = -(scintmxpos + scx/2.);
  double legY = CeilingSurfaceYpos - camlegdy/2.;
  double legZ = cescintz + 1.0*dd4hep::m + (ScintCameraShieldThickness - ScintCameraApertureZ + ScintCameraZ)/2.
                + std::cos(ScintCameraAngle)*(ScintCameraPlatformZ/2. - ScintCameraZ - ScintCameraShieldThickness
                  + (ScintCameraShieldThickness + ScintCameraApertureZ + ScintCameraZ)/2.)
                + std::sin(ScintCameraAngle)*(ScintCameraY/2. + ScintCameraShieldThickness + ScintCameraSupportThickness/2.);

  motherVol.placeVolume(camSupportLogical,
      Transform3D(camSupportRot, Position(legX, legY, legZ))).addPhysVolID("HICSScintCameraSupport", 0);

  // Camera platform
  Box    solidCamPlatform(scx/2. + 5.*dd4hep::cm, ScintCameraSupportThickness/2., ScintCameraPlatformZ/2.);
  Volume camPlatformLogical("ScintCameraScaffoldPlatformLogical", solidCamPlatform, supportMat);
  double platY = ScintCameraYpos
                 - std::cos(ScintCameraAngle)*(ScintCameraY/2.+ScintCameraShieldThickness+ScintCameraSupportThickness/2.)
                 + std::sin(ScintCameraAngle)*(ScintCameraPlatformZ/2.-ScintCameraZ-ScintCameraShieldThickness
                   +(ScintCameraShieldThickness+ScintCameraApertureZ+ScintCameraZ)/2.);
  motherVol.placeVolume(camPlatformLogical,
      Transform3D(cameraRot, Position(legX, platY, legZ))).addPhysVolID("HICSScintCameraScaffoldPlatform", 0);

  return det;
}


// ---------------------------------------------------------------------------
dd4hep::DetElement LxSScreen::PlaceBremsScint(dd4hep::Detector& description,
                                               dd4hep::Volume&   motherVol,
                                               int               detID)
{
  double ScintXpos   = description.constant<double>("ScintXpos");
  double ScintZpos   = description.constant<double>("ScintZpos");
  double ScintAngle  = description.constant<double>("ScintAngle");
  double ScintX      = description.constant<double>("ScintX");
  double ScintCameraAngle = description.constant<double>("ScintCameraAngle");
  double ScintCameraY     = description.constant<double>("ScintCameraY");
  double ScintCameraZ     = description.constant<double>("ScintCameraZ");
  double ScintCameraApertureZ       = description.constant<double>("ScintCameraApertureZ");
  double ScintCameraShieldThickness = description.constant<double>("ScintCameraShieldThickness");
  double ScintCameraYpos            = description.constant<double>("ScintCameraYpos");
  double ScintCameraPlatformZ       = description.constant<double>("ScintCameraPlatformZ");
  double ScintCameraSupportThickness= description.constant<double>("ScintCameraSupportThickness");
  bool   OverlapTest = (description.constant<int>("OverlapTest") != 0);

  // Place BremScintArmMotherLogical
//   PlacedVolume pvBrem = motherVol.placeVolume(fBremScintArmMotherVol,
//       Position(0., ScintXpos, ScintZpos));
  RotationZYX BremScintRot(0.0, 0.0, ScintCameraAngle);
  PlacedVolume pvBrem = motherVol.placeVolume(fBremScintArmMotherVol, Position(0., ScintXpos, ScintZpos));
  pvBrem.addPhysVolID("detector", detID);
  if (OverlapTest) pvBrem.ptr()->CheckOverlaps();

  DetElement det("LxBremsScint", detID);
  det.setPlacement(pvBrem);

  // Clone phosphor DetElement
  DetElement phosphorDE = fPhosphorTemplate.clone("phosphor", 0);
  det.add(phosphorDE);

  // Brem cameras rotation: rotateZ(pi/2), rotateX(-ScintCameraAngle)
  // → RotationZYX(pi/2, 0, -ScintCameraAngle)
  RotationZYX BremCameraRot(M_PI/2., -ScintCameraAngle, 0.0);

  double camZ = ScintZpos + 1.5*dd4hep::m;
  motherVol.placeVolume(fCameraMotherVol,
      Transform3D(BremCameraRot, Position(-60.*dd4hep::cm,
          ScintXpos - std::cos(ScintAngle)*ScintX/4., camZ)))
      .addPhysVolID("BremScintCamera", 0);
  motherVol.placeVolume(fCameraMotherVol,
      Transform3D(BremCameraRot, Position(-60.*dd4hep::cm,
          ScintXpos + std::cos(ScintAngle)*ScintX/4., camZ)))
      .addPhysVolID("BremScintCamera", 1);

  // Brem camera support platform
  Material supportMat = description.material(description.constant<std::string>("ScintCameraSupportMaterial"));

  double camSuppX = -60.*dd4hep::cm
                    + std::cos(ScintCameraAngle)*(ScintCameraY/2.+ScintCameraShieldThickness+ScintCameraSupportThickness/2.)
                    - std::sin(ScintCameraAngle)*(ScintCameraPlatformZ/2.-ScintCameraZ-ScintCameraShieldThickness
                      +(ScintCameraShieldThickness+ScintCameraApertureZ+ScintCameraZ)/2.);
  double camSuppY = ScintXpos;
  double camSuppZ = camZ
                    + std::cos(ScintCameraAngle)*(ScintCameraPlatformZ/2.-ScintCameraZ-ScintCameraShieldThickness
                      +(ScintCameraShieldThickness+ScintCameraApertureZ+ScintCameraZ)/2.)
                    + std::sin(ScintCameraAngle)*(ScintCameraY/2.+ScintCameraShieldThickness+ScintCameraSupportThickness/2.);

  Box    solidCamPlatform(ScintX/2.+5.*dd4hep::cm, ScintCameraSupportThickness/2., ScintCameraPlatformZ/2.);
  Volume camPlatformLogical("BremScintCameraSupportPlatformLogical", solidCamPlatform, supportMat);

  // BremCameraSupportRot: rotateZ(pi/2), rotateX(-ScintCameraAngle) — same as BremCameraRot
  motherVol.placeVolume(camPlatformLogical,
      Transform3D(BremCameraRot, Position(camSuppX, camSuppY, camSuppZ)))
      .addPhysVolID("BremScintCameraScaffoldPlatform", 0);

  return det;
}


// ---------------------------------------------------------------------------
void LxSScreen::BuildSupportAssembly(dd4hep::Detector& description, dd4hep::Volume& motherVol)
{
  double scz1     = description.constant<double>("ScintBaseZ");
  double scz2     = description.constant<double>("ScintPhosphorZ");
  double scz3     = description.constant<double>("ScintFinishZ");
  double ScintXpos= description.constant<double>("ScintXpos");
  double ScintZpos= description.constant<double>("ScintZpos");
  double ScintAngle = description.constant<double>("ScintAngle");
  double ScintX   = description.constant<double>("ScintX");
  double ScintY   = description.constant<double>("ScintY");
  double ScintFrameThickness = description.constant<double>("ScintFrameThickness");
  double CerenkovAngle = description.constant<double>("CerenkovAngle");
  double CerenkovXpos  = description.constant<double>("CerenkovXpos");
  double CerenkovZpos  = description.constant<double>("CerenkovZpos");
  double CerenkovTotalBoxHeight = description.constant<double>("CerenkovTotalBoxHeight");
  double FloorSurfaceYpos = description.constant<double>("FloorSurfaceYpos");
  double ScintCameraAngle = description.constant<double>("ScintCameraAngle");
  double ScintCameraY     = description.constant<double>("ScintCameraY");
  double ScintCameraZ     = description.constant<double>("ScintCameraZ");
  double ScintCameraApertureZ       = description.constant<double>("ScintCameraApertureZ");
  double ScintCameraShieldThickness = description.constant<double>("ScintCameraShieldThickness");
  double ScintCameraPlatformZ       = description.constant<double>("ScintCameraPlatformZ");
  double ScintCameraSupportThickness= description.constant<double>("ScintCameraSupportThickness");
  double cc   = description.constant<double>("CerenkovChannels");
  double ccl  = description.constant<double>("CerenkovChannelLayers");
  double csxf = description.constant<double>("CerenkovStrawXFrequency");
  double cszf = description.constant<double>("CerenkovStrawZFrequency");
  double cslo = description.constant<double>("CerenkovStrawLayerOffset");
  double sir  = description.constant<double>("CerenkovStrawInnerRadius");
  double sglt = description.constant<double>("CerenkovStrawGraphiteLayerThickness");
  double salt = description.constant<double>("CerenkovStrawAlLayerThickness");
  double sklt = description.constant<double>("CerenkovStrawKaptonLayerThickness");
  double splt = description.constant<double>("CerenkovStrawPolyurethaneLayerThickness");
  double cspb = description.constant<double>("CerenkovStrawPlateBuffer");
  double cbt  = description.constant<double>("CerenkovBoxThickness");
  double cspt = description.constant<double>("CerenkovShieldingPlateThickness");
  bool   OverlapTest = (description.constant<int>("OverlapTest") != 0);

  double scdz       = scz1 + scz2 + scz3;
  double cherenkovdy= 2.*(((int)cc/(int)ccl-1)*csxf/2.+(ccl-1)*cslo/2.+sir+2*(sglt+salt+sklt+splt)+cspb+cbt+cspt);
  double cherenkovdz= 2.*((ccl-1)*cszf/2.+sir+2*(sglt+salt+sklt+splt)+cspb+cbt+cspt);

  double scintbypos = ScintXpos - std::cos(ScintAngle)*ScintX/2. - std::sin(ScintAngle)*scdz/2.;
  double cherbypos  = CerenkovXpos - std::cos(CerenkovAngle)*cherenkovdy/2.
                                   - std::sin(CerenkovAngle)*cherenkovdz/2.;

  double tblx   = 2.*std::max(CerenkovTotalBoxHeight, ScintY);
  double tblz   = (CerenkovZpos - ScintZpos)
                  + std::sin(CerenkovAngle)*cherenkovdy/2. + std::cos(CerenkovAngle)*cherenkovdz/2.
                  + std::sin(ScintAngle)*ScintX/2. + std::sin(ScintAngle)*ScintFrameThickness
                  + std::cos(ScintAngle)*ScintFrameThickness/2. + 10.*dd4hep::cm;

  double baseypos = std::min(scintbypos, cherbypos);
  double basezpos = CerenkovZpos + std::sin(CerenkovAngle)*cherenkovdy/2.
                    + std::cos(CerenkovAngle)*cherenkovdz/2. - tblz/2. + 5.*dd4hep::cm;

  double ypestal  = 0.1*dd4hep::m;
  double tblhight = baseypos - ypestal - FloorSurfaceYpos;

  Assembly tablesupport = LxAux::BuildTable(description, "BremDetTable", tblx, tblhight, tblz, 4);
  PlacedVolume pvTablSup = motherVol.placeVolume(tablesupport, Position(0., baseypos, basezpos));

  Volume pedestal = LxAux::BuildPedestal(description, "BremDet", tblx, ypestal, tblz);
  PlacedVolume pvPed = motherVol.placeVolume(pedestal,
      Position(0., 0.5*ypestal + FloorSurfaceYpos, basezpos));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();

  // Brem camera scaffold table
  double cambarposx = -60.*dd4hep::cm
                      + std::cos(ScintCameraAngle)*(ScintCameraY/2.+ScintCameraShieldThickness+ScintCameraSupportThickness/2.)
                      - std::sin(ScintCameraAngle)*(ScintCameraPlatformZ/2.-ScintCameraZ-ScintCameraShieldThickness
                        +(ScintCameraShieldThickness+ScintCameraApertureZ+ScintCameraZ)/2.);
  double cambarposy = ScintXpos;
  double cambarposz = (ScintZpos + 1.5*dd4hep::m)
                      + std::cos(ScintCameraAngle)*(ScintCameraPlatformZ/2.-ScintCameraZ-ScintCameraShieldThickness
                        +(ScintCameraShieldThickness+ScintCameraApertureZ+ScintCameraZ)/2.)
                      + std::sin(ScintCameraAngle)*(ScintCameraY/2.+ScintCameraShieldThickness+ScintCameraSupportThickness/2.);

  double cambardy = 2.*(ScintX/2. + 5.*dd4hep::cm);
  double cambartblx = 4.*(ScintCameraPlatformZ*std::sin(ScintCameraAngle)
                          +ScintCameraSupportThickness*std::cos(ScintCameraAngle));
  double cambartblhight = cambarposy - cambardy/2. - FloorSurfaceYpos;
  double cambartblz = 3.*(ScintCameraPlatformZ*std::cos(ScintCameraAngle)
                          +ScintCameraSupportThickness*std::sin(ScintCameraAngle));

  Assembly cambartable = LxAux::BuildTable(description, "BremScintCameraSupportTable",
                                           cambartblx, cambartblhight, cambartblz, 3);
  PlacedVolume pvCambarTbl = motherVol.placeVolume(cambartable, Position(cambarposx, cambarposy-cambardy/2., cambarposz));
}


// ---------------------------------------------------------------------------
dd4hep::DetElement LxSScreen::Place(dd4hep::Detector&  description,
                                     dd4hep::Volume&    motherVol,
                                     const std::string& detName,
                                     int                detID)
{
  auto it = fPlaceMap.find(detName);
  if (it == fPlaceMap.end())
    dd4hep::except("LxSScreen::Place", "Unknown detector name '%s'", detName.c_str());
  return it->second(description, motherVol, detID);
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------
static Ref_t create_LxSScreen(dd4hep::Detector& description,
                               xml_h             e,
                               dd4hep::SensitiveDetector sd)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  sd.setType("tracker");

  // Retrieve or create the shared LxSScreen instance
  LxSScreen* scr = nullptr;
  try {
    scr = description.extension<LxSScreen>();
    printout(INFO, "LxSScreen", "Reusing existing LxSScreen instance for '%s'.", detName.c_str());
  } catch (...) {
    printout(INFO, "LxSScreen", "Creating new LxSScreen instance for '%s'.", detName.c_str());
    scr = new LxSScreen();
    scr->Build(description, sd);
    description.addExtension<LxSScreen>(scr);
  }

  // Envelope assembly
  Volume       fLogicWorld = description.worldVolume();
  Assembly     envelope(detName + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.attr<int>(xml_tag_t("system_id")));

  dd4hep::DetElement sdet(detName, detID);
  sdet.setPlacement(envPV);

  // Dispatch to the appropriate placement method
  dd4hep::DetElement detDE = scr->Place(description, envelope, detName, detID);
  sdet.add(detDE);

  // Build support on last detector (LxBremsScint) — both detectors share one table
  bool buildSupport = x_det.hasAttr(xml_tag_t("build_support"))
                      && x_det.attr<int>(xml_tag_t("build_support"));
  if (buildSupport)
    scr->BuildSupportAssembly(description, envelope);

  return sdet;
}

DECLARE_DETELEMENT(LxSScreen, create_LxSScreen)
