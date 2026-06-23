//
/// \brief Implementation of the LxTargetChamber class (DD4hep version)
//

#include <cmath>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TGeoManager.h"

#include "LxAux.h"
#include "LxTargetChamber.h"

using namespace dd4hep;

static const std::string kChamberVolName = "logicTargetChamberContainer";


// ---------------------------------------------------------------------------
dd4hep::Volume LxTargetChamber::GetChamberVolume(dd4hep::Detector& description)
{
  // Check TGeoManager — return cached volume if already built
  TGeoVolume* existing = description.manager().GetVolume(kChamberVolName.c_str());
  if (existing) {
    printout(INFO, "LxTargetChamber::GetChamberVolume",
             "Reusing cached chamber volume '%s'.", kChamberVolName.c_str());
    return Volume(existing);
  }
  printout(INFO, "LxTargetChamber::GetChamberVolume",
           "Building chamber volume '%s'.", kChamberVolName.c_str());
  return BuildChamberVolume(description);
}


// ---------------------------------------------------------------------------
dd4hep::Volume LxTargetChamber::BuildChamberVolume(dd4hep::Detector& description)
{
  Material chamberMaterial = description.material(
      description.constant<std::string>("TargetChamberMaterial"));
  Material vacuumMaterial  = description.material(
      description.constant<std::string>("BeamPipeVacuumMaterial"));

  double TargetChamberX         = description.constant<double>("TargetChamberX");
  double TargetChamberY         = description.constant<double>("TargetChamberY");
  double TargetChamberZ         = description.constant<double>("TargetChamberZ");
  double TargetChamberThickness = description.constant<double>("TargetChamberThickness");
  double TargetChamberXPos      = description.constant<double>("TargetChamberXPos");
  double ChamberTargetVolX      = description.constant<double>("ChamberTargetVolX");
  double ChamberTargetVolY      = description.constant<double>("ChamberTargetVolY");
  double ChamberTargetVolZ      = description.constant<double>("ChamberTargetVolZ");
  double BPipeR                 = description.constant<double>("BPipeR");
  double BPipeThickness         = description.constant<double>("BPipeThickness");
  bool   OverlapTest            = (description.constant<int>("OverlapTest") != 0);

  double tcshiftx    = TargetChamberXPos;
  double tcwallthick = TargetChamberThickness;

  // -----------------------------------------------------------------------
  // Container: chamber outer box minus target volume cut-out
  // -----------------------------------------------------------------------
  Box solidTargetChamberContainer1(TargetChamberX/2.0, TargetChamberY/2.0, TargetChamberZ/2.0);
  Box solidChamberTargetVolume(ChamberTargetVolX/2.0, ChamberTargetVolY/2.0, ChamberTargetVolZ/2.0);

  SubtractionSolid solidTargetChamberContainer("solidTargetChamberContainer",
      solidTargetChamberContainer1, solidChamberTargetVolume,
      Position(-tcshiftx, 0.0, 0.0));

  Volume logicTargetChamberContainer(kChamberVolName, solidTargetChamberContainer, vacuumMaterial);

  // -----------------------------------------------------------------------
  // Chamber wall: hollow box with beam pipe holes on ±Z faces
  // -----------------------------------------------------------------------
  Box solidTargetChamberBoxOut(TargetChamberX/2.0, TargetChamberY/2.0, TargetChamberZ/2.0);
  Box solidTargetChamberBoxIn (TargetChamberX/2.0 - tcwallthick,
                               TargetChamberY/2.0 - tcwallthick,
                               TargetChamberZ/2.0 - tcwallthick);
  Tube solidBeamPipeCut(0.0, BPipeR - BPipeThickness, tcwallthick, 0.0, 2.0*M_PI);

  SubtractionSolid solidTargetChamberBox1("solidTargetChamberBox1",
      solidTargetChamberBoxOut, solidTargetChamberBoxIn);

  // Beam pipe hole on +Z face
  SubtractionSolid solidTargetChamberBox2("solidTargetChamberBox2",
      solidTargetChamberBox1, solidBeamPipeCut,
      Position(-tcshiftx, 0.0, TargetChamberZ/2.0));

  // Beam pipe hole on -Z face
  SubtractionSolid solidTargetChamberBox("solidTargetChamberBox",
      solidTargetChamberBox2, solidBeamPipeCut,
      Position(-tcshiftx, 0.0, -TargetChamberZ/2.0));

  Volume logicTargetChamberBox("logicTargetChamberBox", solidTargetChamberBox, chamberMaterial);

  PlacedVolume pvBox = logicTargetChamberContainer.placeVolume(logicTargetChamberBox,
      Position(0.0, 0.0, 0.0));
  pvBox.addPhysVolID("TargetChamberBox", 0);
  if (OverlapTest) pvBox.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Target support assembly (frame + holder + motor)
  // -----------------------------------------------------------------------
  Assembly targetSupportAssembly = ConstructTargetSupportAssembly(description);
//   LxAux::AddAssemblyVolumes(logicTargetChamberContainer, targetSupportAssembly,
//                             Position(-tcshiftx, 0.0, 0.0));
  PlacedVolume pvSupportAssy = logicTargetChamberContainer.placeVolume(targetSupportAssembly,Position(-tcshiftx, 0.0, 0.0));

  return logicTargetChamberContainer;
}


// ---------------------------------------------------------------------------
dd4hep::Volume LxTargetChamber::GetTargetVolume(dd4hep::Detector&  description,
                                                 const std::string& volname)
{
  Material vacuumMaterial = description.material(
      description.constant<std::string>("BeamPipeVacuumMaterial"));

  double ChamberTargetVolX = description.constant<double>("ChamberTargetVolX");
  double ChamberTargetVolY = description.constant<double>("ChamberTargetVolY");
  double ChamberTargetVolZ = description.constant<double>("ChamberTargetVolZ");

  std::string sname = "solid" + volname + "Container";
  std::string lname = "logic" + volname + "Container";

  Box    solidChamberTargetVolume(sname, ChamberTargetVolX/2.0, ChamberTargetVolY/2.0, ChamberTargetVolZ/2.0);
  Volume logicChamberTargetVolume(lname, solidChamberTargetVolume, vacuumMaterial);
  return logicChamberTargetVolume;
}


// ---------------------------------------------------------------------------
void LxTargetChamber::ConstructSupport(dd4hep::Detector&  description,
                                       dd4hep::Volume&    motherVol,
                                       const std::string& tcname,
                                       double             zpos)
{
  double TargetChamberX   = description.constant<double>("TargetChamberX");
  double TargetChamberY   = description.constant<double>("TargetChamberY");
  double TargetChamberZ   = description.constant<double>("TargetChamberZ");
  double TargetChamberXPos= description.constant<double>("TargetChamberXPos");
  double FloorSurfaceYpos = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest      = (description.constant<int>("OverlapTest") != 0);

  double ypestal  = 0.5*dd4hep::m;
  double ylevel   = TargetChamberY/2.0;
  double tblhight = -ylevel - ypestal - FloorSurfaceYpos;
  double tblx     = 1.5*TargetChamberX;
  double tblz     = 2.0*TargetChamberZ;

  Assembly tablesupport = LxAux::BuildTable(description, tcname + "Support",
                                            tblx, tblhight, tblz, 1);
//   LxAux::AddAssemblyVolumes(motherVol, tablesupport,
//                             Position(TargetChamberXPos, -ylevel, zpos));
  PlacedVolume pvTableAssy = motherVol.placeVolume(tablesupport, Position(TargetChamberXPos, -ylevel, zpos));

  Volume pedestal = LxAux::BuildPedestal(description, tcname + "Support", tblx, ypestal, tblz);
  PlacedVolume pvPed = motherVol.placeVolume(pedestal,
      Position(TargetChamberXPos, 0.5*ypestal + FloorSurfaceYpos, zpos));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTargetChamber::ConstructTargetSupportAssembly(dd4hep::Detector& description)
{
  Material targetMotorMaterial = description.material(
      description.constant<std::string>("TargetChamberMotorMaterial"));
  Material targetFrameMaterial = description.material(
      description.constant<std::string>("TargetFrameMaterial"));

  double TargetChamberMotorX    = description.constant<double>("TargetChamberMotorX");
  double TargetChamberMotorY    = description.constant<double>("TargetChamberMotorY");
  double TargetChamberMotorZ    = description.constant<double>("TargetChamberMotorZ");
  double TargetChamberY         = description.constant<double>("TargetChamberY");
  double TargetChamberThickness = description.constant<double>("TargetChamberThickness");
  double ChamberTargetVolX      = description.constant<double>("ChamberTargetVolX");
  double ChamberTargetVolY      = description.constant<double>("ChamberTargetVolY");
  double TargetFrameD           = description.constant<double>("TargetFrameD");
  double TargetFrameZ           = description.constant<double>("TargetFrameZ");
  double TargetFrameHolderR     = description.constant<double>("TargetFrameHolderR");
  double TargetFrameHolderGroveY= description.constant<double>("TargetFrameHolderGroveY");

  // Motor
  Box    solidTCMotorH(TargetChamberMotorX/2.0, TargetChamberMotorY/2.0, TargetChamberMotorZ/2.0);
  Volume logicTCMotorH("logicTCMotorH", solidTCMotorH, targetMotorMaterial);

  // Target frame (L-shaped: box minus cut-out)
  double framex = ChamberTargetVolX + TargetFrameD;
  double framey = ChamberTargetVolY + 2.0*TargetFrameD;

  Box solidTargetFrame1(framex/2.0, framey/2.0, TargetFrameZ/2.0);
  Box solidTargetFrameCut(framex/2.0, ChamberTargetVolY/2.0, TargetFrameZ);

  SubtractionSolid solidTargetFrame("solidTargetFrame", solidTargetFrame1, solidTargetFrameCut,
      Position(-TargetFrameD, 0.0, 0.0));
  Volume logicTargetFrame("logicTargetFrame", solidTargetFrame, targetFrameMaterial);

  // Frame holder (tube with groove cut)
  double tfholdy = (TargetChamberY - framey)/2.0 - TargetChamberThickness
                   - TargetChamberMotorY + TargetFrameHolderGroveY;

  Tube solidTFrameHolder1(0.0, TargetFrameHolderR, tfholdy/2.0, 0.0, 2.0*M_PI);
  Box  solidTFrameHolderCut(2.0*TargetFrameHolderR, TargetFrameZ/2.0, TargetFrameHolderGroveY);

  SubtractionSolid solidTFrameHolder("solidTFrameHolder", solidTFrameHolder1, solidTFrameHolderCut,
      Position(0.0, 0.0, tfholdy/2.0));
  Volume logicTFrameHolder("logicTFrameHolder", solidTFrameHolder, targetFrameMaterial);

  // Assembly
  Assembly targetFrameAssembly("TargetFrameAssembly");

  targetFrameAssembly.placeVolume(logicTargetFrame,
      Position(TargetFrameD/2.0, 0.0, 0.0));

  // G4RotationMatrix(G4ThreeVector(-1,0,0), pi/2): axis-angle around -X by pi/2
  // → RotationZYX(0, 0, -pi/2)
  Position trahold(0.0, TargetFrameHolderGroveY - (framey + tfholdy)/2.0, 0.0);
  targetFrameAssembly.placeVolume(logicTFrameHolder,
      Transform3D(RotationZYX(0.0, 0.0, -M_PI/2.0), trahold));

  double motorx = 0.5*(TargetChamberMotorX - framex + TargetFrameD);
  Position tramotor(motorx, trahold.Y() - (TargetChamberMotorY + tfholdy)/2.0, 0.0);
  targetFrameAssembly.placeVolume(logicTCMotorH, tramotor);

  return targetFrameAssembly;
}
