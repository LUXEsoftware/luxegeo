//
/// \brief Implementation of the LxBremsDump class (DD4hep version)
//

#include <cmath>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxBremsDump.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
void LxBremsDump::Construct(dd4hep::Detector&   description,
                             dd4hep::DetElement& sdet,
                             xml_h&              e)
{
  Volume fLogicWorld = description.worldVolume();

  xml_comp_t   x_det(e);
  Assembly     envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);

  ConstructBeamDump(description, envelope);
  ConstructShielding(description, envelope);
}


// ---------------------------------------------------------------------------
void LxBremsDump::ConstructBeamDump(dd4hep::Detector& description,
                                    dd4hep::Volume&   motherVol)
{
  Material beamDumpMaterial   = description.material(description.constant<std::string>("BeamDumpMaterial"));
  Material beamPipeMaterial   = description.material(description.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial     = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));
  Material dumpInsertMaterial = description.material(description.constant<std::string>("BeamDumpInsertMaterial"));

  double BeamDumpR          = description.constant<double>("BeamDumpR");
  double BeamDumpZ          = description.constant<double>("BeamDumpZ");
  double BeamDumpAngle      = description.constant<double>("BeamDumpAngle");
  double BeamDumpXpos       = description.constant<double>("BeamDumpXpos");
  double BeamDumpZpos       = description.constant<double>("BeamDumpZpos");
  double BeamDumpFrontXpos  = description.constant<double>("BeamDumpFrontXpos");
  double BeamDumpInsertR    = description.constant<double>("BeamDumpInsertR");
  double BeamDumpInsertZ    = description.constant<double>("BeamDumpInsertZ");
  double BPipeR             = description.constant<double>("BPipeR");
  double BPipeThickness     = description.constant<double>("BPipeThickness");
  bool   OverlapTest        = (description.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Main dump body
  // -----------------------------------------------------------------------
  double dx_dump_pipe = BeamDumpFrontXpos - BeamDumpR * std::cos(BeamDumpAngle) - BPipeR;

  Solid solidBeamDump;
  if (dx_dump_pipe < 0.0) {
    // Dump overlaps beam pipe region — subtract beam pipe hole
    Tube solidBeamDump0(0.0, BeamDumpR, BeamDumpZ/2.0, 0.0, 2.0*M_PI);
    double lbpipecut = BeamDumpZ / std::cos(BeamDumpAngle)
                       + 2.0*BPipeR * std::tan(BeamDumpAngle);
    Tube solidBeamPipeHole(0.0, BPipeR, 0.5*lbpipecut, 0.0, 2.0*M_PI);

    // Commented: alternative rotation around Y
    // Transform3D(RotationZYX(0.0, BeamDumpAngle, 0.0),
    //             Position(BeamDumpXpos/std::cos(BeamDumpAngle), 0.0, 0.0))

    // G4RotationMatrix(G4ThreeVector(-1,0,0), angle): axis-angle around -X by angle
    // → rotation around X by -angle → RotationZYX(0, 0, -angle)
    solidBeamDump = SubtractionSolid("solidBeamDump", solidBeamDump0, solidBeamPipeHole,
        Transform3D(RotationZYX(0.0, 0.0, -BeamDumpAngle),
                    Position(0.0, BeamDumpXpos / std::cos(BeamDumpAngle), 0.0)));
  } else {
    solidBeamDump = Tube(0.0, BeamDumpR, BeamDumpZ/2.0, 0.0, 2.0*M_PI);
  }

  Volume logicBeamDump("logicBeamDump", solidBeamDump, beamDumpMaterial);

  // -----------------------------------------------------------------------
  // Dump insert (front tungsten/copper plug)
  // -----------------------------------------------------------------------
  Tube   solidBeamDumpInsert(0.0, BeamDumpInsertR, BeamDumpInsertZ/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamDumpInsert("logicBeamDumpInsert", solidBeamDumpInsert, dumpInsertMaterial);

  PlacedVolume pvInsert = logicBeamDump.placeVolume(logicBeamDumpInsert,
      Position(0.0, 0.0, (BeamDumpInsertZ - BeamDumpZ)/2.0));
  pvInsert.addPhysVolID("BeamDumpInsert", 0);
  if (OverlapTest) pvInsert.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Place dump into world
  // Commented: alternative rotation around Y
  // motherVol.placeVolume(logicBeamDump,
  //     Transform3D(RotationZYX(0.0, BeamDumpAngle, 0.0),
  //                 Position(-BeamDumpXpos, 0.0, BeamDumpZpos)));

  // G4RotationMatrix(G4ThreeVector(-1,0,0), angle) → RotationZYX(0, 0, -angle)
  PlacedVolume pvDump = motherVol.placeVolume(logicBeamDump,
      Transform3D(RotationZYX(0.0, 0.0, BeamDumpAngle),
                  Position(0.0, -BeamDumpXpos, BeamDumpZpos)));
  pvDump.addPhysVolID("BeamDumpAssembly", 0);
  if (OverlapTest) pvDump.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Beam pipe section alongside dump
  // -----------------------------------------------------------------------
  double lbpipe = (2.0*BeamDumpR * std::tan(BeamDumpAngle) + BeamDumpZ)
                  * std::cos(BeamDumpAngle);

  Tube   solidBeamPipeNextToDump(BPipeR - BPipeThickness, BPipeR, 0.5*lbpipe, 0.0, 2.0*M_PI);
  Volume logicBeamPipeNextToDump("logicBeamPipeNextToDump", solidBeamPipeNextToDump, beamPipeMaterial);
  PlacedVolume pvPipe = motherVol.placeVolume(logicBeamPipeNextToDump,
      Position(0.0, 0.0, BeamDumpZpos));
  pvPipe.addPhysVolID("BeamPipeNextToDump", 0);
  if (OverlapTest) pvPipe.ptr()->CheckOverlaps();

  Tube   solidBeamPipeNextToDumpVac(0.0, BPipeR - BPipeThickness, 0.5*lbpipe, 0.0, 2.0*M_PI);
  Volume logicBeamPipeNextToDumpVac("logicBeamPipeNextToDumpVac", solidBeamPipeNextToDumpVac, vacuumMaterial);
  PlacedVolume pvPipeVac = motherVol.placeVolume(logicBeamPipeNextToDumpVac,
      Position(0.0, 0.0, BeamDumpZpos));
  pvPipeVac.addPhysVolID("BeamPipeNextToDumpVac", 0);
  if (OverlapTest) pvPipeVac.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
void LxBremsDump::ConstructShielding(dd4hep::Detector& description,
                                     dd4hep::Volume&   motherVol)
{
  Material shieldingMaterial    = description.material(description.constant<std::string>("ShieldingMaterial"));
  Material beamPipeMaterial     = description.material(description.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial       = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));
  Material shieldingAbsMaterial = description.material(description.constant<std::string>("ShieldingAbsorberMaterial"));
  Material environmentMaterial  = description.material(description.constant<std::string>("EnvironmentMaterial"));

  double BeamDumpR              = description.constant<double>("BeamDumpR");
  double BeamDumpZ              = description.constant<double>("BeamDumpZ");
  double BeamDumpAngle          = description.constant<double>("BeamDumpAngle");
  double BeamDumpXpos           = description.constant<double>("BeamDumpXpos");
  double BeamDumpZpos           = description.constant<double>("BeamDumpZpos");
  double BeamDumpFrontXpos      = description.constant<double>("BeamDumpFrontXpos");
  double BPipeR                 = description.constant<double>("BPipeR");
  double BPipeThickness         = description.constant<double>("BPipeThickness");
  double ShieldingX             = description.constant<double>("ShieldingX");
  double ShieldingY             = description.constant<double>("ShieldingY");
  double ShieldingZ             = description.constant<double>("ShieldingZ");
  double ShieldingDeepZ         = description.constant<double>("ShieldingDeepZ");
  double ShieldingDeepMargine   = description.constant<double>("ShieldingDeepMargine");
  double ShieldingDXBeamWall    = description.constant<double>("ShieldingDXBeamWall");
  double FloorSurfaceYpos       = description.constant<double>("FloorSurfaceYpos");
  double ShieldingAbsorberX     = description.constant<double>("ShieldingAbsorberX");
  double ShieldingAbsorberZ     = description.constant<double>("ShieldingAbsorberZ");
  double ShieldingAbsorberSlitX = description.constant<double>("ShieldingAbsorberSlitX");
  double ShieldingAbsorberSlitZ = description.constant<double>("ShieldingAbsorberSlitZ");
  double ShieldingAbsorberTopZ  = description.constant<double>("ShieldingAbsorberTopZ");
  bool   OverlapTest            = (description.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Geometry calculations
  // -----------------------------------------------------------------------
  double ydeep     = 2.0 * (BeamDumpR + ShieldingDeepMargine);
  double bdumpxproj= 2.0*BeamDumpR * std::cos(BeamDumpAngle)
                     + BeamDumpZ * std::sin(BeamDumpAngle);

  double xdeep, xdeeppos;
  if (BeamDumpFrontXpos + BPipeR < BeamDumpR * std::cos(BeamDumpAngle)) {
    xdeep    = bdumpxproj + 2.0*ShieldingDeepMargine;
    xdeeppos = BeamDumpXpos;
  } else {
    xdeep    = BeamDumpXpos + BPipeR + 0.5*bdumpxproj + 2.0*ShieldingDeepMargine;
    xdeeppos = 0.5*xdeep - BPipeR - ShieldingDeepMargine;
  }

  double ddd = 0.5*(2.0*BeamDumpR * std::sin(BeamDumpAngle)
                    + BeamDumpZ * std::cos(BeamDumpAngle))
               + ShieldingDeepMargine;
  double zdeeppos  = 0.5 * (ShieldingDeepZ - ShieldingZ);
  double zshieldpos= BeamDumpZpos - ShieldingDeepZ + 0.5*ShieldingZ + ddd;

  double walldx  = ShieldingX/2.0 - ShieldingDXBeamWall;
  double floordy = -(ShieldingY/2.0 + FloorSurfaceYpos);

  // -----------------------------------------------------------------------
  // Main shielding block (with beam pipe hole and deep dump cavity)
  // -----------------------------------------------------------------------
  Box  solidShielding0(ShieldingX/2.0, ShieldingY/2.0, ShieldingZ/2.0);
  Tube solidShieldingH(0.0, BPipeR, 0.55*ShieldingZ, 0.0, 2.0*M_PI);

  // Commented: original xdeep/ydeep ordering (geometry was in XZ plane)
  // Box solidShieldingDeep(xdeep/2.0, ydeep/2.0, ShieldingDeepZ/2.0 + ShieldingDeepMargine);
  // Transform3D transformd(RotationZYX(), Position(-xdeeppos + walldx, 0.0, zdeeppos));

  // Active version: dump is in YZ plane so X and Y are swapped in the deep cavity
  Box solidShieldingDeep(ydeep/2.0, xdeep/2.0, ShieldingDeepZ/2.0 + ShieldingDeepMargine);

  SubtractionSolid solidShielding1("solidShielding1", solidShielding0, solidShieldingH,
      Position(walldx, floordy, 0.0));

  SubtractionSolid solidShielding("solidShielding", solidShielding1, solidShieldingDeep,
      Position(walldx, -xdeeppos + floordy, zdeeppos));

  Volume logicShielding("logicShielding", solidShielding, shieldingMaterial);

  PlacedVolume pvShield = motherVol.placeVolume(logicShielding,
      Position(-walldx, -floordy, zshieldpos));
  pvShield.addPhysVolID("Shielding", 0);
  if (OverlapTest) pvShield.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Front lead absorber bar (inside shielding logical volume)
  // -----------------------------------------------------------------------
  double shabsorby = (floordy - xdeeppos) + ShieldingY/2.0 - xdeep/2.0;

  Box    solidShieldAbsorber(ShieldingAbsorberX/2.0, shabsorby/2.0, ShieldingAbsorberZ/2.0);
  Volume logicShieldAbsorber("logicShieldAbsorber", solidShieldAbsorber, shieldingAbsMaterial);

  PlacedVolume pvAbs = logicShielding.placeVolume(logicShieldAbsorber,
      Position(walldx, 0.5*(shabsorby - ShieldingY),
               0.5*(ShieldingAbsorberZ - ShieldingZ)));
  pvAbs.addPhysVolID("ShieldingAbsorber", 0);
  if (OverlapTest) pvAbs.ptr()->CheckOverlaps();

  // Air slit in front absorber bar (inside absorber logical volume)
  Box    solidShieldAbsSlit(ShieldingAbsorberSlitX/2.0, shabsorby/2.0, ShieldingAbsorberSlitZ/2.0);
  Volume logicShieldAbsSlit("logicShieldAbsSlit", solidShieldAbsSlit, environmentMaterial);

  PlacedVolume pvSlit = logicShieldAbsorber.placeVolume(logicShieldAbsSlit,
      Position(0.0, 0.0, 0.5*(ShieldingAbsorberSlitZ - ShieldingAbsorberZ)));
  pvSlit.addPhysVolID("ShieldingAbsorberSlit", 0);
  if (OverlapTest) pvSlit.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Top lead absorber bar (covers gap between shielding and dump)
  // placed directly into world
  // -----------------------------------------------------------------------
  double shabsorbtopy = ShieldingDeepMargine + BeamDumpR/2.0;

  Box solidShieldAbsorberTop1(ShieldingAbsorberX/2.0, shabsorbtopy/2.0, ShieldingAbsorberTopZ/2.0);
  Box solidShieldAbsorberTopCut(ShieldingAbsorberSlitX/2.0, shabsorbtopy, ShieldingAbsorberSlitZ/2.0);

  SubtractionSolid solidShieldAbsorberTop("solidShieldAbsorberTop",
      solidShieldAbsorberTop1, solidShieldAbsorberTopCut,
      Position(0.0, 0.0, 0.5*(ShieldingAbsorberSlitZ - ShieldingAbsorberTopZ)));

  Volume logicShieldAbsorberTop("logicShieldAbsorberTop", solidShieldAbsorberTop, shieldingAbsMaterial);

  double shabstopposy = shabsorby + shabsorbtopy/2.0 + FloorSurfaceYpos;
  double shabstopposz = zshieldpos + 0.5*(ShieldingAbsorberTopZ - ShieldingZ);

  PlacedVolume pvAbsTop = motherVol.placeVolume(logicShieldAbsorberTop,
      Position(0.0, shabstopposy, shabstopposz));
  pvAbsTop.addPhysVolID("ShieldAbsorberTop", 0);
  if (OverlapTest) pvAbsTop.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Beam pipe through shielding
  // -----------------------------------------------------------------------
  double lshildpipe   = ShieldingZ - ShieldingDeepZ + ShieldingDeepMargine;
  double shildpipezpos= zshieldpos + 0.5*(ShieldingZ - lshildpipe);

  Tube   solidShieldingPipe(BPipeR - BPipeThickness, BPipeR, lshildpipe/2.0, 0.0, 2.0*M_PI);
  Volume logicShieldingPipe("logicShieldingPipe", solidShieldingPipe, beamPipeMaterial);
  PlacedVolume pvSPipe = motherVol.placeVolume(logicShieldingPipe,
      Position(0.0, 0.0, shildpipezpos));
  pvSPipe.addPhysVolID("ShieldingPipe", 0);
  if (OverlapTest) pvSPipe.ptr()->CheckOverlaps();

  Tube   solidShieldingPipeVac(0.0, BPipeR - BPipeThickness, lshildpipe/2.0, 0.0, 2.0*M_PI);
  Volume logicShieldingPipeVac("logicShieldingPipeVac", solidShieldingPipeVac, vacuumMaterial);
  PlacedVolume pvSPipeVac = motherVol.placeVolume(logicShieldingPipeVac,
      Position(0.0, 0.0, shildpipezpos));
  pvSPipeVac.addPhysVolID("ShieldingPipeVac", 0);
  if (OverlapTest) pvSPipeVac.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxBremsDump(dd4hep::Detector& description,
                                 xml_h             e,
                                 dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxBremsDump bd;
  bd.Construct(description, sdet, e);

  return sdet;
}

DECLARE_DETELEMENT(LxBremsDump, create_LxBremsDump)
