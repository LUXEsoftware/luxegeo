//
/// \brief DD4hep plugin for placing LxBeamPipes sections.
///
/// Reads pipe section names from XML <pipe> child elements and
/// delegates to LxBeamPipes::Construct() for each one.
///
/// Example XML:
///   <detector id="5" name="BeamPipes" type="LxBeamPipes">
///     <pipe name="BeamPipeToDump"/>
///     <pipe name="BeamPipeToIP"/>
///     <pipe name="GammaVacuumChamber"/>
///   </detector>
//

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxBeamPipes.h"

using namespace dd4hep;


static Ref_t create_LxBeamPipes(dd4hep::Detector& description,
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

  // Pass world volume directly — beam pipe sections are placed into the world,
  // not into a sub-assembly, matching original Geant4 behaviour.
  LxBeamPipes pipes(description, fLogicWorld);

  // Iterate over <pipe name="..."/> child elements
  for (xml_coll_t c(x_det, xml_tag_t("pipe")); c; ++c) {
    xml_comp_t  x_pipe(c);
    std::string pipeName = x_pipe.attr<std::string>(xml_tag_t("name"));
    printout(INFO, "LxBeamPipes", "Constructing pipe section '%s'.", pipeName.c_str());
    pipes.Construct(pipeName);
  }

  return sdet;
}

DECLARE_DETELEMENT(LxBeamPipes, create_LxBeamPipes)
