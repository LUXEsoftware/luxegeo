//
/// \brief Implementation of the LxXTD8 class (DD4hep version)
//

#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxXTD8.h"

using namespace dd4hep;


void LxXTD8::Construct(dd4hep::Detector&   description,
                       dd4hep::DetElement& sdet,
                       xml_h&              e)
{
//   description.manager().GetElementTable()->Print();

  // World (mother) volume — equivalent of fLogicWorld in the Geant4 version
  Volume fLogicWorld = description.worldVolume();

  // DD4hep requires every DetElement to have a valid placement.
  // An Assembly envelope with an identity placement in the world satisfies
  // this requirement without introducing any extra geometry.
  xml_comp_t x_det(e);
  Assembly   envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);

  // All sub-volumes are placed inside the assembly
  Volume& motherVol = envelope;

  CreateMaterial();

  // Read geometry constants from the <define> section of the XML.
  // Variable names match lxs->MemberName in the original code.
  double TunnelRin        = description.constant<double>("TunnelRin");
  double TunnelThickness  = description.constant<double>("TunnelThickness");
  double TunnelLength     = description.constant<double>("TunnelLength");
  double TunnelWallXpos   = description.constant<double>("TunnelWallXpos");
  double FloorX           = description.constant<double>("FloorX");
  double FloorSurfaceYpos = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest      = (description.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Primary solids
  // Tube  parameters: { rmin, rmax, hz, startPhi, endPhi }
  // Box   parameters: { dx,   dy,   dz }   (full lengths; halved in the loop)
  // -----------------------------------------------------------------------
  std::map<std::string, std::tuple<std::string, std::vector<double>>> ugb;

  ugb["TunnelA"] = std::make_tuple("Tube", std::vector<double>(
      { TunnelRin, TunnelRin + TunnelThickness,
        TunnelLength, 0.0, 0.5*M_PI }));
  ugb["TunnelB"] = std::make_tuple("Tube", std::vector<double>(
      { TunnelRin, TunnelRin + TunnelThickness,
        TunnelLength, 0.5*M_PI, 1.0*M_PI }));
  ugb["TunnelC"] = std::make_tuple("Tube", std::vector<double>(
      { TunnelRin, TunnelRin + TunnelThickness,
        TunnelLength, 1.0*M_PI, 1.5*M_PI }));
  ugb["TunnelD"] = std::make_tuple("Tube", std::vector<double>(
      { TunnelRin, TunnelRin + TunnelThickness,
        TunnelLength, 1.5*M_PI, 2.0*M_PI }));
  ugb["Floor1"]  = std::make_tuple("Tube", std::vector<double>(
      { 0.0, TunnelRin, TunnelLength, 0.0, 2.0*M_PI }));
  ugb["FloorCut"]= std::make_tuple("Box",  std::vector<double>(
      { 3.0*TunnelRin, 2.0*TunnelRin, 2.0*TunnelLength }));

  // -----------------------------------------------------------------------
  // Composite solids: { operation, solid1, solid2 }
  // operation 0 = subtraction
  // Ordered vector of pairs: entries are processed in insertion order, which
  // is essential when a later cut operand is itself the result of an earlier
  // composite (e.g. successive subtractions from the same solid).
  // -----------------------------------------------------------------------
  std::vector<std::pair<std::string, std::tuple<int, std::string, std::string>>> compositsolid;
  compositsolid.push_back({"Floor", std::make_tuple(0, "Floor1", "FloorCut")});

  // -----------------------------------------------------------------------
  // Logic volume materials
  // Ordered vector of pairs to guarantee deterministic placement sequence.
  // -----------------------------------------------------------------------
  std::vector<std::pair<std::string, std::string>> material;
  material.push_back({"TunnelA", "ShieldingConcrete"});
  material.push_back({"TunnelB", "ShieldingConcrete"});
  material.push_back({"TunnelC", "ShieldingConcrete"});
  material.push_back({"TunnelD", "ShieldingConcrete"});
  material.push_back({"Floor",   "ShieldingConcrete"});

  // -----------------------------------------------------------------------
  // Transformations: { rx, ry, rz, tx, ty, tz }
  // Used both for placements and for composite-solid cut positions.
  // -----------------------------------------------------------------------
  double yrot      = 0.0;
  double floorcuty = TunnelRin - std::sqrt(TunnelRin*TunnelRin
                                           - FloorX*FloorX / 4.0);
  double tdx = TunnelWallXpos - TunnelRin;
  double tdy = TunnelRin - (floorcuty - FloorSurfaceYpos);
  double tdz = 0.0;

  std::map<std::string, std::tuple<double,double,double,double,double,double>> trns;
  trns["TunnelA"]             = std::make_tuple(0.0, yrot, 0.0, tdx, tdy, tdz);
  trns["TunnelB"]             = std::make_tuple(0.0, yrot, 0.0, tdx, tdy, tdz);
  trns["TunnelC"]             = std::make_tuple(0.0, yrot, 0.0, tdx, tdy, tdz);
  trns["TunnelD"]             = std::make_tuple(0.0, yrot, 0.0, tdx, tdy, tdz);
  trns["Floor"]               = std::make_tuple(0.0, yrot, 0.0, tdx, tdy, tdz);
  trns["Cut_Floor1_FloorCut"] = std::make_tuple(0.0, 0.0,  0.0, 0.0, floorcuty, 0.0);

  // -----------------------------------------------------------------------
  // Loop: build primary solids
  // -----------------------------------------------------------------------
  std::map<std::string, Solid> volsolid;

  for (const auto& ssolid : ugb) {
    const std::string&         vname     = ssolid.first;
    const std::string&         solidtype = std::get<0>(ssolid.second);
    const std::vector<double>& parv      = std::get<1>(ssolid.second);

    if (solidtype == "Box") {
      volsolid[vname] = Box(parv.at(0)/2.0, parv.at(1)/2.0, parv.at(2)/2.0);
    } else if (solidtype == "Tube") {
      // parv: { rmin, rmax, length, startPhi, endPhi }
      // DD4hep Tube: (rmin, rmax, hz, startPhi, endPhi)
      volsolid[vname] = Tube(parv.at(0), parv.at(1), parv.at(2)/2.0,
                             parv.at(3), parv.at(4));
    } else {
      printout(WARNING, "LxXTD8::Construct",
               "Solid type %s is not supported! Ignoring.", solidtype.c_str());
    }
  }

  // -----------------------------------------------------------------------
  // Loop: build composite solids
  // -----------------------------------------------------------------------
  for (const auto& csolid : compositsolid) {
    const std::string& cname  = csolid.first;
    const int&         optype = std::get<0>(csolid.second);
    const std::string& solv1  = std::get<1>(csolid.second);
    const std::string& solv2  = std::get<2>(csolid.second);

    if (optype == 0) {  // subtraction
      std::string trnsname = std::string("Cut_") + solv1 + std::string("_") + solv2;
      double rx, ry, rz, tx, ty, tz;
      std::tie(rx, ry, rz, tx, ty, tz) = trns[trnsname];

      Transform3D cutTrf(RotationZYX(rz, ry, rx), Position(tx, ty, tz));
      volsolid[cname] = SubtractionSolid(volsolid[solv1], volsolid[solv2], cutTrf);
    } else {
      printout(WARNING, "LxXTD8::Construct",
               "Composite solid operation %d is not supported! Ignoring.", optype);
    }
  }

  // -----------------------------------------------------------------------
  // Loop: create logical volumes and place them into the assembly
  // -----------------------------------------------------------------------
  int volID = 1;
  for (const auto& volm : material) {
    const std::string& vname = volm.first;

    Material vmt = description.material(volm.second);
    Volume   lv(std::string("logic") + vname, volsolid[vname], vmt);

    double rx, ry, rz, tx, ty, tz;
    std::tie(rx, ry, rz, tx, ty, tz) = trns[vname];

    Transform3D placeTrf(RotationZYX(rz, ry, rx), Position(tx, ty, tz));
    PlacedVolume pv = motherVol.placeVolume(lv, placeTrf);
    pv.addPhysVolID(vname, volID++);

    if (OverlapTest) pv.ptr()->CheckOverlaps();
  }

}


void LxXTD8::CreateMaterial()
{
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxXTD8(dd4hep::Detector& description,
                            xml_h             e,
                            dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxXTD8 lxXTD8;
  lxXTD8.Construct(description, sdet, e);

  return sdet;
}

DECLARE_DETELEMENT(LxXTD8, create_LxXTD8)
