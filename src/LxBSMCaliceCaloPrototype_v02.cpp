//====================================================================
//  DD4hep Geometry driver for Sampling Calo BOX prototype
//  adapted for LUXE BSM detector
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:  $
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDSegmentation/TiledLayerGridXY.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

#include "LxAux.h"

using namespace std;

using dd4hep::_toString;
using dd4hep::Box;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::Layer;
using dd4hep::Layering;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Segmentation;
using dd4hep::SensitiveDetector;
using dd4hep::Volume;

Volume CaliceConstructCasing(dd4hep::Detector& description, const double Calo_dim_y, const double Calo_dim_z)
{
  Material ecalCasingMaterial = description.material(
      description.constant<std::string>("Calice_Ecal_casing_material"));
  Material envMaterial = description.material(description.constant<std::string>("EnvironmentMaterial"));
  bool   OverlapTest               = (description.constant<int>("OverlapTest") != 0);

  double Ecal_casing_thickness = description.constant<double>("Calice_Ecal_casing_thickness");
  double Ecal_casing_box_x   = description.constant<double>("Calice_Ecal_casing_box_x");
  double Ecal_casing_box_y   = Calo_dim_y + 2.0 * Ecal_casing_thickness;
  double Ecal_casing_box_z   = Calo_dim_z + Ecal_casing_thickness;

  // Hollow box: outer box minus inner box offset by thickness along Z
  // so the back face is open (no wall at +Z end)
  Box solidCaliceECalCasing1(Ecal_casing_box_x/2., Ecal_casing_box_y/2., Ecal_casing_box_z/2.);
  Box solidCaliceECalCasingCut(Ecal_casing_box_x/2. - Ecal_casing_thickness,
                               Ecal_casing_box_y/2. - Ecal_casing_thickness,
                               Ecal_casing_box_z/2.);

  dd4hep::SubtractionSolid solidCaliceECalCasing("solidCaliceECalCasing",
                                         solidCaliceECalCasing1, solidCaliceECalCasingCut,
                                         Position(0., 0., Ecal_casing_thickness));

  Volume logicCaliceECalCasing("logicCaliceECalCasing", solidCaliceECalCasing, ecalCasingMaterial);

//   Container
  Box    solidCaliceECalContainer(Ecal_casing_box_x/2.0, Ecal_casing_box_y/2.0, Ecal_casing_box_z/2.0);
  Volume logicCaliceECalContainer("logicCaliceECalContainer", solidCaliceECalContainer, envMaterial);

  PlacedVolume phvCaliceECalContainer = logicCaliceECalContainer.placeVolume(logicCaliceECalCasing, Position(0.0, 0.0, 0.0));
  if (OverlapTest) phvCaliceECalContainer.ptr()->CheckOverlaps();

  return logicCaliceECalContainer;
}



void ConstrucBSMSupportAssembly(dd4hep::Detector& description, Volume& envvol, const double shift_x,
                                const double Calo_dim_y, const double Calo_dim_z)
{
  double ypestal         = 1.0*dd4hep::m;
  double Ecal_casing_thickness = description.constant<double>("Calice_Ecal_casing_thickness");
  double Ecal_casing_box_x = description.constant<double>("Calice_Ecal_casing_box_x");
  double FloorSurfaceYpos  = description.constant<double>("FloorSurfaceYpos");
  double BSMCaloZPos       = description.constant<double>("BSMCaloZPos");
  bool   OverlapTest       = (description.constant<int>("OverlapTest") != 0);
  double Ecal_casing_box_y = Calo_dim_y + 2.0 * Ecal_casing_thickness;
  double Ecal_casing_box_z = Calo_dim_z + Ecal_casing_thickness;

  double shift_z = BSMCaloZPos + 0.5*Ecal_casing_box_z;

  double ylevel   = Ecal_casing_box_y/2.;
  double tblhight = -ylevel - ypestal - FloorSurfaceYpos;
  double tblx     = 1.1*Ecal_casing_box_x;
  double tblz     = 1.1*Ecal_casing_box_z;

  dd4hep::Assembly tablesupport = LxAux::BuildTable(description, "BSMCaloTable", tblx, tblhight, tblz, 3);
  PlacedVolume pvTable = envvol.placeVolume(tablesupport, Position(shift_x, -ylevel, shift_z));
  if (OverlapTest) pvTable.ptr()->CheckOverlaps();

  Volume pedestal = LxAux::BuildPedestal(description, "BSMCalo",
                                       1.3*Ecal_casing_box_x, ypestal, 1.3*Ecal_casing_box_z);
  PlacedVolume pvPed = envvol.placeVolume(pedestal,
    Position(shift_x, 0.5*ypestal + FloorSurfaceYpos, shift_z));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();
}

static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens) {

  xml_det_t x_det = element;
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());

  Layering layering(x_det);
  xml_dim_t dim = x_det.dimensions();

  // --- create an envelope volume and position it into the world ---------------------

//   Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, element, sdet);
  dd4hep::Assembly   envelope(det_name + "_assembly");
  Volume fLogicWorld = theDetector.worldVolume();
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, dd4hep::Transform3D());
  sdet.setPlacement(envPV);

  dd4hep::xml::setDetectorTypeFlag(element, sdet);

  if (theDetector.buildType() == BUILD_ENVELOPE)
    return sdet;

  //-----------------------------------------------------------------------------------

//   Material air = theDetector.air();
  Material air = theDetector.material(theDetector.constant<std::string>("EnvironmentMaterial"));

  sens.setType("calorimeter");

  //====================================================================
  //
  // Read all the dimensions from compact.xml, user can update the value.
  // Use them to build a calo box prototye.
  //
  //====================================================================

  double Calo_dim_x = dim.x();
  double Calo_dim_y = dim.y();
  double Calo_dim_z = dim.z();

  printout(dd4hep::DEBUG, "building SamplingCaloBoxPrototype_v01",
           "Calo_dim_x : %e    Calo_dim_y: %e    Calo_dim_z: %e ", Calo_dim_x, Calo_dim_y, Calo_dim_z);

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();

  // unused: std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  // unused: double cell_sizeX      = cellSizeVector[0];
  // unused: double cell_sizeY      = cellSizeVector[1];

  // check if we have a TiledLayerGridXY segmentation :
  dd4hep::DDSegmentation::TiledLayerGridXY* tileSeg =
      dynamic_cast<dd4hep::DDSegmentation::TiledLayerGridXY*>(seg.segmentation());

  // access the layer identifier via the segmentation.
  // the layer identifier is defined in the compact xml file.
  // and use it to set the volumeID layer value later here.
  string identifierLayer = tileSeg->fieldNameLayer(); // "K" or "layer" or "..."

  //====================================================================
  //
  // general calculated parameters
  //
  //====================================================================

  // calorimeter dimensions
  double cal_hx = Calo_dim_x / 2.0;
  double cal_hy = Calo_dim_y / 2.0;
  double cal_hz = Calo_dim_z / 2.0;

  bool   OverlapTest       = (theDetector.constant<int>("OverlapTest") != 0);
  Material envMaterial = theDetector.material(theDetector.constant<std::string>("EnvironmentMaterial"));
  Box solidCaloBox(cal_hx, cal_hy, cal_hz);
  Volume logicCaloBox("logicCaliceStackContainer", solidCaloBox, envMaterial);

  //====================================================================
  //
  // build sampling layers in the CaloBox
  //
  //====================================================================

  int layer_num = 0;
  int layerType = 0;

  double layer_pos_z = -cal_hz;

  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    int repeat = x_layer.repeat();                // Get number of times to repeat this layer.
    const Layer* lay = layering.layer(layer_num); // Get the layer from the layering engine.
    double layer_thickness = lay->thickness();
    string layer_type_name = _toString(layerType, "layerType%d");

    // Loop over repeats for this layer.
    for (int j = 0; j < repeat; j++) {
      string layer_name = _toString(layer_num, "layer%d");
      DetElement layer(layer_name, layer_num);

      // Layer box & volume
      Volume layer_vol(layer_type_name, Box(cal_hx, cal_hy, layer_thickness / 2), air);

      // Create the slices (sublayers) within the layer.
      double slice_pos_z = -(layer_thickness / 2);
      int slice_number = 0;

      for (xml_coll_t k(x_layer, _U(slice)); k; ++k) {
        xml_comp_t x_slice = k;
        string slice_name = _toString(slice_number, "slice%d");
        double slice_thickness = x_slice.thickness();
        Material slice_material = theDetector.material(x_slice.materialStr());
        DetElement slice(layer, slice_name, slice_number);

        slice_pos_z += slice_thickness / 2;
        // Slice volume & box
        Volume slice_vol(slice_name, Box(cal_hx, cal_hy, slice_thickness / 2), slice_material);

        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }

        // Set region, limitset, and vis.
        slice_vol.setAttributes(theDetector, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        // slice PlacedVolume
        PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
        if (OverlapTest) slice_phv.ptr()->CheckOverlaps();
        slice_phv.addPhysVolID("slice", slice_number);
        slice.setPlacement(slice_phv);

        // Increment Z position for next slice.
        slice_pos_z += slice_thickness / 2;
        // Increment slice number.
        ++slice_number;
      }

      // Set region, limitset, and vis.
      layer_vol.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Layer position in Z within the stave.
      layer_pos_z += layer_thickness / 2;
      // Layer physical volume.
      PlacedVolume layer_phv = logicCaloBox.placeVolume(layer_vol, Position(0, 0, layer_pos_z));
      if (OverlapTest) layer_phv.ptr()->CheckOverlaps();
      layer_phv.addPhysVolID(identifierLayer, layer_num);
      layer.setPlacement(layer_phv);

      // Increment the layer Z position.
      layer_pos_z += layer_thickness / 2;
      // Increment the layer number.
      ++layer_num;
    }

    ++layerType;
  }

  Volume logicCaliceECalContainer = CaliceConstructCasing(theDetector, Calo_dim_y, Calo_dim_z);
  double shift_x = 0.0;
  double Ecal_casing_thickness = theDetector.constant<double>("Calice_Ecal_casing_thickness");
  dd4hep::RotationZYX crrot(0.0, M_PI, 0.0);
  PlacedVolume ecalstack_phv = logicCaliceECalContainer.placeVolume(logicCaloBox,
                                  dd4hep::Transform3D(crrot, Position(-shift_x, 0.0, 0.5*Ecal_casing_thickness)));

  double Ecal_casing_box_z   = Calo_dim_z + Ecal_casing_thickness;
  double BSMCaloZPos = theDetector.constant<double>("BSMCaloZPos");
  double shift_z = BSMCaloZPos + 0.5*Ecal_casing_box_z;

  PlacedVolume caloBox_phv = envelope.placeVolume(logicCaliceECalContainer,
                                                  dd4hep::Transform3D(crrot, Position(shift_x, 0, shift_z)));
  if (OverlapTest) caloBox_phv.ptr()->CheckOverlaps();
//   caloBox_phv.addPhysVolID("CaliceECal", 0);

  ConstrucBSMSupportAssembly(theDetector, envelope, shift_x, Calo_dim_y, Calo_dim_z);

  return sdet;
}

DECLARE_DETELEMENT(CaloPrototype_v02, create_detector)


