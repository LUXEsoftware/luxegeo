#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  string det_name = x_det.nameStr();
  Material air = description.air();
  DetElement sdet(det_name, x_det.id());
  Assembly assembly(det_name + "_assembly");
  PlacedVolume pv;
  sens.setType("tracker");

  int n = 0;
  for (xml_coll_t i(x_det, _U(layer)); i; ++i, ++n) {
    xml_comp_t x_layer = i;
    string l_name = det_name + _toString(n, "_layer%d");
    double x = x_layer.x();
    double y = x_layer.y();
    double z = x_layer.z();
    double dx = x_layer.dx();
    double dy = x_layer.dy();
    double dz = x_layer.dz();

    DetElement layer(sdet, _toString(n, "layer%d"), x_layer.id());
    Box l_box(0.5 * dx, 0.5 * dy, 0.5 * dz);
    Volume l_vol(l_name, l_box, air);

    int im = 0, nmod = 0;
    for (xml_coll_t j(x_layer, _U(stave)); j; ++j, ++im) {
      xml_comp_t x_stave = j;
      Material mat = description.material(x_stave.materialStr());
      string stave_name = l_name + _toString(im, "_stave%d");
      double thickness = x_stave.thickness();
      double length = x_stave.length();
      double width = x_stave.width();
      double xoffset = x_stave.x_offset();
      double zoffset = x_stave.z_offset();

      Box stave_box(0.5 * length, 0.5 * width, 0.5 * thickness);
      Volume stave_vol(stave_name, stave_box, mat);
      stave_vol.setAttributes(description, x_stave.regionStr(), x_stave.limitsStr(), x_stave.visStr());
      pv = l_vol.placeVolume(stave_vol, Transform3D(Position(xoffset, 0.0, zoffset)));

      int nModules = x_stave.nmodules();
      double gap = x_stave.gap();
      xml_comp_t x_sensor = x_layer.child(_U(sensor));
      Material sen_mat = description.material(x_sensor.materialStr());
      double sen_thickness = x_sensor.thickness();
      double sen_length = x_sensor.length();
      double sen_width = x_sensor.width();

      Assembly moduleAssembly("module");

      for (int imod = 0; imod < nModules; imod++) {
        nmod++;
        string sen_name = l_name + _toString(nmod, "_sensor%d");
        Box sen_box(0.5 * sen_length, 0.5 * sen_width, 0.5 * sen_thickness);
        Volume sen_vol(sen_name, sen_box, sen_mat);
        if (x_sensor.isSensitive()) {
          sen_vol.setSensitiveDetector(sens);
        }
        sen_vol.setAttributes(description, x_sensor.regionStr(), x_sensor.limitsStr(), x_sensor.visStr());

        DetElement sensorDE(layer, sen_name, nmod);
        pv = moduleAssembly.placeVolume(
            sen_vol, Transform3D(Position(-0.5 * length + 0.5 * sen_length + imod * (sen_length + gap), 0.0, 0.0)));
        pv.addPhysVolID("sensor", nmod);
        sensorDE.setPlacement(pv);
      }
      pv = l_vol.placeVolume(moduleAssembly,
                             Transform3D(Position(xoffset, 0.0, zoffset - 0.5 * sen_thickness - 0.5 * thickness)));
    }
    l_vol.setVisAttributes(description, x_layer.visStr());
    assembly.setVisAttributes(description, x_det.visStr());
    pv = assembly.placeVolume(l_vol, Transform3D(Position(x, y, z)));
    pv.addPhysVolID("layer", n);
    layer.setPlacement(pv);
  }

  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", sdet.id());
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(LUXETracker, create_detector)
