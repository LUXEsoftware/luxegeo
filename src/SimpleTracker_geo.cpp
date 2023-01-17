#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  Material   air       = description.air();
  DetElement sdet        (det_name,x_det.id());
  Assembly   assembly    (det_name+"_assembly");
  PlacedVolume pv;
  int n = 0;
  sens.setType("tracker");

  for(xml_coll_t i(x_det,_U(layer)); i; ++i, ++n)  {
    xml_comp_t x_layer = i;
    string  l_name = det_name+_toString(n,"_layer%d");
    double  x    = x_layer.x();
    double  y    = x_layer.y();
    double  z    = x_layer.z();
    double  dx    = x_layer.dx();
    double  dy    = x_layer.dy();
    double  dz    = x_layer.dz();

    DetElement layer(sdet,_toString(n,"layer%d"),x_layer.id());
    Box    l_box (0.5*dx, 0.5*dy, 0.5*dz);
    Volume  l_vol(l_name,l_box,air);
    int im = 0;

    for(xml_coll_t j(x_layer,_U(slice)); j; ++j, ++im)  {
      xml_comp_t x_slice = j;
      Material mat = description.material(x_slice.materialStr());
      string s_name= l_name+_toString(im,"_slice%d");
      double thickness = x_slice.thickness();
      double length = x_slice.length();
      double width = x_slice.width();

      Box   s_box(0.5*length, 0.5*width, 0.5*thickness);
      Volume s_vol(s_name, s_box, mat);

      if ( x_slice.isSensitive() ) {
        s_vol.setSensitiveDetector(sens);
      }
      // Set Attributes
      s_vol.setAttributes(description,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
      pv = l_vol.placeVolume(s_vol);
      // Slices have no extra id. Take the ID of the layer!
      pv.addPhysVolID("slice",im);
    }
    l_vol.setVisAttributes(description,x_layer.visStr());
      
    pv = assembly.placeVolume(l_vol, Transform3D(Position(x,y,z)));
    pv.addPhysVolID("layer",n);
    layer.setPlacement(pv);
  }

  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system",sdet.id());
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(SimpleTracker,create_detector)
