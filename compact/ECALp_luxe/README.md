LUXE ECAL-P - DD4hep
====================

## Structure
This directory contains the definition of the detector model geometries and configurations.
The following descriptions used by all detectors are stored in `..` directory:
* chemical elements,
* materials and extra materials,
* common definitions on the geometry of ECALs, and
* common definitions on the visulisation of ECALs.

The names of compact files for LUXE ECAL-P are formatted as `ECALp_luxe_{variant}.xml` and put under directory `ECALp_luxe`.

The detector consists of multiple layers that are stacked along the longitudinal _z_ direction.
Each layer is then further divided into several slices, which includes a tungsten plate, air gaps, and several silicon and electronic modules.
A module can be moved horizontally within the boundary of a layer,
and more than one module can be placed at the same _z_ position.

Current variant `ECALp_luxe_2512.xml`, partially based on the updated design for ECAL-P TB2025.
- Detector type `LUXEEcal`, id `2000`.
- 20 layers/20 X0 in total.
- Layer composition:
    - Pure tungsten (3.55 mm),
    - Compact silicon sandwich (CSIS), which consists of
        - Carbon fibre support (0.25 mm, to be confirmed [TBC]),
        - Regular glue (0.05 mm air, TBC),
        - Fanout film (0.1 mm copper, TBC),
        - Conductive glue (0.05 mm air, TBC),
        - Si wafer (0.32 mm),
        - Conductive glue (0.05 mm air, TBC),
        - High-voltage film (0.1 mm copper, TBC), and
    - Air gap to make 4.5 mm longitudinal pitch.
- Per-slice module placement to implement the dead areas between CSISs.
    - This variant has 6 CSISs.


_**TODO**_
- The readouts in the compact files are not properly implemented.
- The historic definitions should be cleaned out.


## Visualization
To test and visualize the detector, run
```shell
geoDisplay detector_compactfile.xml # or any other geometry
```
To check materials, distances and possible overlappings, run
```shell
materialScan detector_compactfile.xml x0 y0 z0 x1 y1 z1
```
It will display a list of materials encountered by a straight line from `(x0,y0,z0)` to `(x1,y1,z1)`.


## Contribution
This subrepository is based on [Gyros (permission required)](https://gitlab.desy.de/luxe-ecal/performance/gyros/-/tree/main/generation/geometry) by the LUXE ECAL group.
