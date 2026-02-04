LUXE ECAL-E - DD4hep
====================

## Structure
This directory contains the definition of the detector model geometries and configurations.
The following descriptions used by all detectors are stored in `..` directory:
* chemical elements,
* materials and extra materials,
* common definitions on the geometry of ECALs, and
* common definitions on the visulisation of ECALs.

The names of compact files for LUXE ECAL-E are formatted as `ECALe_luxe_{variant}.xml` and put under directory `ECALe_luxe`.

The detector consists of multiple layers that are stacked along the longitudinal _z_ direction.
Each layer is then further divided into several slices, which includes a tungsten plate, air gaps, and several silicon and electronic modules.
A module can be moved horizontally within the boundary of a layer,
and more than one module can be placed at the same _z_ position.

Current variant `ECALe_luxe_2512.xml`, partially based on this [interim report for DRD-Calo SiW ECAL](https://indico.cern.ch/event/1551941/contributions/6683124).
- Detector type `LUXEEcal`, id `2500`.
- 15 layers/18 X0 in total.
- Layer composition:
    - Pure tungsten (4.2 mm),
    - Carbon fibre support (1.5 mm, to be confirmed [TBC]),
    - Active sensor unit (ASU), which consists of
        - High-voltage film (0.1 mm copper, TBC),
        - Conductive glue (0.1 mm air, TBC),
        - Si wafer (0.5 mm),
        - Conductive glue (0.1 mm air, TBC),
        - PCB board (1.7 mm),
        - On board ASIC (1.2 mm PCB, TBC), and
    - Air gap to make 15 mm longitudinal pitch.
- Per-slice module placement to implement the dead areas between ASUs and wafers.
    - This variant has 2 ASUs and 8 wafers.


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

For previous variants `ECALe_luxe_v{0-2}.xml`, please also check the said Gyros repository.
