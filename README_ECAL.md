LUXE ECAL Gyros - Generation
============================
<!-- (lines omitted) -->

## Structure
### `compact`
This directory contains the definition of the detector model geometries and configurations. The following descriptions are used by all detectors:
* chemical elements,
* materials and extra materials,
* common definitions on the geometry of ECALs, and
* common definitions on the visulisation of ECALs.

Each detector should have its own subdirectory, where stores its description (called “compact file”) of various versions.
The name of compact files should have this format: `{name of detector}_{variation}.xml` and must be put in a directory called `{name of detector}`.


#### `ECALe_luxe` compact files
Summaries for different versions:
- v0: Initial compact file of ECALe_luxe derives from the CALICE `SiWECAL-TB202206`
- v1: W thickness changed to 4.2 mm ~~(previously 6.3 mm)~~
- v2: REQUIRE DETECTOR TYPE `LUXEEcal`. Dead area implemented between staves; layer distance reduced to 10 mm ~~(previously 15 mm)~~

_**IMPORTANT**_ Until v1, there is no pixelization described in the compact file. Namely, each layer corresponds to one single active area.
The pixelization needs to be done posteriorly (e.g. via `PixelizationProcessor`).
In v2, the dead areas, where there are no silicon, is implemented.
But still, the pixelization inside a stave needs to be done posteriorly.

### `src` detector libraries
The definition of dd4hep detector type.
The dd4hep framework requires a library to understand how the corresponding detector is constructed.
For `ECALe_luxe_v0` and `v1`, the type of detector `CaloPrototype_v02` are used.
This detector type has been implemented as a part of dd4hep in [`k4geo/detector/CaloTB`](https://github.com/key4hep/k4geo/blob/main/detector/CaloTB).
The detector consists of multiple layers, which are composited by slices of different materials.

From `ECALe_luxe_v2` and beyond, a customised library, as a part of the future `luxegeo`, is built based on `CaloPrototype_v02`.
It enables subunits called staves inside a layer slice.
Namely, some slices are no longer monolithic but have gaps inside, which represent the bordering space between sensor modules.
For more details, please refer to the [`luxegeo/README.md`](https://github.com/LUXEsoftware/luxegeo/blob/main/README.md) file.


<!-- (lines omitted) -->

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
This subrepository is partially based on the [work](https://gitlab.cern.ch/calice/calice_dd4heptestbeamsim/-/tree/master/2017_SiECAL_DESY/) of Daniel Jeans.
