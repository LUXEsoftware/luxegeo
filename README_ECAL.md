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

The detector consists of multiple layers that are stacked along the longitudinal _z_ direction.
Each layer is then further divided into several slices, which includes a tungsten plate, air gaps, and several silicon and electronic modules.
A module can be moved horizontally within the boundary of a layer,
and more than one module can be placed at the same _z_ position.

#### `ECALe_luxe` compact files
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

For previous versions and variants:
- v2 (`ECALe_luxe_v2.xml`):
  - Detector type `LUXEEcal`, id `2500`.
  - Per-slice module placement (8 modules per slice) with dead areas
  - layer spacing reduced to 10 mm via `Ecal_w_slab_gap500_10mmgap`.
- v1 (`ECALe_luxe_v1.xml`):
  - Detector type `LUXEEcal`, 15 layers
  - Tungsten per layer: 2 × 2.1 mm; gap constant `Ecal_w_slab_gap500_4p2mmW`
  - Envelope `z = Ecal_dim_z/2.0`
  - v1 variants:
    - `ECALe_luxe_v1_2p1mm.xml`: 1 × 2.1 mm per layer; uses `Ecal_w_slab_gap500_2p1mmW`
    - `ECALe_luxe_v1_5x4p2mm_10x2p1mm.xml`: first 5 layers at 2 × 2.1 mm (`_4p2mmW`), next 10 layers at 1 × 2.1 mm (`_2p1mmW`)
    - `ECALe_luxe_v1_10x2p1mm_5x4p2mm.xml`: first 10 layers at 1 × 2.1 mm (`_2p1mmW`), next 5 layers at 2 × 2.1 mm (`_4p2mmW`)
    - `ECALe_luxe_v1_compact.xml`: consists of only W and Si layers
- v0 (`ECALe_luxe_v0.xml`):
  - Initial compact file of ECALe_luxe derives from the CALICE `SiWECAL-TB202206` TB prototype.
  - Detector type `CaloPrototype_v02`, 15 layers
  - Tungsten per layer: 3 × 2.1 mm; gap constant `Ecal_w_slab_gap500`
  - Envelope `z` set numerically (`14.25`)

_**IMPORTANT**_ Until v1, there is no pixelization described in the compact file. Namely, each layer corresponds to one single active area.
The pixelization needs to be done posteriorly (e.g. via `PixelizationProcessor`).
In v2, the dead areas, where there are no silicon, is implemented.
But still, the pixelization inside a module needs to be done posteriorly.

#### `ECALp_luxe` compact files
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


### `src` detector libraries
The definition of dd4hep detector type.
The dd4hep framework requires a library to understand how the corresponding detector is constructed.
For `ECALe_luxe_v0` and `v1`, the type of detector `CaloPrototype_v02` are used.
This detector type has been implemented as a part of dd4hep in [`k4geo/detector/CaloTB`](https://github.com/key4hep/k4geo/blob/main/detector/CaloTB).
The detector consists of multiple layers, which are composited by slices of different materials.

From `ECALe_luxe_v2` and beyond, a customised library, as a part of the future `luxegeo`, is built based on `CaloPrototype_v02`.
It enables subunits called modules inside a layer slice.
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
