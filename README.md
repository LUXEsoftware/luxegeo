# luxegeo

[![key4hep](https://github.com/LUXEsoftware/luxegeo/actions/workflows/key4hep.yml/badge.svg)](https://github.com/LUXEsoftware/luxegeo/actions/workflows/key4hep.yml)

Implementation of the detector geometry of the LUXE experiment using [DD4hep](https://dd4hep.web.cern.ch/dd4hep/).

## Requirements
luxegeo needs to be built against DD4hep with Geant4 enabled and a c++ compiler
supporting c++17. It is tested regularly in a continuous integration (CI)
workflow against [Key4hep](https://key4hep.web.cern.ch).

## Download and Installation
### Download
To get the latest version of luxegeo simply clone this repository

```bash
git clone https://github.com/LUXESoftware/luxegeo
```

or download and unpack a release tarball from github (see [Releases](https://github.com/LUXEsoftware/luxegeo/releases)).

### Setting up the dependencies
There are several ways of setting up the necessary dependencies for building and
installing luxegeo, the main is to use a Key4hep software stack, which requires
a `/cvmfs/` mount and CentOS7.

#### Using a Key4hep software stack
A Key4hep software stack comes with all dependencies fulfilled. Hence, sourcing
the setup script is enough, e.g.

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
```

Replacing `sw.hsf.org` by `sw-nightlies.hsf.org` will setup the bleeding edge
version of the Key4hep stack that is built nightly. This setup is less stable,
but might offer some fixes that have not yet landed in the proper releases on
`sw.hsf.org`.

<details>
  <summary>Instructions for other options</summary>

#### Using an LCG software stack
Similar to the Key4hep stack it is possible to use an LCG software stack that
also comes with all depenencies fulfilled. These releases are available or
different linux flavors and compilers, and again sourcing the approriate setup
script is enough, e.g.

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
```

#### Using an existing DD4hep installation
If you have an existing installation of DD4hep on your system, it should be
enough to initialize that via

```bash
source </path/to/dd4hep/installation>/bin/thisdd4hep.sh
```

</details>

### Building and installing luxegeo
luxegeo is a CMake project, so building it follows the usual steps

```bash
cd luxegeo
mkdir build && cd build
cmake -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_INSTALL_PREFIX=../install \
  ..
make install
```

luxegeo comes with a few tests to verify that the build and installation process
has been successful, they can be run via (in the `build` folder)

```bash
ctest --output-on-failure
```

### Running a simple example using a particle gun
After luxegeo has been installed, it is possible to run a simple simulation. In
the following we assume that you are in the `luxegeo` base directory, and that
the build artifacts are installed to `install` (this is what you get from
following the [build instructions](#building-and-installing-luxegeo)) The first
step is to setup luxegeo.

```bash
source install/bin/thisluxegeo.sh
export luxegeo_DIR=$(pwd)/install/share/luxegeo/compact
```

We then use `ddsim` to run a simulation using the LUXE tracker. The example uses
a particle gun to shoot 100 positrons in one shot and produces an output in
EDM4hep format. (Note, that you need DD4hep that has been built with EDM4hep
enabled)

```bash
ddsim --compactFile ${luxegeo_DIR}/LUXE/Tracker/LUXETracker.xml \
  --numberOfEvents 1 \
  --enableGun \
  --gun.multiplicity 100 \
  --gun.particle e+ \
  --outputFile positrons_tracker_edm4hep.root
```
