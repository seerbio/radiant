#!/bin/bash
#
# build_deb.sh
#
# Run from a directory containing an already-built version of Pythia,
# e.g. the `/app/` directory of the container built by `Dockerfile`.
#
# You MUST define `$package_dir`, which is the location where the
# contents of the DEB will be assembled. They will not be cleaned up.
#
# The created DEB should be located at `${package_dir}.deb`.

# "Strict mode"
set -eu -o pipefail

# Required argument (as env var)
package_dir=${package_dir}

pythia_bin=${pythia_bin:=${package_dir}/usr/local/bin/PythiaDIACpp}
pythia_lib=${pythia_lib:=${package_dir}/usr/lib}

mkdir -p ${pythia_bin}
mkdir -p ${pythia_lib}

echo "Copying DEB contents to package_dir…"

cp bin/IstrosLibraryBuilder $pythia_bin/IstrosLibraryBuilder
cp bin/LycaonMzMLTransmorgitron $pythia_bin/LycaonMzMLTransmorgitron
cp bin/OrpheusFastaShredder $pythia_bin/OrpheusFastaShredder
cp bin/PythiaDIA $pythia_bin/PythiaDIA
cp bin/rnn_linear_charge_w_precursors_nce_1.hdf5.json $pythia_bin/rnn_linear_charge_w_precursors_nce_1.hdf5.json
cp bin/rnn_linear_charge_w_precursors_nce_2.hdf5.json $pythia_bin/rnn_linear_charge_w_precursors_nce_2.hdf5.json
cp bin/rnn_linear_charge_w_precursors_nce_3.hdf5.json $pythia_bin/rnn_linear_charge_w_precursors_nce_3.hdf5.json
cp bin/rnn_linear_charge_w_precursors_nce_4.hdf5.json $pythia_bin/rnn_linear_charge_w_precursors_nce_4.hdf5.json
cp AlgorithmsLib/libAlgorithmsLib.so $pythia_lib/libAlgorithmsLib.so
cp ChemLib/libChemLib.so $pythia_lib/libChemLib.so
cp EigenLib/libEigenLib.so $pythia_lib/libEigenLib.so
cp FileReadersLib/libFileReadersLib.so $pythia_lib/libFileReadersLib.so
cp MachineLrnLib/libMachineLrnLib.so $pythia_lib/libMachineLrnLib.so
cp UtilsLib/libUtilsLib.so $pythia_lib/libUtilsLib.so
cp WorkFlowsLib/libWorkFlowsLib.so $pythia_lib/libWorkFlowsLib.so
cp /usr/lib/libarrow.so.1100 $pythia_lib/libarrow.so.1100
cp /usr/lib/libparquet.so.1100 $pythia_lib/libparquet.so.1100

mkdir -p ${package_dir}/DEBIAN/
cp control $package_dir/DEBIAN/control

echo "Building PythiaDIA DEB…"

dpkg-deb --build $package_dir
