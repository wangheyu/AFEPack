# FAST AFEPack examples

This directory contains a copy of the AFEPack examples from
`~/Projects/FEM/examples/afepack` and their EasyMesh inputs from
`~/Projects/FEM/meshes/afepack`.

The original FEM project is not used or modified by this Makefile.

## Layout

- `src/`: AFEPack example sources.
- `meshes/afepack/`: EasyMesh `.d` inputs.
- `build/`: generated meshes and binaries.

## Build

```sh
make -C examples/FAST
```

Useful targets:

```sh
make -C examples/FAST list
make -C examples/FAST info
make -C examples/FAST meshes
make -C examples/FAST poisson_convergence_afepack
make -C examples/FAST run-poisson_convergence_afepack
```

The Makefile defaults to the AFEPack library, headers, templates, and
EasyMesh executable built in this repository. Override these when needed:

```sh
make -C examples/FAST \
  AFEPACK_INCLUDE_DIR=$HOME/local/AFEPack/include \
  AFEPACK_LIB_DIR=$HOME/local/AFEPack/lib \
  AFEPACK_TEMPLATE_DIR=$HOME/local/AFEPack/include/AFEPack/template/triangle \
  EASYMESH=$HOME/local/AFEPack/bin/easymesh
```
