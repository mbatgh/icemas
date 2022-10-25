
# icemas

Icemas (Integration of the Classical Equations of motion of Molecules of Arbitrary Shape) is a 
C++ code to perform molecular dynamics (MD) simulations. The program has most of the basic features to be found
in any MD package. If you look for high performance and/or a number of special features icemas is not your choice.
Quite a few free software packages are available nowadays (e.g. NAMD, Gromacs, OpenMM, etc) most of which
include functionalities and features that Icemas does not have. Icemas does, however, have the advantage of a very tidy and
clearly structured source code (achieved by extensive usage of the C++ concept of operator overloading). It therefore might
be useful for pedagogical purposes, i.e., for people that are new to the field and want to have a look at a readable source
code. It could also be helpful as a basis for testing new algorithms, as the structure of the code is transparent, so that it
is easy to make modifications, without causing undesired effects.

## Features

- energy expression including Lennard-Jones, Coulomb, bond, angle, and dihedral angle terms
- Electrostatic long range interactions via Ewald summation
- NVE and NVT (Nose Hoover) ensembles
- SHAKE
- neighbour-list

## Prerequisites

- A C++ compiler (recent versions of gcc on most Linux distros should work)

## Installation

On Linux (tested on Debian and Ubuntu, but should work on most distros)
```
git clone https://github.com/mbatgh/icemas.git
cd icemas/src
make
sudo make install
```
The binary *icemas* should now be in /usr/local/bin

## Usage

Have a look at the content of the [examples](examples) folder for a quick start.

Generally three input-files are required:

1. A parameter-file called *pars.name*, where *name* is an arbitrary ID chosen by the user.
A default-parameter-file with short explanations of all possible entries are explained can be found
in the [examples](examples) folder.
2. A system-file (name of which must be given in the parameter-file). It contains information such
as number and type of molecules and force-field-parameters. Examples can be found in the [system](system) folder.
3. a file (name provided in the parameter-file) containing the initial-coordinates. Examples, and a tool
to generate such files are provided in the [coords](coords) folder.

Given these files the simulation is started on the command-line, by typing
```
icemas name [name] ...
```
If more than one names are provided as arguments icemas will do each of the named runs in succession.

Depending on the data in the parameter-file, Icemas will
write several output-files in the diectory, it was started in.

### *data.name*
    this file contains the thermodynamic data. In order of columns:
      - simulation-time [ps]
      - temperature [K]
      - Pressure [MPa]
      - Virial [kJ/mol]
      - kinetic Energy [kJ/mol]
      - pot. Energy (VdW-interactions) [kJ/mol]
      - pot. Energy (real part of Ewald-sum) [kJ/mol]
      - pot. Energy (imag. part of Ewald-sum) [kJ/mol]
      - pot. Energy (due to bond-angles) [kJ/mol]
      - pot. Energy (due to dihedral-angles) [kJ/mol]
      - pot. Energy (due to improper-dihedral-angles) [kJ/mol]    
      - total Energy [kJ/mol]

### *r_out.name*
The final coordinates

### *angdis.name, dhangdis.name, idhangdis.name*
Files, containig distribution-functions for the (proper dihedral,
improper dihedral) angles.

### *phsp.name*
    this file contains data quantifying the extent of a possible phase-separation. 

### *pcf.N.name*
Pair-correlation-functions, where N is a number (1 up to the
total number of pairs specified in 'pars.name')

## Limitations

- Runs on CPU only (no GPU support)
- No automated generation of force field parameters (writing a simple script to transform topology files
  from other packages, such as Gromacs, Charmm, or Amber should be straight forward, though)

## Contact

Please send any questions/bug-reports/suggestions to me. My full name is rare enough, so that,
together with the place where I live (Graz), google will promptly provide you with my email address.

