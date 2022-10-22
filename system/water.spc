# water system file
#

   1     # Nr. of different species' of molecules
   2     # Nr. of different Lennard-Jones-types (including "NO" LJ-type !)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# for the first species:
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# nmb. of molecules  
# nmb. of sites (per molec.)  
# nmb. of angles
# nmb. of dihedral angles
# nmb. of improper dihedral angles
# 
# NoM NoS NoB NoA NoD NoI
#

   330  3   3   0   0   0

#
# for each LJ-type
#
# Number  - epsilon [kJ/mol] - sigma [nm]
#

   0  0.0      0.0
   1  0.65017  0.31656

#
# for each species :
# 
# for each site:
#
# number  LJ-type   mass [amu]  charge [e0]
#

   0  0  1.0   0.4238
   1  1  16.0  -0.8476
   2  0  1.0   0.4238

#
# for each bond
#
# number index1 index2 bondlenght [nm]
#

   0  0  1  0.1   
   1  1  2  0.1
   2  0  2  0.1632981

#
# for each angle
# 
# number index1 index2 indxe3 angle force 
# 

   
# 
# for each dihedral angle
#
# number index1 index2 indxe3 index4 angle force
#

#
# if there is more than one type of molecule
# data for the second species of molecules follows here:
# 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# for second species:
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ...
