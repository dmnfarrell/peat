# $Id: TOPOLOGY.H,v 1.51 2002/12/28 09:22:54 vriend Exp $
#
# This is the topology file for WHAT IF if you do want to use protons.
# First you give the amino acids.
# Thereafter DNA/RNA
# The file ends with the drugs for which a topology has been defined.
# Lines starting with a # are comments, and are neglected by WHAT IF
#
# Per amino acid the following information is given:
#  1) An asterisk
#  2) The name in 3 and 1 letter code (A4,A1 left justified)
#  3) Twelve integers (free format) that indicate:
#      1) AT   The number of atoms including protons
#      2) BND  The number of intra-residue bonds
#      3) CHI  The number of rotable side-chain torsion angles (torsion angles that the
#              user can rotate with the TORS option in the GRAFIC window)
#      4) TRS  The number of measurable side-chain torsion angles (total number of torsion
#              angles that gets listed in the SHOCHI option)
#      5) CRD1  -> Standard atomic coordinates are present ( 0 -> not ). Historical
#              artefact, you should always use 1, and give the coordinates.
#      6) #BL  Number of bondlengths (historical artefact, should be same as BND).
#      7) HGR  Number of hydrogen bond 'groups' (eg COO group in ASP or ring in HIS).
#      8) HAT  Number of heavy atoms (non-protons) in the residue.
#      9) ANG  Number of bond angles
#     10) 1-3  1 means 1-3 connections are given, 0 means they are not given.
#              (Historical artefact, always use 1 and give the 1-3 interactions)
#     11) PAR  Number of atomic parameter sets per class (max : 4)
#     12) PC   Number of atomic parameter classes (max : 2)
#     13) HANG Number of angles given under 9 that do not include protons
#     14) TYP  Residue family type:
#              1 = protein
#              3 = nucleic acid
#              6 = co-factor with topology
#  4) The atom names (14A4; Mendelev symbol in left two charaters; side chain depth and
#     branching indicators in characters 3 and 4)
#  5) Pairs of bonded atoms (10 pairs of integers per line, free format)
#     Use multiple lines if needed.
#  6) Per rotatable side-chain: the two atoms that define the axis followed
#     by all atoms that should move around if the users chooses to rotate around
#     this bond with the TORS option (20 integers per line, free format)
#  7) Per measurable torsion angle the four atoms (20 free format 
#     integers per line; that is 5 torsion angles per line)
#  8) Per atom the following information (AT lines):
#      1) Standard X coordinate \
#      2) Standard Y coordinate |-> free format within columns 1-29
#      3) Standard Z coordinate /
#      4) T or F for hydrogen bondable or not. Fixed in column 30.
#      5) + or - for donor or acceptor (blank means that this atom can either be
#         a donor or an acceptor, depending on the environment. Keep this character also
#         blank in case the atom cannot be involved in H-bonds) Fixed in column 31
#      6) Hydrogen bond donor type number. Fixed column 33.
#         1 = backbone nitrogen
#         2 = Arg Ne or Trp side chain N
#         3 = Asn side chain NH2
#         4 = Ser side chain Og
#         5 = Lysine side chain N
#         6 = Histidine side chain N
#         7 = Oxygen in Asp or Glu (only used when these residues are protonated, as these are donor
#             type indicators.
#      7) Number of hydrogen bond group this donor/acceptor is member of. Fixed column 35.
#      8) Numbers of the planes this atom is part of. Fixed column 36-39. In amino acids
#         and nucleic acids only column 36 will be needed because they all have at most
#         one planar group, but in a drug with two phenyl rings that are connected via a flexible 
#         linker, you need to use column 36 for the one ring and column 37 for the second
#      9) H angle type, fixed column 40 (not yet used)
#  9) Per bond length the standard bond lengths with sigma (one per line, free format; BND lines)
# 10) Bond angles with sigma, preceded by numbers of the atoms making the bond (ANG lines)
# 11) AT lines; one for each atom in the residue. Each line contains (free format integers):
#     If the number of 1-3 bonded atoms is X, you type here X+1
#     X numbers indicating those 1-3 bonded atoms
# 12) One line for each hydrogen bond group (I5,30F5.2):
#     1) Number of atoms in the hydrogen group 
#     2) Penalties for having 0...N hydrogens in the group.
# 13) PAR times PC lines with atomic parameter values. Each of these lines contains 
#     AT free format real numbers. At present parameter types 1 and 2 in parameter group 1 are
#     used to store partial atomic charges. I.e. 6 parameter sets are still free for use.
#
#============================================================================
#
#                  Asp as an annotated example
#                  Remove the first two charaters from each
#                  line in this example....
#
#============================================================================
# *                                                           This asterix indicates start of residue
# ASP D                                                       Three and one letter code A4,A1
# # AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP   Names of the 14 parameters
#    12  11  2   2   1   11  1   8  17   1   2   1  8    1    Values of the 14 parameters
#    |   |   |   |   |   |   |   |  |    |   |   |  |    \--->Residue type 1 (amino acid)
#    |   |   |   |   |   |   |   |  |    |   |   |  \-------->8 of the angles include only heavy atoms
#    |   |   |   |   |   |   |   |  |    |   |   \----------->Only one group of parameters is being used
#    |   |   |   |   |   |   |   |  |    |   \--------------->and this group holds two sets of params
#    |   |   |   |   |   |   |   |  |    \------------------->Must always be 1
#    |   |   |   |   |   |   |   |  \------------------------>17 bond angle values are given
#    |   |   |   |   |   |   |   \--------------------------->The first 8 atoms are non-proton
#    |   |   |   |   |   |   \------------------------------->Asp has 1 hydrogen bond group (COO)
#    |   |   |   |   |   \----------------------------------->Should be same as BND
#    |   |   |   |   \--------------------------------------->Should always be 1
#    |   |   |   \------------------------------------------->The side chain contains two torsion angles
#    |   |   |                                                that can be measured
#    |   |   \----------------------------------------------->The side chain contains two torsion angles
#    |   |                                                    that the user can rotate around with TORS
#    |   \--------------------------------------------------->ASP contains 11 bonds
#    \------------------------------------------------------->ASP has 12 atoms
# #===== Atom names
# #...||||....||||....||||....||||....||||....||||....||||    Marks the N A4 positions for atoms
#  N   CA  C   O   CB  CG  OD1 OD2 H   HA 2HB 3HB             The atoms. First heavy atoms, protons at the
#                                                             end. 14 atoms per line. Use multiple lines
#                                                             if needed.
# #===== Bonds
# 1 2  2 3  3 4  2 5  5 6  6 7  6 8  1 9  2 10  5 11          The atom pairs that are covalently connected.
# 5 12                                                        E.g. 1 2 means that N is bound to CA. Give
#                                                             maximally 10 pairs per line.Use multiple lines
#                                                             if needed (like the 5 and 12 here).
# #===== Rotatable bonds
#   2  5    6  7  8 11 12                                     This means the user can rotate around the bond
#                                                             2 5 which is the CA-CB bond. And if she does, the 
#                                                             atoms 6,7,8,11,12 being (CG.OD1,OD2,2HB,3HB) will
#                                                             move around.
#                                                             
#   5  6    7  8
# #===== Torsion angles
# 1 2 5 6   2 5 6 7                                           Two times four numbers indicate two torsion 
#                                                             angles that will be listed by SHOCHI. The first
#                                                             is the angle Chi-1 from N to CG via CA and CB, the
#                                                             second is Chi-2 from CA to OD1 via CB and CG.
# #===== Coordinates, hydrogen bonding and planarity
#   1.20131  0.84661  0.00001  T+ 1                           N
#   0.00001  0.00001  0.00001  F                              CA
#  -1.25031  0.88111  0.00001  F                              C
#  -2.18531  0.66031 -0.78411  T-                             O
#   0.04341 -1.01681  1.12341  F     1                        CB
#   1.22931 -1.96551  1.05271  F     1                        CG
#   2.05771 -1.86411  0.12721  T  7 11                        OD1
#   1.31911 -2.81041  1.95881  T  7 11                        OD2
#      |        |        |     || | |\----------------------->Planar group indicator (CB,CG,OD1,OD2 are planar)
#      |        |        |     || | \------------------------>OD2 sits in hydrogen bond group 1 together with OD1
#      |        |        |     || \-------------------------->If OD2 is protonated and donates, it is donortype 7
#      |        |        |     |\---------------------------->blank indicates that OD2 can donate and/or accept
#      |        \        \     \----------------------------->T indicates that OD2 can make Hbonds
#      \----------------------------------------------------->The atomic coordinates
#   1.20131  1.84661  0.00001  F         1                    H (the 1 at the end of the line isn't being used yet)
#   0.00002 -0.51093 -0.85961  F                              HA
#   0.08630 -0.52385  1.99241  F                              2HB
#  -0.79593 -1.55935  1.08923  F                              3HB
# #===== Standard bonds
# 1.458 0.019                                                 The first bond (1-2 is N-CA) is 1.458+/-0.019 A
# 1.525 0.021
# 1.231 0.020
# 1.530 0.020
# 1.516 0.025
# 1.249 0.019 
# 1.249 0.019
# 1.000 0.050
# 1.000 0.050
# 1.000 0.050
# 1.000 0.050
# #===== Standard angles
#  1  2  3  111.2  2.8                                        The first angle (1-2-3, N-CA-C) is 111.2+/-2.8 degrees
#  2  3  4  120.8  1.7
#  1  2  5  110.5  1.7
#  3  2  5  110.1  1.9
#  2  5  6  112.6  1.0
#  5  6  8  118.4  2.3
#  5  6  7  118.4  2.3
#  8  6  7  122.9  2.4
#   1  2 10  109.5    5.0                                     Angles involving protons aren't known too well yet,
#   2  1  9  120.0    5.0                                     so they get ideal values and large standard
#   2  5 11  109.5    5.0                                     deviations.
#   2  5 12  109.5    5.0
#   3  2 10  109.5    5.0
#   5  2 10  109.5    5.0
#   6  5 11  109.5    5.0
#   6  5 12  109.5    5.0
#  11  5 12  109.5    5.0
# #===== 1-3 connections
#    4   3   5  10                                            This is a bit stupid. The first number is 4 because
#    6   9   4   6  11  12                                    atom 1 (the backbone N) has three atoms that are
#    4   1   5  10                                            1-3 connected (plus one, see above). These atoms are
#    2   2                                                    3 (C), 5 (CB) and 10 (HA).
#    6   1   3  10   7   8
#    4   2  11  12
#    3   5   8
#    3   5   7
#    2   2
#    4   1   3   5
#    4   2   6  12
#    4   2   6  11
# #===== Hydrogen bond donor groups
#     2 0.00 0.75 9.99                                       Here one line per donor group. There is only one
#                                                            group in ASP (OD1,OD2). This group consists of 2 atoms.
#                                                            The group can contain 0,1,2 protons and the penalties
#                                                            are: no proton=0.0; 1 proton=0.75; 2 protons=9.99
#                                                            (9.99 means that you need to gain 9.99 ideal Hbonds
#                                                            to compensate for 2 protons to happen)
# #===== Parameter groups
#  0.1 0.0 0.38 -0.48 0.0 0.27 -0.635 -0.635                 First parameter set one real per atom
#  0.1 0.0 0.38 -0.48 0.0 0.53 -0.18  -0.18                  Second parameter set one real per atom
# #
#
#============================================================================
#
#                End of ASP example
#
#==1=========================================================================
*
ALA A   71.079
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  10  9   0   0   1   9   0   5  14   1   1   1   4    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  H   HA 1HB 2HB 3HB
#===== Bonds
1 2  2 3  3 4  2 5  1 6  2 7  5 8  5 9  5 10
#===== Rotatable bonds
#===== Torsion angles
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1      
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.02021 -0.92681  1.20941  F
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.55670 -0.83069  F
  0.83343 -1.50748  1.17080  F
  0.03968 -0.38117  2.04720  F
 -0.79922 -1.49996  1.20384  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.521 0.033
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.4  1.5
 3  2  5  110.5  1.5
  1  2  7  109.5    5.0
  2  1  6  120.0    5.0
  2  5  8  109.5    5.0
  2  5  9  109.5    5.0
  2  5 10  109.5    5.0
  3  2  7  109.5    5.0
  5  2  7  109.5    5.0
  8  5  9  109.5    5.0
  8  5 10  109.5    5.0
  9  5 10  109.5    5.0
#===== 1-3 connections
   4   3   5   7
   6   6   4   8   9  10
   4   1   5   7
   2   2
   4   1   3   7
   2   2
   4   1   3   5
   4   2   9  10
   4   2   8  10
   4   2   8   9
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0
#==2=========================================================================
*
CYS C  103.145
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
   11 10   1   1   1   10  0   6  15   1   1   1  5    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  SG  H   HA 2HB 3HB  HG
#===== Bonds
1 2  2 3  3 4  2 5  5 6  1 7  2 8  5 9 5 10 6 11 
#===== Rotatable bonds
  2  5    6  9 10 11
#===== Torsion angles
1 2 5 6
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.01801 -0.95831  1.18351  F
  1.43091 -2.08501  1.21211  T-
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.54187 -0.84045  F
  0.03295 -0.41677  2.02406  F
 -0.81812 -1.50609  1.15494  F
  2.01426 -1.90827  0.41935  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.808 0.033 1.822 0.020   (second is for S in S-S)
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.300 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  114.4  2.3 114.4  2.3 (second is for S in S-S)
  1  2  8  109.5    5.0
  2  1  7  120.0    5.0
  2  5  9  109.5    5.0
  2  5 10  109.5    5.0
  3  2  8  109.5    5.0
  5  2  8  109.5    5.0
  5  6 11  109.5    5.0
  6  5  9  109.5    5.0
  6  5 10  109.5    5.0
  9  5 10  109.5    5.0
#===== 1-3 connections
   4   3   5   8
   6   7   4   6   9  10
   4   1   5   8
   2   2
   5   1   3   8  11
   4   2   9  10
   2   2
   4   1   3   5
   4   2   6  10
   4   2   6   9
   2   5
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0
#==3=========================================================================
*
ASP D  115.104
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
   12  11  2   2   1   11  1   8  17   1   2   1  8    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  OD1 OD2 H   HA 2HB 3HB 
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8  1 9  2 10  5 11
5 12 
#===== Rotatable bonds
  2  5    6  7  8 11 12
  5  6    7  8
#===== Torsion angles
1 2 5 6   2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.04341 -1.01681  1.12341  F     1
  1.22931 -1.96551  1.05271  F     1
  2.05771 -1.86411  0.12721  T  7 11
  1.31911 -2.81041  1.95881  T  7 11
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.51093 -0.85961  F
  0.08630 -0.52385  1.99241  F
 -0.79593 -1.55935  1.08923  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.516 0.025
1.249 0.019 
1.249 0.019
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  112.6  1.0
 5  6  8  118.4  2.3
 5  6  7  118.4  2.3
 8  6  7  122.9  2.4
  1  2 10  109.5    5.0
  2  1  9  120.0    5.0
  2  5 11  109.5    5.0
  2  5 12  109.5    5.0
  3  2 10  109.5    5.0
  5  2 10  109.5    5.0
  6  5 11  109.5    5.0
  6  5 12  109.5    5.0
 11  5 12  109.5    5.0
#===== 1-3 connections
   4   3   5  10
   6   9   4   6  11  12
   4   1   5  10
   2   2
   6   1   3  10   7   8
   4   2  11  12
   3   5   8
   3   5   7
   2   2
   4   1   3   5
   4   2   6  12
   4   2   6  11
#===== Hydrogen bond donor groups
    2 0.00 0.75 9.99
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.27 -0.635 -0.635
 0.1 0.0 0.38 -0.48 0.0 0.53 -0.18  -0.18 
#==4=========================================================================
*
GLU E  129.116
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  15  14   3   3   1   14  1   9  23   1   2   1  9    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD  OE1 OE2 H   HA 2HB 3HB 2HG 
3HG 
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  7 8  7 9  1 10 2 11
5 12 5 13 6 14 6 15  
#===== Rotatable bonds
  2  5    6  7  8  9 12 13 14 15
  5  6    7  8  9 14 15
  6  7    8  9
#===== Torsion angles
1 2 5 6  2 5 6 7  5 6 7 8
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
 -0.05931 -0.87451  1.25831  F
  1.08171 -1.84981  1.50081  F     1
  0.98981 -2.68741  2.73981  F     1
  0.06741 -2.62371  3.53931  T  7 11
  1.97021 -3.44941  2.87091  T  7 11
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.58204 -0.81314  F
 -0.09876 -0.26025  2.04643  F
 -0.90302 -1.40913  1.21004  F
  1.12423 -2.46941  0.71704  F
  1.92852 -1.32051  1.55339  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.520 0.030
1.516 0.025
1.249 0.019 
1.249 0.019
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  114.1  2.0
 5  6  7  112.6  1.7
 6  7  9  118.4  2.3
 6  7  8  118.4  2.3
 9  7  8  122.9  2.4
  1  2 11  109.5    5.0
  2  1 10  120.0    5.0
  2  5 12  109.5    5.0
  2  5 13  109.5    5.0
  3  2 11  109.5    5.0
  5  2 11  109.5    5.0
  5  6 14  109.5    5.0
  5  6 15  109.5    5.0
  6  5 12  109.5    5.0
  6  5 13  109.5    5.0
  7  6 14  109.5    5.0
  7  6 15  109.5    5.0
 12  5 13  109.5    5.0
 14  6 15  109.5    5.0
#===== 1-3 connections
   4   3   5  11
   6  10   4   6  12  13
   4   1   5  11
   2   2
   7   1   3  11   7  14  15
   6   2  12  13   8   9
   4   5  14  15
   3   6   9
   3   6   8
   2   2
   4   1   3   5
   4   2   6  13
   4   2   6  12
   4   5   7  15
   4   5   7  14
#===== Hydrogen bond donor groups
    2 0.00 0.75 9.99
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.27 -0.635 -0.635
 0.1 0.0 0.38 -0.48 0.0 0.0 0.53 -0.18  -0.18 
#==5=========================================================================
*
PHE F  147.177
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  20  20   2   2   1  20   0  11  32   1   1   1  13    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  H   HA 2HB 
3HB  HD1 HD2 HE1 HE2 HZ
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8  7 9  8 10  9 11
10 11 1 12 2 13 5 14 5 15 7 16 8 17 9 18 10 19 11 20
#===== Rotatable bonds
  2  5    6  7  8  9 10 11 14 15 16 17 18 19 20
  5  6    7  8  9 10 11 16 17 18 19 20
#===== Torsion angles
1 2 5 6  2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.01741 -1.03521  1.11351  F     1
  1.15571 -2.01001  1.09071  F     1
  2.34071 -1.71821  1.76431  F     1
  1.04051 -3.22411  0.41651  F     1
  3.40411 -2.60301  1.76471  F     1
  2.09341 -4.13761  0.39651  F     1
  3.28821 -3.80771  1.06231  F     1
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.50342 -0.86403  F
  0.05168 -0.54592  1.98496  F
 -0.83339 -1.55777  1.05808  F
  2.42376 -0.85202  2.25707  F     1   1
  0.18936 -3.44351 -0.06037  F     1   1
  4.24464 -2.38455  2.26047  F     1   1
  2.00124 -5.01084 -0.08200  F     1   1
  4.06212 -4.44034  1.03327  F     1   1
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.502 0.023
1.384 0.021
1.384 0.021
1.382 0.030
1.382 0.030
1.382 0.030
1.382 0.030
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  113.8  1.0
 5  6  7  120.7  1.7
 5  6  8  120.7  1.7
 6  7  9  120.7  1.7
 6  8 10  120.7  1.7
 7  9 11  120.0  1.8
 8 10 11  120.0  1.8
 9 11 10  120.0  1.8
 7  6  8  120.0  1.8
  1  2 13  109.5    5.0
  2  1 12  120.0    5.0
  2  5 14  109.5    5.0
  2  5 15  109.5    5.0
  3  2 13  109.5    5.0
  5  2 13  109.5    5.0
  6  5 14  109.5    5.0
  6  5 15  109.5    5.0
  6  7 16  120.0    5.0
  6  8 17  120.0    5.0
  7  9 18  120.0    5.0
  8 10 19  120.0    5.0
  9  7 16  120.0    5.0
  9 11 20  120.0    5.0
 10  8 17  120.0    5.0
 10 11 20  120.0    5.0
 11  9 18  120.0    5.0
 11 10 19  120.0    5.0
 14  5 15  109.5    5.0
#===== 1-3 connections
   4   3   5  13
   6  12   4   6  14  15
   4   1   5  13
   2   2
   6   1   3  13   7   8
   8   2  14  15   9  16  10  17
   5   5   8  11  18
   5   5   7  11  19
   5   6  16  10  20
   5   6  17   9  20
   5   7  18   8  19
   2   2
   4   1   3   5
   4   2   6  15
   4   2   6  14
   3   6   9
   3   6  10
   3   7  11
   3   8  11
   3   9  10
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#==6=========================================================================
*
GLY G   57.052
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP 
   7   6   0   0   1   6   0   4   8   1   1   1  2    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   H  2HA 3HA
#===== Bonds
1 2  2 3  3 4  1 5  2 6  2 7
#===== Rotatable bonds
#===== Torsion angles
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.57742 -0.81643  F
  0.00002 -0.57742  0.81645  F
#===== Standard bonds
1.451 0.016
1.516 0.015
1.231 0.020
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  112.5  2.9
 2  3  4  120.8  2.1
  1  2  6  109.5    5.0
  1  2  7  109.5    5.0
  2  1  5  120.0    5.0
  3  2  6  109.5    5.0
  3  2  7  109.5    5.0
  6  2  7  109.5    5.0
#===== 1-3 connections
   4   3   6   7
   3   5   4
   4   1   6   7
   2   2
   2   2
   4   1   3   7
   4   1   3   6
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48
#==7=========================================================================
*
HIS H  137.141
# AT BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR  PC HNG TYP
  18  18   2   2   1  18   1  10  29   1   1   4  12   1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  ND1 CD2 CE1 NE2 H   HA 2HB 3HB 
 HD1 HD2 HE1 HE2
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8  7 9  8 10  9 10   
1 11 2 12 5 13 5 14 7 15 8 16 9 17 10 18
#===== Rotatable bonds
  2  5    6  7  8  9 10 13 14 15 16 17 18
  5  6    7  8  9 10 15 16 17 18
#===== Torsion angles
1 2 5 6  2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
#2345678901234567890123456789012345678901
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
 -0.06221 -0.96331  1.21101  F     1
  1.10221 -1.90821  1.26631  F     1
  2.27861 -1.64251  1.92131  T  6 11
  1.19481 -3.18861  0.83801  F     1
  3.04441 -2.71771  1.89061  F     1
  2.41061 -3.67351  1.24241  T  6 11
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.54656 -0.83741  F
 -0.07252 -0.42016  2.05058  F
 -0.90477 -1.49846  1.15006  F
  2.51778 -0.77286  2.35317  F     1   1
  0.50152 -3.68987  0.32025  F     1   1
  3.95645 -2.79400  2.29355  F     1   1
  2.75788 -4.59606  1.07417  F     1   1
#===== Standard bonds
1.458 0.019 
1.525 0.021 
1.231 0.020  
1.530 0.020 
1.497 0.014 
1.371 0.017 
1.356 0.011 
1.319 0.013
1.374 0.021 
1.345 0.020 
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  113.8  1.0
 5  6  7  121.6  1.5
 6  7  9  105.6  1.0
 7  9 10  111.7  1.3
 9 10  8  106.9  1.3
10  8  6  106.5  1.0
 8  6  7  106.1  1.0
 5  6  8  129.1  1.3
  1  2 12  109.5    5.0
  2  1 11  120.0    5.0
  2  5 13  109.5    5.0
  2  5 14  109.5    5.0
  3  2 12  109.5    5.0
  5  2 12  109.5    5.0
  6  5 13  109.5    5.0
  6  5 14  109.5    5.0
  6  7 15  125.0    5.0
  6  8 16  125.0    5.0
  7  9 17  125.0    5.0
  8 10 18  125.0    5.0
  9  7 15  125.0    5.0
  9 10 18  125.0    5.0
 10  8 16  125.0    5.0
 10  9 17  125.0    5.0
 13  5 14  109.5    5.0
#===== 1-3 connections
   4   3   5  12
   6  11   4   6  13  14
   4   1   5  12
   2   2
   6   1   3  12   7   8
   8   2  13  14   9  15  10  16
   5   5   8  10  17
   5   5   7   9  18
   5   6  15   8  18
   5   6  16   7  17
   2   2
   4   1   3   5
   4   2   6  14
   4   2   6  13
   3   6   9
   3   6  10
   3   7  10
   3   8   9
#===== Hydrogen bond donor groups
    2 1.00 0.00 0.20
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0  0.00  0.19  0.13  0.26 -0.58
 0.1 0.0 0.38 -0.48 0.0  0.13 -0.58  0.00  0.26  0.19
 0.1 0.0 0.38 -0.48 0.0 -0.05  0.68  0.00 -0.24  0.61
 0.1 0.0 0.38 -0.48 0.0  0.0   0.19  0.13  0.26 -0.58 
#==8=========================================================================
*
ILE I  113.160
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  19  18   2   2   1  18   0   8  32   1   1   1  8    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG1 CG2 CD1 H   HA  HB 2HG13HG11HG2
2HG23HG21HD12HD13HD1
#===== Bonds
1 2  2 3  3 4  2 5  5 6  5 7  6 8  1 9  2 10   5 11   
6 12   6 13   7 14   7 15   7 16 8 17  8 18 8 19
#===== Rotatable bonds
  2  5    6  7  8 11 12 13 14 15 16 17 18 19
  5  6    8 12 13 17 18 19
#===== Torsion angles
1 2 5 6  2 5 6 8
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
 -0.0070  -1.0190   1.1360   F           
  1.2250  -1.9010   1.0720   F           
 -1.2700  -1.8500   1.0970   F           
  1.5200  -2.6100   2.3730   F           
  2.1000   0.4130   0.0120   F         1 
  0.0050  -0.4900  -0.8700   F           
  0.0120  -0.5200   2.0000   F           
  1.0950  -2.5830   0.3550   F           
  2.0110  -1.3320   0.8340   F           
 -1.3150  -2.3570   0.2360   F           
 -1.2660  -2.4900   1.8650   F           
 -2.0630  -1.2470   1.1740   F           
  0.7600  -3.2180   2.6000   F           
  2.3600  -3.1450   2.2800   F           
  1.6340  -1.9350   3.1010   F           
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.520 0.030
1.513 0.039
1.513 0.039
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  7  110.5  1.7
 2  5  6  110.4  1.7
 7  5  6  110.7  3.0
 5  6  8  113.8  2.1
  1  2 10  109.5    5.0
  2  1  9  120.0    5.0
  2  5 11  109.5    5.0
  3  2 10  109.5    5.0
  5  2 10  109.5    5.0
  5  6 12  109.5    5.0
  5  6 13  109.5    5.0
  5  7 14  109.5    5.0
  5  7 15  109.5    5.0
  5  7 16  109.5    5.0
  6  5 11  109.5    5.0
  6  8 17  109.5    5.0
  6  8 18  109.5    5.0
  6  8 19  109.5    5.0
  7  5 11  109.5    5.0
  8  6 12  109.5    5.0
  8  6 13  109.5    5.0
 12  6 13  109.5    5.0
 14  7 15  109.5    5.0
 14  7 16  109.5    5.0
 15  7 16  109.5    5.0
 17  8 18  109.5    5.0
 17  8 19  109.5    5.0
 18  8 19  109.5    5.0
#===== 1-3 connections
   4   3   5  10
   6   9   4   6   7  11
   4   1   5  10
   2   2
  10   1   3  10   8  12  13  14  15  16
   7   2   7  11  17  18  19
   4   2   6  11
   4   5  12  13
   2   2
   4   1   3   5
   4   2   6   7
   4   5   8  13
   4   5   8  12
   4   5  15  16
   4   5  14  16
   4   5  14  15
   4   6  18  19
   4   6  17  19
   4   6  17  18
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0 0.0
#==9=========================================================================
*
LYS K  128.170
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  22  21   4   4   1  21   0   9  38   1   1   1  8    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD  CE  NZ  H   HA 2HB 3HB 2HG 
3HG 2HD 3HD 2HE 3HE  HZ1 HZ2 HZ3
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  7 8  8 9  1 10   2 11    
5 12   5 13   6 14   6 15   7 16   7 17 8 18 8 19 9 20 9 21 
9 22
#===== Rotatable bonds
  2  5    6  7  8  9 12 13 14 15 16 17 18 19 20 21 22
  5  6    7  8  9 14 15 16 17 18 19 20 21 22
  6  7    8  9 16 17 18 19 20 21 22
  7  8    9 18 19 20 21 22
#===== Torsion angles
1 2 5 6  2 5 6 7  5 6 7 8  6 7 8 9
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
 -0.02981 -0.92231  1.21851  F
  1.15411 -1.87741  1.31451  F
  1.05961 -2.76251  2.54641  F
  2.22441 -3.71901  2.66271  F
  2.09011 -4.56191  3.87831  T+ 5
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.55997 -0.82850  F
 -0.04044 -0.35367  2.04103  F
 -0.86802 -1.46620  1.17866  F
  1.17117 -2.45633  0.49931  F
  1.99860 -1.34395  1.36201  F
  1.04046 -2.18011  3.35909  F
  0.21311 -3.29254  2.49636  F
  2.24756 -4.30832  1.85514  F
  3.07489 -3.19591  2.71789  F
  2.06950 -3.97681  4.68901  F
  1.24209 -5.08931  3.82619  F
  2.86833 -5.18705  3.93808  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.520 0.030
1.520 0.030
1.520 0.030
1.489 0.030
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  114.1  2.0
 5  6  7  111.3  2.3
 6  7  8  111.3  2.3
 7  8  9  111.9  3.2
  1  2 11  109.5    5.0
  2  1 10  120.0    5.0
  2  5 12  109.5    5.0
  2  5 13  109.5    5.0
  3  2 11  109.5    5.0
  5  2 11  109.5    5.0
  5  6 14  109.5    5.0
  5  6 15  109.5    5.0
  6  5 12  109.5    5.0
  6  5 13  109.5    5.0
  6  7 16  109.5    5.0
  6  7 17  109.5    5.0
  7  6 14  109.5    5.0
  7  6 15  109.5    5.0
  7  8 18  109.5    5.0
  7  8 19  109.5    5.0
  8  7 16  109.5    5.0
  8  7 17  109.5    5.0
  8  9 20  109.5    5.0
  8  9 21  109.5    5.0
  8  9 22  109.5    5.0
  9  8 18  109.5    5.0
  9  8 19  109.5    5.0
 12  5 13  109.5    5.0
 14  6 15  109.5    5.0
 16  7 17  109.5    5.0
 18  8 19  109.5    5.0
 20  9 21  109.5    5.0
 20  9 22  109.5    5.0
 21  9 22  109.5    5.0
#===== 1-3 connections
   4   3   5  11
   6  10   4   6  12  13
   4   1   5  11
   2   2
   7   1   3  11   7  14  15
   7   2  12  13   8  16  17
   7   5  14  15   9  18  19
   7   6  16  17  20  21  22
   4   7  18  19
   2   2
   4   1   3   5
   4   2   6  13
   4   2   6  12
   4   5   7  15
   4   5   7  14
   4   6   8  17
   4   6   8  16
   4   7   9  19
   4   7   9  18
   4   8  21  22
   4   8  20  22
   4   8  20  21
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0 0.127 0.873
#=10=========================================================================
*
LEU L  113.160
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  19  18   2   2   1  18   0   8  32   1   1   1  8    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD1 CD2 H   HA 2HB 3HB  HG 1HD1
2HD13HD11HD22HD23HD2
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8  1 9  2 10  5 11 
5 12 6 13 7 14 7 15 7 16 8 17 8 18 8 19
#===== Rotatable bonds
  2  5    6  7  8 11 12 13 14 15 16 17 18 19
  5  6    7  8 13 14 15 16 17 18 19
#===== Torsion angles
1 2 5 6  2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.09921 -0.94171  1.20611  F
  1.29671 -1.87931  1.25351  F
  1.24881 -2.74481  2.50801  F
  1.27451 -2.75501  0.00631  F
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.55133 -0.83427  F
  0.12668 -0.37528  2.02976  F
 -0.72638 -1.50588  1.21715  F
  2.13755 -1.33806  1.25984  F
  1.26321 -2.15867  3.31809  F
  0.41009 -3.28936  2.50543  F
  2.04172 -3.35396  2.52316  F
  0.42666 -3.28486 -0.01389  F
  1.32570 -2.17675 -0.80793  F
  2.05660 -3.37796  0.02247  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.530 0.020
1.521 0.033
1.521 0.033
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  116.3  3.5
 5  6  7  110.7  3.0
 5  6  8  110.7  3.0
 7  6  8  110.8  2.2
  1  2 10  109.5    5.0
  2  1  9  120.0    5.0
  2  5 11  109.5    5.0
  2  5 12  109.5    5.0
  3  2 10  109.5    5.0
  5  2 10  109.5    5.0
  5  6 13  109.5    5.0
  6  5 11  109.5    5.0
  6  5 12  109.5    5.0
  6  7 14  109.5    5.0
  6  7 15  109.5    5.0
  6  7 16  109.5    5.0
  6  8 17  109.5    5.0
  6  8 18  109.5    5.0
  6  8 19  109.5    5.0
  7  6 13  109.5    5.0
  8  6 13  109.5    5.0
 11  5 12  109.5    5.0
 14  7 15  109.5    5.0
 14  7 16  109.5    5.0
 15  7 16  109.5    5.0
 17  8 18  109.5    5.0
 17  8 19  109.5    5.0
 18  8 19  109.5    5.0
#===== 1-3 connections
   4   3   5  10
   6   9   4   6  11  12
   4   1   5  10
   2   2
   7   1   3  10   7   8  13
  10   2  11  12  14  15  16  17  18  19
   4   5   8  13
   4   5   7  13
   2   2
   4   1   3   5
   4   2   6  12
   4   2   6  11
   4   5   7   8
   4   6  15  16
   4   6  14  16
   4   6  14  15
   4   6  18  19
   4   6  17  19
   4   6  17  18
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0 0.0
#=11=========================================================================
* 
MET M  131.199
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  17  16   3   3   1  16   0   8  27   1   1   1  7    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  SD  CE  H   HA 2HB 3HB 2HG 3HG 
1HE 2HE 3HE
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  7 8 1 9 2 10 5 11 
5 12 6 13 6 14 8 15 8 16 8 17
#===== Rotatable bonds
  2  5    6  7  8 11 12 13 14 15 16 17
  5  6    7  8 13 14 15 16 17
  6  7    8 15 16 17
#===== Torsion angles
1 2 5 6  2 5 6 7  5 6 7 8
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.02191 -0.96601  1.18101  F
  1.22441 -1.87211  1.13071  F
  1.13121 -2.81761 -0.43071  T-
 -0.43601 -3.66141 -0.18901  F
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.53910 -0.84222  F
  0.04684 -0.43973  2.03096  F
 -0.80717 -1.52472  1.15931  F
  2.06335 -1.32829  1.15151  F
  1.21370 -2.49801  1.91054  F
 -0.38411 -4.23273  0.63007  F
 -1.16481 -2.98567 -0.07848  F
 -0.63154 -4.23555 -0.98408  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.520 0.030
1.803 0.034
1.791 0.059
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  114.1  2.0
 5  6  7  112.7  3.0
 6  7  8  100.9  2.2
  1  2 10  109.5    5.0
  2  1  9  120.0    5.0
  2  5 11  109.5    5.0
  2  5 12  109.5    5.0
  3  2 10  109.5    5.0
  5  2 10  109.5    5.0
  5  6 13  109.5    5.0
  5  6 14  109.5    5.0
  6  5 11  109.5    5.0
  6  5 12  109.5    5.0
  7  6 13  109.5    5.0
  7  6 14  109.5    5.0
  7  8 15  109.5    5.0
  7  8 16  109.5    5.0
  7  8 17  109.5    5.0
 11  5 12  109.5    5.0
 13  6 14  109.5    5.0
 15  8 16  109.5    5.0
 15  8 17  109.5    5.0
 16  8 17  109.5    5.0
#===== 1-3 connections
   4   3   5  10
   6   9   4   6  11  12
   4   1   5  10
   2   2
   7   1   3  10   7  13  14
   5   2  11  12   8
   7   5  13  14  15  16  17
   2   6
   2   2
   4   1   3   5
   4   2   6  12
   4   2   6  11
   4   5   7  14
   4   5   7  13
   4   7  16  17
   4   7  15  17
   4   7  15  16
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0 0.0
#=12=========================================================================
*
ASN N  114.104
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  14  13   2   2   1  13   0   8  20   1   1   1  8    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  OD1 ND2 H   HA 2HB 3HB 1HD22HD2
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8 1 9 2 10 5 11 
5 12 8 13 8 14
#===== Rotatable bonds
  2  5    6  7  8 11 12 13 14
  5  6    7  8 13 14
#===== Torsion angles
1 2 5 6  2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.04971 -1.02491  1.12921  F     1
  1.23931 -1.96531  1.04351  F     1
  2.06011 -1.86711  0.12571  T-    1
  1.34351 -2.87011  2.01451  T+ 3  1
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.51012 -0.86008  F
  0.09602 -0.53470  1.99958  F
 -0.78775 -1.57053  1.09831  F
  0.65702 -2.90544  2.74079  F     1   1
  2.10827 -3.51442  2.01759  F     1   1
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.516 0.025
1.231 0.020
1.328 0.021
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  112.6  1.0
 5  6  8  116.4  1.5
 5  6  7  120.8  2.0
 8  6  7  122.6  1.0
  1  2 10  109.5    5.0
  2  1  9  120.0    5.0
  2  5 11  109.5    5.0
  2  5 12  109.5    5.0
  3  2 10  109.5    5.0
  5  2 10  109.5    5.0
  6  5 11  109.5    5.0
  6  5 12  109.5    5.0
  6  8 13  120.0    5.0
  6  8 14  120.0    5.0
 11  5 12  109.5    5.0
 13  8 14  120.0    5.0
#===== 1-3 connections
   4   3   5  10
   6   9   4   6  11  12
   4   1   5  10
   2   2
   6   1   3  10   7   8
   6   2  11  12  13  14
   3   5   8
   3   5   7
   2   2
   4   1   3   5
   4   2   6  12
   4   2   6  11
   3   6  14
   3   6  13
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.38 -0.48 0.1
#=13=========================================================================
*
PRO P   97.117
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  14  14   0   0   1  14   0   7  26   1   1   1  8    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD  HA 2HB 3HB 2HG 3HG 2HD 3HD
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  7 1 2 8 5 9 5 10 
6 11 6 12 7 13 7 14
#===== Rotatable bonds
#===== Torsion angles
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  F
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.12371 -0.78261  1.31391  F
  0.89491  0.13851  2.22061  F
  1.87411  0.86171  1.30571  F
  0.00002 -0.61913 -0.78527  F
 -0.77928 -0.98128  1.69487  F
  0.62091 -1.63849  1.17159  F
  0.28046  0.79103  2.66409  F
  1.38550 -0.38403  2.91794  F
  2.02825  1.79934  1.61730  F
  2.74773  0.37698  1.26302  F
#===== Standard bonds
1.466 0.015
1.525 0.021
1.231 0.020
1.530 0.020
1.492 0.050
1.503 0.034
1.473 0.014
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.8  2.5
 2  3  4  120.8  1.7
 1  2  5  103.0  1.1
 3  2  5  110.1  1.9
 2  5  6  104.5  1.9
 5  6  7  106.1  3.2
 6  7  1  103.2  1.5
 7  1  2  112.0  1.4
  1  2  8  109.5    5.0
  1  7 13  109.5    5.0
  1  7 14  109.5    5.0
  2  5  9  109.5    5.0
  2  5 10  109.5    5.0
  3  2  8  109.5    5.0
  5  2  8  109.5    5.0
  5  6 11  109.5    5.0
  5  6 12  109.5    5.0
  6  5  9  109.5    5.0
  6  5 10  109.5    5.0
  6  7 13  109.5    5.0
  6  7 14  109.5    5.0
  7  6 11  109.5    5.0
  7  6 12  109.5    5.0
  9  5 10  109.5    5.0
 11  6 12  109.5    5.0
 13  7 14  109.5    5.0
#===== 1-3 connections
   7   3   5   8   6  13  14
   6   7   4   6   9  10
   4   1   5   8
   2   2
   7   1   3   8   7  11  12
   7   2   9  10   1  13  14
   5   2   5  11  12
   4   1   3   5
   4   2   6  10
   4   2   6   9
   4   5   7  12
   4   5   7  11
   4   1   6  14
   4   1   6  13
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.38 -0.48 0.0 0.0 0.0 0.0
#=14=========================================================================
*
GLN Q  128.131
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  17  16   3   3   1  16   0   9  26   1   1   1  9    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD  OE1 NE2 H   HA 2HB 3HB 2HG 
3HG 1HE22HE2
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  7 8  7 9 1 10  2 11  
5 12  5 13  6 14  6 15  9 16  9 17
#===== Rotatable bonds
  2  5    6  7  8  9 12 13 14 15 16 17
  5  6    7  8  9 14 15 16 17
  6  7    8  9 16 17
#===== Torsion angles
1 2 5 6  2 5 6 7  5 6 7 8
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
 -0.01791 -0.98631  1.17851  F
  1.15911 -1.93941  1.18751  F     1
  1.11651 -2.89591  2.35421  F     1
  0.21441 -2.88171  3.19101  T-    1
  2.14441 -3.75211  2.42701  T+ 3  1
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.53282 -0.84621  F
 -0.00715 -0.46213  2.03005  F
 -0.85981 -1.52363  1.12868  F
  1.15090 -2.46785  0.33858  F
  2.00356 -1.40635  1.23996  F
  1.69020 -4.10292  3.24593  F     1   1
  3.00798 -4.15443  2.12307  F     1   1
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.520 0.030
1.516 0.025
1.231 0.020
1.328 0.021
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  114.1  2.0
 5  6  7  112.6  1.7
 6  7  9  116.4  1.5
 6  7  8  120.8  2.0
 9  7  8  122.6  1.0
  1  2 11  109.5    5.0
  2  1 10  120.0    5.0
  2  5 12  109.5    5.0
  2  5 13  109.5    5.0
  3  2 11  109.5    5.0
  5  2 11  109.5    5.0
  5  6 14  109.5    5.0
  5  6 15  109.5    5.0
  6  5 12  109.5    5.0
  6  5 13  109.5    5.0
  7  6 14  109.5    5.0
  7  6 15  109.5    5.0
  7  9 16  120.0    5.0
  7  9 17  120.0    5.0
 12  5 13  109.5    5.0
 14  6 15  109.5    5.0
 16  9 17  120.0    5.0
#===== 1-3 connections
   4   3   5  11
   6  10   4   6  12  13
   4   1   5  11
   2   2
   7   1   3  11   7  14  15
   6   2  12  13   8   9
   6   5  14  15  16  17
   3   6   9
   3   6   8
   2   2
   4   1   3   5
   4   2   6  13
   4   2   6  12
   4   5   7  15
   4   5   7  14
   3   7  17
   3   7  16
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.38 -0.48 0.1
#=15=========================================================================
*
ARG R  156.188
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  24  23   4   5   1  23   0  11  38   1   1   1  11    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD  NE  CZ  NH1 NH2 H   HA 2HB 
3HB 2HG 3HG 2HD 3HD  HE 1HH12HH11HH22HH2
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  7 8  8 9  9 10  9 11 
1 12  2 13  5 14  5 15  6 16  6 17 7 18  7 19  8 20  10 21  
10 22  11 23  11 24 
#===== Rotatable bonds
  2  5    6  7  8  9 10 11 14 15 16 17 18 19 20 21 22 23 24
  5  6    7  8  9 10 11 16 17 18 19 20 21 22 23 24
  6  7    8  9 10 11 18 19 20 21 22 23 24
  7  8    9 10 11 20 21 22 23 24
#===== Torsion angles
1 2 5 6  2 5 6 7  5 6 7 8  6 7 8 9  7 8 9 10
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
 -0.01911 -0.93181  1.22161  F
  1.15841 -1.89161  1.32361  F
  1.06001 -2.75511  2.53351  F     1
  2.18381 -3.67671  2.63811  T+ 2  1
  2.32141 -4.54891  3.63131  F     1
  1.54051 -4.50711  4.70541  T+ 3  1
  3.29551 -5.45391  3.57311  T+ 3  1
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.55794 -0.82986  F
 -0.02699 -0.36331  2.04425  F
 -0.85818 -1.47450  1.18364  F
  1.17421 -2.47328  0.51034  F
  2.00538 -1.36209  1.37098  F
  1.04350 -2.17225  3.34592  F
  0.21235 -3.28345  2.48526  F
  2.88654 -3.65163  1.92710  F     1   1
  0.81678 -3.82015  4.77098  F     1   1
  1.67649 -5.16261  5.44827  F     1   1
  3.90266 -5.48465  2.77911  F     1   1
  3.42028 -6.10317  4.32338  F     1   1
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.520 0.030
1.520 0.030
1.460 0.018
1.329 0.014
1.326 0.018
1.326 0.018
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  114.1  2.0
 5  6  7  114.1  2.0
 6  7  8  110.5  1.7
 7  8  9  121.7  1.8
 8  9 10  120.0  1.9
 8  9 11  120.0  1.9
10  9 11  119.7  1.8
  1  2 13  109.5    5.0
  2  1 12  120.0    5.0
  2  5 14  109.5    5.0
  2  5 15  109.5    5.0
  3  2 13  109.5    5.0
  5  2 13  109.5    5.0
  5  6 16  109.5    5.0
  5  6 17  109.5    5.0
  6  5 14  109.5    5.0
  6  5 15  109.5    5.0
  6  7 18  109.5    5.0
  6  7 19  109.5    5.0
  7  6 16  109.5    5.0
  7  6 17  109.5    5.0
  7  8 20  120.0    5.0
  8  7 18  109.5    5.0
  8  7 19  109.5    5.0
  9  8 20  120.0    5.0
  9 10 21  120.0    5.0
  9 10 22  120.0    5.0
  9 11 23  120.0    5.0
  9 11 24  120.0    5.0
 14  5 15  109.0    5.0
 16  6 17  109.5    5.0
 18  7 19  109.5    5.0
 21 10 22  120.0    5.0
 23 11 24  120.0    5.0
#===== 1-3 connections
   4   3   5  13
   6  12   4   6  14  15
   4   1   5  13
   2   2
   7   1   3  13   7  16  17
   7   2  14  15   8  18  19
   6   5  16  17   9  20
   6   6  18  19  10  11
   7   7  20  21  22  23  24
   3   8  11
   3   8  10
   2   2
   4   1   3   5
   4   2   6  15
   4   2   6  14
   4   5   7  17
   4   5   7  16
   4   6   8  19
   4   6   8  18
   3   7   9
   3   9  22
   3   9  21
   3   9  24
   3   9  23
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.09 0.13 0.34 0.22 0.22
#=16=========================================================================
*
SER S   87.078
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  11  10   1   1   1  10   0   6  15   1   1   1  5    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  OG  H   HA 2HB 3HB  HG
#===== Bonds
1 2  2 3  3 4  2 5  5 6 1 7 2 8 5 9 5 10 6 11
#===== Rotatable bonds
  2  5    6  9 10 11
#===== Torsion angles
1 2 5 6
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.01411 -0.97081  1.17401  F
  0.05551 -0.29631  2.41811  T  4
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.53613 -0.84412  F
 -0.81229 -1.53277  1.13850  F
  0.81989 -1.55812  1.09793  F
  0.06696  0.69130  2.26160  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.417 0.020
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  111.1  2.0
  1  2  8  109.5    5.0
  2  1  7  120.0    5.0
  2  5  9  109.5    5.0
  2  5 10  109.5    5.0
  3  2  8  109.5    5.0
  5  2  8  109.5    5.0
  5  6 11  109.5    5.0
  6  5  9  109.5    5.0
  6  5 10  109.5    5.0
  9  5 10  109.5    5.0
#===== 1-3 connections
   4   3   5   8
   6   7   4   6   9  10
   4   1   5   8
   2   2
   5   1   3   8  11
   4   2   9  10
   2   2
   4   1   3   5
   4   2   6  10
   4   2   6   9
   2   5
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.15 -0.15
#=17=========================================================================
*
THR T  101.105
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  14  13   1   1   1  13   0   7  21   1   1   1  7    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  OG1 CG2 H   HA  HB  HG11HG22HG23HG2
#===== Bonds
1 2  2 3  3 4  2 5  5 6  5 7 1 8 2 9 5 10 6 11 
7 12 7 13 7 14
#===== Rotatable bonds
  2  5    6  7 10 11 12 13 14
#===== Torsion angles
1 2 5 6
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.00141 -0.98081  1.24201  F
  1.22641 -1.77501  1.13621  T  4
 -1.25231 -1.86601  1.32361  F
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.54867 -0.83602  F
  0.05360 -0.43183  2.07621  F
  1.73487 -1.49651  0.32141  F
 -2.06445 -1.28748  1.39931  F
 -1.31845 -2.42710  0.49850  F
 -1.18973 -2.45937  2.12611  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.433 0.016
1.521 0.033
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  7  110.5  1.7
 2  5  6  109.6  1.5
 7  5  6  109.3  2.0
  1  2  9  109.5    5.0
  2  1  8  120.0    5.0
  2  5 10  109.5    5.0
  3  2  9  109.5    5.0
  5  2  9  109.5    5.0
  5  6 11  109.5    5.0
  5  7 12  109.5    5.0
  5  7 13  109.5    5.0
  5  7 14  109.5    5.0
  6  5 10  109.5    5.0
  7  5 10  109.5    5.0
 12  7 13  109.5    5.0
 12  7 14  109.5    5.0
 13  7 14  109.5    5.0
#===== 1-3 connections
   4   3   5   9
   6   8   4   6   7  10
   4   1   5   9
   2   2
   8   1   3   9  11  12  13  14
   4   2   7  10
   4   2   6  10
   2   2
   4   1   3   5
   4   2   6   7
   2   5
   4   5  13  14
   4   5  12  14
   4   5  12  13
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.15 -0.15 0.0
#=18=========================================================================
*
VAL V   99.133
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  16  15   1   1   1  15   0   7  26   1   1   1  7    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG1 CG2 H   HA  HB 1HG12HG13HG11HG2
2HG23HG2
#===== Bonds
1 2  2 3  3 4  2 5  5 6  5 7 1 8 2 9 5 10 6 11 
6 12 6 13 7 14 7 15 7 16
#===== Rotatable bonds
  2  5    6  7 10 11 12 13 14 15 16
#===== Torsion angles
1 2 5 6
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.06511 -0.98861  1.20431  F
 -1.15211 -1.90321  1.27151  F
  1.34881 -1.79721  1.19551  F
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.53798 -0.84294  F
  0.06796 -0.43811  2.03914  F
 -1.97945 -1.35027  1.37031  F
 -1.21022 -2.44298  0.43171  F
 -1.06530 -2.51647  2.05661  F
  1.40311 -2.32718  0.34923  F
  2.13264 -1.17885  1.25247  F
  1.35674 -2.41854  1.97902  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.521 0.033
1.521 0.033
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  110.5  1.7
 2  5  7  110.5  1.7
 6  5  7  110.8  2.2
  1  2  9  109.5    5.0
  2  1  8  120.0    5.0
  2  5 10  109.5    5.0
  3  2  9  109.5    5.0
  5  2  9  109.5    5.0
  5  6 11  109.5    5.0
  5  6 12  109.5    5.0
  5  6 13  109.5    5.0
  5  7 14  109.5    5.0
  5  7 15  109.5    5.0
  5  7 16  109.5    5.0
  6  5 10  109.5    5.0
  7  5 10  109.5    5.0
 11  6 12  109.5    5.0
 11  6 13  109.5    5.0
 12  6 13  109.5    5.0
 14  7 15  109.5    5.0
 14  7 16  109.5    5.0
 15  7 16  109.5    5.0
#===== 1-3 connections
   4   3   5   9
   6   8   4   6   7  10
   4   1   5   9
   2   2
  10   1   3   9  11  12  13  14  15  16
   4   2   7  10
   4   2   6  10
   2   2
   4   1   3   5
   4   2   6   7
   4   5  12  13
   4   5  11  13
   4   5  11  12
   4   5  15  16
   4   5  14  16
   4   5  14  15
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0
#=19=========================================================================
*
TRP W  186.213
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  24  25   2   2   1  25   0  14  41   1   1   1 20    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 
 H   HA  HB2 HB3 HD1 HE1 HE3 HZ2 HZ3 HH2 
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8  7 9  9 10  8 10  
8 11  10 12  11 13  12 14 13 14 1 15 2 16 5 17 5 18 7 19 
9 20 11 21 12 22 13 23 14 24
#===== Rotatable bonds
  2  5    6  7  8  9 10 11 12 13 14 17 18 19 20 21 22 23 24
  5  6    7  8  9 10 11 12 13 14 19 20 21 22 23 24
#===== Torsion angles
1 2 5 6  2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.04591 -1.01101  1.18031  F     1
  1.19811 -1.94811  1.17901  F     1
  2.42381 -1.76731  1.76271  F     1
  1.24141 -3.25121  0.57561  F     1
  3.21721 -2.87001  1.56001  T+ 2  1
  2.52081 -3.78721  0.82521  F     1
  0.31981 -3.99721 -0.15771  F     1
  2.89771 -5.03591  0.38271  F     1
  0.68801 -5.23691 -0.60511  F     1
  1.95811 -5.75211 -0.34400  F     1
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.52634 -0.85026  F
  0.07567 -0.48814  2.03220  F
 -0.79285 -1.55479  1.15230  F
  2.70169 -0.94860  2.26520  F     1   1
  4.15239 -2.98303  1.89567  F     1   1
 -0.59126 -3.63397 -0.35269  F     1   1
  3.80539 -5.40835  0.57613  F     1   1
  0.03474 -5.78367 -1.12885  F     1   1
  2.19680 -6.65979 -0.68918  F     1   1
#===== Standard bonds
1.458 0.019  
1.525 0.021
1.231 0.020 
1.530 0.020 
1.498 0.031 
1.365 0.025 
1.433 0.018  
1.374 0.021 
1.370 0.011  
1.409 0.017 
1.398 0.016  
1.394 0.021 
1.382 0.030 
1.368 0.019  
1.400 0.025  
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  113.6  1.9
 5  6  7  126.9  1.5
 5  6  8  126.8  1.4
 7  6  8  106.3  1.6
 6  7  9  110.2  1.3
 7  9 10  108.9  1.8
 9 10 12  130.1  1.5
 9 10  8  107.4  1.3
12 10  8  122.4  1.0
10 12 14  117.5  1.3
12 14 13  121.5  1.3
14 13 11  121.1  1.3
13 11  8  118.6  1.3
11  8  6  133.9  1.0
11  8 10  118.8  1.0
 6  8 10  107.2  1.2
  1  2 16  109.5    5.0
  2  1 15  120.0    5.0
  2  5 17  109.5    5.0
  2  5 18  109.5    5.0
  3  2 16  109.5    5.0
  5  2 16  109.5    5.0
  6  5 17  109.5    5.0
  6  5 18  109.5    5.0
  6  7 19  125.0    5.0
  7  9 20  125.0    5.0
  8 11 21  120.0    5.0
  9  7 19  125.0    5.0
 10  9 20  125.0    5.0
 10 12 22  120.0    5.0
 11 13 23  120.0    5.0
 12 14 24  120.0    5.0
 13 11 21  120.0    5.0
 13 14 24  120.0    5.0
 14 12 22  120.0    5.0
 14 13 23  120.0    5.0
 17  5 18  109.5    5.0
#===== 1-3 connections
   4   3   5  16
   6  15   4   6  17  18
   4   1   5  16
   2   2
   6   1   3  16   7   8
   8   2  17  18   9  19  10  11
   5   5   8  10  20
   7   5   7   9  12  13  21
   5   6  19   8  12
   7   6  11   7  20  14  22
   5   6  10  14  23
   5   8   9  13  24
   5   8  21  12  24
   5  10  22  11  23
   2   2
   4   1   3   5
   4   2   6  18
   4   2   6  17
   3   6   9
   3   7  10
   3   8  13
   3  10  14
   3  11  14
   3  12  13
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 -0.14 0.0 0.0 0.14 0.0 0.0 0.0 0.0 0.0
#=20=========================================================================
*
TYR Y  163.176
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  21  21   2   2   1  21   0  12  33   1   1   1  15    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OH  H   HA  
 HB2 HB3 HD1 HD2 HE1 HE2 HH
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8  7 9  8 10   9 11  
10 11  11 12 1 13 2 14 5 15 5 16 7 17 8 18 9 19 10 20 
12 21
#===== Rotatable bonds
  2  5    6  7  8  9 10 11 12 15 16 17 18 19 20 21
  5  6    7  8  9 10 11 12 17 18 19 20 21
#===== Torsion angles
1 2 5 6  2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.01201 -0.96311  1.19561  F     1
  1.19321 -1.91361  1.23451  F     1
  2.38131 -1.56971  1.89081  F     1
  1.13581 -3.15221  0.59721  F     1
  3.48011 -2.42711  1.88941  F     1
  2.22791 -4.01221  0.58291  F     1
  3.40351 -3.64371  1.21931  F     1
  4.46881 -4.49421  1.14111  T  4  1
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.54328 -0.83953  F
  0.02334 -0.41820  2.03403  F
 -0.82548 -1.50865  1.16414  F
  2.44252 -0.69355  2.36892  F     1   1
  0.29023 -3.42815  0.14020  F     1   1
  4.31946 -2.16896  2.36779  F     1   1
  2.16592 -4.89388  0.11512  F     1   1
  4.21480 -5.29802  0.60318  F         1
#===== Standard bonds
1.458 0.019  
1.525 0.021
1.231 0.020 
1.530 0.020 
1.512 0.022
1.389 0.021
1.389 0.021
1.382 0.030
1.382 0.030
1.378 0.024
1.378 0.024
1.376 0.021
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  113.6  1.9
 5  6  7  120.8  1.5
 5  6  8  120.8  1.5
 6  7  9  121.2  1.5
 6  8 10  121.2  1.5
 7  9 11  119.6  1.8
 8 10 11  119.6  1.8
 9 11 10  120.3  2.0
 7  6  8  118.1  1.5
 9 11 12  119.9  3.0
10 11 12  119.9  3.0
  1  2 14  109.5    5.0
  2  1 13  120.0    5.0
  2  5 15  109.5    5.0
  2  5 16  109.5    5.0
  3  2 14  109.5    5.0
  5  2 14  109.5    5.0
  6  5 15  109.5    5.0
  6  5 16  109.5    5.0
  6  7 17  120.0    5.0
  6  8 18  120.0    5.0
  7  9 19  120.0    5.0
  8 10 20  120.0    5.0
  9  7 17  120.0    5.0
 10  8 18  120.0    5.0
 11  9 19  120.0    5.0
 11 10 20  120.0    5.0
 11 12 21  120.0    5.0
 15  5 16  109.5    5.0
#===== 1-3 connections
   4   3   5  14
   6  13   4   6  15  16
   4   1   5  14
   2   2
   6   1   3  14   7   8
   8   2  15  16   9  17  10  18
   5   5   8  11  19
   5   5   7  11  20
   5   6  17  10  12
   5   6  18   9  12
   6   7  19   8  20  21
   3   9  10
   2   2
   4   1   3   5
   4   2   6  16
   4   2   6  15
   3   6   9
   3   6  10
   3   7  11
   3   8  11
   2  11
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0 0.0 0.0 0.0 0.15 -0.15
#=21=========================================================================
*
HYP O  113.117
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  15  15   0   0   1   8   0   8   0   1   1   1  0    1
#===== Atom names
#1   2   3   4   5   6   7   8   9  10  11  12  13  14
 N   CA  C   O   CB  CG  CD  OD  HA 2HB 3HB 1HD 2HD  HG
 HOD
#===== Bonds
1 2   2 3   3 4   2 5   5 6   6 7   7 1   8 6   2 9   5 10  
5 11  6 14  7 12  7 13  8 15
#===== Rotatable bonds
#===== Torsion angles
#===== Coordinates, hydrogen bonding and planarity
  15.963   9.809   -4.565    F
  14.885   9.685   -3.555    F
  13.586   9.162   -4.229    F
  13.534   9.088   -5.442    T-
  14.693  11.075   -2.951    F
  15.519  12.046   -3.857    F
  16.249  11.218   -4.936    F
  16.547  12.547   -3.014    T  4
  15.239   9.006   -2.800  
  15.051  11.066   -1.930  
  13.644  11.335   -2.930 
  17.311  11.401   -4.924  
  15.851  11.417   -5.920  
  16.824  13.398   -3.359  
  14.949  12.858   -4.277  
#===== Standard bonds
1.453 0.020
1.530 0.020
1.230 0.020
1.530 0.020
1.530 0.020
1.530 0.020
1.453 0.020
1.400 0.020
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.8  2.5
 2  3  4  120.8  1.7
 1  2  5  103.0  1.1
 3  2  5  110.1  1.9
 2  5  6  104.5  1.9
 5  6  7  106.1  3.2
 6  7  1  103.2  1.5
 7  1  2  112.0  1.4
 5  6  8  109.5  5.0
 7  6  8  109.5  5.0
  1  2  9  109.5    5.0
  1  7 12  109.5    5.0
  1  7 13  109.5    5.0
  2  5 10  109.5    5.0
  2  5 11  109.5    5.0
  3  2  9  109.5    5.0
  5  2  9  109.5    5.0
  5  6 14  109.5    5.0
  6  5 10  109.5    5.0
  6  5 11  109.5    5.0
  6  7 12  109.5    5.0
  6  7 13  109.5    5.0
  7  6 14  109.5    5.0
 10  5 11  109.5    5.0
 14  6  8  109.5    5.0
 13  7 12  109.5    5.0
  6  8 15  109.5    5.0
#===== 1-3 connections
   4   3   5   6
   4   4   6   7
   3   1   5
   2   2
   5   1   3   7   8
   3   1   2
   4   2   5   8
   3   5   7
   4   1   3   5
   4   2   6  11
   4   2   6  10
   4   1   6  13
   4   1   6  12
   4   5   7   8
   2   6
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.38 -0.48 0.0 0.0 0.1 0.0 -0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#=22=========================================================================
*
AIB B   85.106
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
   6   5   0   0   1   5   0   6   0   1   1   1  0    1
#===== Atom names
 N   CA  C   O   CB1 CB2
#===== Bonds
1 2  2 3  3 4  2 5  2 6
#===== Rotatable bonds
#===== Torsion angles
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.02021 -0.92681  1.20941  F
  0.02021 -0.92681  1.20941  F
#===== Standard bonds
1.453 0.020
1.530 0.020
1.230 0.020
1.530 0.020
1.530 0.020
#===== Standard angles
#===== 1-3 connections
   4   3   5   6
   2   4
   4   1   5   6
   2   2
   4   1   3   6
   4   1   3   5
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0
#=23=========================================================================
*
PCA P  115.132
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  13  13   0   0   1  13   0   8  23   1   1   1  8    1
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  CD  OE  HA 2HB 3HB 2HG 3HG 
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  7 1  7 8  2 9  5 10 
5 11  6 12  6 13  
#===== Rotatable bonds
#===== Torsion angles
#===== Coordinates, hydrogen bonding and planarity
  1.68600 -6.63900  2.89000  F     1
  2.21200 -5.56500  2.07800  F     1
  1.33800 -5.27800  0.85500  F
  0.16100 -5.64300  0.83400  T-
  3.57200 -6.12700  1.63900  F
  4.00000 -7.04700  2.75800  F     1
  2.77200 -7.29000  3.62800  F     1
  2.69700 -8.25100  4.40900  T  7  1
  2.33800 -4.73500  2.62700
  4.23700 -5.38600  1.51900
  3.47000 -6.64700  0.78800
  4.33800 -7.91300  2.38400
  4.71500 -6.60500  3.30100
#===== Standard bonds
1.466 0.050
1.525 0.050
1.231 0.050
1.530 0.050
1.492 0.050
1.503 0.050
1.333 0.050
1.231 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.8  2.5
 2  3  4  120.8  1.7
 1  2  5  103.0  1.1
 3  2  5  110.1  1.9
 2  5  6  104.5  1.9
 5  6  7  106.1  3.2
 6  7  8  120.8  1.7
 1  7  8  122.6  5.0
 6  7  1  103.2  1.5
 7  1  2  112.0  1.4
  1  2  9  109.5    5.0
  2  5 10  109.5    5.0
  2  5 11  109.5    5.0
  3  2  9  109.5    5.0
  5  2  9  109.5    5.0
  5  6 12  109.5    5.0
  5  6 13  109.5    5.0
  6  5 10  109.5    5.0
  6  5 11  109.5    5.0
  7  6 12  109.5    5.0
  7  6 13  109.5    5.0
 11  5 10  109.5    5.0
 13  6 12  109.5    5.0
#===== 1-3 connections
   6   6   8   3   5   9
   6   6  11  10   4   7
   4   1   5   9
   2   2
   7   1   3   9   7  12  13
   6   2  10  11   1   8  
   5   2   5  12  13
   3   1   6
   4   1   3   5
   4   2   6  11
   4   2   6  10
   4   5   7  13
   4   5   7  12
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.38 -0.48 0.0 0.0 0.0 0.0
#=24=========================================================================
*
UNK U
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
   6   5   1   1   1   5   0   6   5   1   1   1   0   1
#===== Atom names
 N   CA  C   O   CB  CG1    
#===== Bonds
1 2  2 3  3 4  2 5  5 6  
#===== Rotatable bonds
2 5   6 
#===== Torsion angles
1 2 5 6
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.00001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.06511 -0.98861  1.20431  F
 -1.15211 -1.90321  1.27151  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.521 0.033
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  110.5  1.7
#===== 1-3 connections
3   3   5
3   4   6
3   1   5
2   2
3   1   3
2   2
#===== Hydrogen bond donor groups
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.0 0.0
#============================================================================
*
DINOI
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  21  23   4   4   1  23   0  21  33   1   0   1  33    3
#===== Atom names
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N9  C8  N7  
 C5  C6  O6  N1  C2  N3  C4
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 21  13 14  14 15  15 21  20 21  19 20  18 19
15 16 16 18  16 17 
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 21
 10 12   13 14 15 16 17 18 19 20 21
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
-29.727  -28.046    7.356    F
-26.462  -26.954    9.625    F
-25.812  -26.086   10.079    T-
-30.010  -29.481    7.134    T-
-29.011  -27.356    6.258    T-
-28.928  -27.851    8.724    T-
-29.014  -26.588    9.410    F
-27.742  -26.350   10.208    F
-27.879  -26.983   11.515    T-
-27.092  -28.163   11.569    F     1
-26.472  -28.372   10.187    F
-27.971  -29.281   12.011    F     1
-28.806  -30.080   11.263    F     1
-29.448  -30.984   11.958    T-    1
-29.007  -30.769   13.265    F     1
-29.349  -31.440   14.466    F     1
-30.126  -32.380   14.629    T-    1
-28.666  -30.897   15.566    T+ 2  1
-27.771  -29.844   15.510    F     1
-27.453  -29.218   14.382    T-    1
-28.107  -29.733   13.309    F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.431 0.013
1.593 0.010 
1.440 0.016 
1.511 0.008 
1.528 0.010 
1.446 0.011 
1.420 0.013 
1.521 0.014 
1.518 0.010 
1.459 0.009  
1.374 0.007 
1.375 0.008 
1.305 0.006 
1.388 0.006 
1.379 0.007  
1.350 0.007  
1.323 0.008 
1.373 0.008 
1.419 0.010 
1.391 0.007  
1.237 0.009 
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  114.7  1.5
 2  8  9  105.6  1.0
 2 11 10  102.7  1.4
 3  2  8  110.3  2.2
 3  2 11  110.6  2.7
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.4  1.6
 8  9 10  109.7  1.4
 8  2 11  103.2  1.0
 9 10 11  106.1  1.1
 9 10 12  107.8  0.8
11 10 12  114.2  1.6
10 12 13  127.0  1.3
10 12 21  126.5  1.3
12 13 14  111.1  0.5
12 21 15  105.4  0.4
12 21 20  126.0  0.6
13 14 15  104.3  0.5
13 12 21  106.4  0.4
14 15 16  130.4  0.6
14 15 21  110.8  0.4
15 16 17  128.6  0.6
15 16 18  111.5  0.5
15 21 20  128.6  0.5
16 18 19  125.1  0.6
16 15 21  118.8  0.6
17 16 18  119.9  0.6
18 19 20  123.9  0.6
19 20 21  111.9  0.5
#===== 1-3 connections
2   7
4   7   9  10
3   8  11
3   5   6
3   4   6
4   4   5   8
4   1   2   9
5   3   6  10  11
5   2   7  11  12
5   2   8  13  21
5   3   8   9  12
6   9  11  14  15  20
4  10  15  21
4  12  16  21
6  12  13  17  18  20
4  14  19  21
3  15  18
4  15  17  20
3  16  21
4  12  15  18 
6  10  13  14  16  19
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
DADEA  313.21
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  32  34   4   4   1  34   0  21  59   1   0   1  33    3
#===== Atom names
#---++++----++++----++++----++++----++++----++++----++++----
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N9  C8  N7
 C5  C6  N6  N1  C2  N3  C4  H5''H5' H4' H3' H2''H2' H1'
 H8  H61 H62 H2 
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 21  13 14  14 15  15 21  20 21  19 20  18 19
15 16 16 18  16 17   7 22   7 23   8 24   2 25  11 27  11 26  10 28
13 29 17 30  17 31  19 32
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 21 24 25 26 27 28 29 30 31 32
 10 12   13 14 15 16 17 18 19 20 21 29 30 31 32
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
 -4.997   8.861 -18.658      F
 -5.770   8.950 -14.203      F
 -6.401   9.305 -12.982      T-
 -4.952  10.134 -19.417      T-
 -3.886   7.902 -18.859      T-
 -5.054   9.191 -17.096      T-
 -6.100   9.967 -16.538      F
 -6.741   9.202 -15.370      F
 -7.203   7.945 -15.850      T-
 -6.691   6.892 -15.050      F     1
 -5.448   7.469 -14.369      F
 -6.446   5.690 -15.888      F     1
 -5.269   5.281 -16.459      F     1
 -5.383   4.298 -17.307      T-    1
 -6.737   3.992 -17.263      F     1
 -7.519   3.032 -17.932      F     1
 -7.002   2.158 -18.798      T+ 3  1
 -8.834   2.995 -17.678      T-    1
 -9.342   3.842 -16.787      F     1
 -8.722   4.778 -16.066      T-    1
 -7.395   4.812 -16.377      F     1
 -6.874  10.172 -17.278      F
 -5.701  10.919 -16.189      F
 -7.595   9.783 -15.016      F
 -4.861   9.534 -14.310      F
 -4.589   7.438 -15.029      F
 -5.203   6.947 -13.445      F
 -7.415   6.661 -14.272      F
 -4.326   5.765 -16.269      F     1
 -5.995   2.146 -18.911      F     1
 -7.543   1.397 -19.189      F     1
-10.401   3.759 -16.613      F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.431 0.013
1.593 0.010 
1.440 0.016  
1.511 0.008 
1.528 0.010 
1.446 0.011 
1.420 0.013 
1.521 0.014 
1.518 0.010 
1.462 0.010 
1.373 0.008 
1.374 0.006 
1.311 0.007 
1.388 0.006 
1.383 0.007 
1.344 0.006 
1.331 0.009 
1.339 0.009 
1.406 0.009 
1.351 0.007 
1.335 0.008 
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  114.7  1.5
 2  8  9  105.6  1.0
 2 11 10  102.7  1.4
 3  2  8  110.3  2.2
 3  2 11  110.6  2.7
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.4  1.6
 8  9 10  109.7  1.4
 8  2 11  103.2  1.0
 9 10 11  106.1  1.1
 9 10 12  107.8  0.8
11 10 12  114.2  1.6
10 12 13  127.7  1.8
10 12 21  126.3  1.8
12 13 14  113.8  0.5
12 21 15  105.8  0.4
12 21 20  127.4  0.8
13 14 15  103.9  0.5
13 12 21  105.8  0.4
14 15 16  132.3  0.7
14 15 21  110.7  0.5
15 16 17  123.7  0.8
15 16 18  117.7  0.5
15 21 20  126.8  0.7
16 18 19  118.6  0.6
16 15 21  117.0  0.5
17 16 18  118.6  0.6
18 19 20  129.3  0.5
19 20 21  110.6  0.5
  2  8 24  109.5    5.0
  2 11 26  109.5    5.0
  2 11 27  109.5    5.0
  3  2 25  109.5    5.0
  6  7 22  109.5    5.0
  6  7 23  109.5    5.0
  7  8 24  109.5    5.0
  8  2 25  109.5    5.0
  8  7 22  109.5    5.0
  8  7 23  109.5    5.0
  9  8 24  109.5    5.0
  9 10 28  109.5    5.0
 10 11 26  109.5    5.0
 10 11 27  109.5    5.0
 11  2 25  109.5    5.0
 11 10 28  109.5    5.0
 12 10 28  109.5    5.0
 12 13 29  120.0    5.0
 14 13 29  120.0    5.0
 16 17 30  120.0    5.0
 16 17 31  120.0    5.0
 18 19 32  120.0    5.0
 20 19 32  120.0    5.0
 22  7 23  109.5    5.0
 26 11 27  109.5    5.0
 30 17 31  120.0    5.0
#===== 1-3 connections
2   7
7   7   9  10  24  26  27
4   8  11  25
3   5   6
3   4   6
6   4   5   8  22  23
5   1   2   9  24
8   3   6  10  11  22  23  25
7   2   7  11  12  24  28
7   2   8  13  21  26  27
7   3   8   9  12  25  28
8   9  11  14  15  20  28  29
4  10  15  21
5  12  16  21  29
6  12  13  17  18  20
6  14  19  21  30  31
3  15  18
5  15  17  20  32
3  16  21
5  12  15  18  32
6  10  13  14  16  19
4   6   8  23
4   6   8  22
4   2   7   9
4   3   8  11
4   2  10  27
4   2  10  26
4   9  11  12
3  12  14
3  16  31
3  16  30
3  18  20
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
DTHYT  304.190
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  32  33   4   4   1  33   0  20  59   1   0   1  30    3
#===== Atom names
#---++++----++++----++++----++++----++++----++++----++++----
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N1  C6  C5
 C4  O4  N3  C2  O2  C7  H5''H5' H4' H3' H2''H2' H1' H3 
 H71 H72 H73 H6 
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 18  13 14  18 19  17 18  15 17  14 15  14 20
15 16  7 21   7 22   8 23   2 24  11 25  11 26  10 27  17 28  20 29
20 30 20 31  13 32
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 23 24 25 26 27 28 29 30 31 32
 10 12   13 14 15 16 17 18 19 20 28 29 30 31 32
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
 10.615   -8.801   -0.901    F
 11.767   -4.239   -2.288    F
 12.262   -3.632   -3.472    T-
 11.626   -9.642   -1.564    T-
 10.361   -9.023    0.547    T-
 11.017   -7.267   -1.154    T-
 11.283   -6.811   -2.478    F
 10.814   -5.366   -2.722    F
  9.573   -5.125   -2.094    T-
  9.439   -3.727   -1.936    F     1
 10.837   -3.333   -1.475    F
  8.395   -3.421   -0.919    F     1
  8.180   -4.267    0.139    F     1
  7.137   -4.109    0.982    F     1
  6.195   -3.032    0.774    F     1
  5.153   -2.888    1.403    T-    1
  6.542   -2.150   -0.222    T+ 2  1
  7.637   -2.262   -1.055    F     1
  7.921   -1.363   -1.837    T-    1
  6.987   -5.046    2.168    F     1
 10.731   -7.430   -3.189    F
 12.347   -6.906   -2.698    F
 10.642   -5.260   -3.795    F
 12.570   -4.616   -1.656    F
 10.961   -3.570   -0.421    F
 11.011   -2.276   -1.648    F
  9.218   -3.277   -2.906    F
  5.961   -1.325   -0.357    F     1
  6.117   -4.777    2.768    F
  6.871   -6.073    1.827    F
  7.878   -4.979    2.795    F
  8.855   -5.083    0.303    F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.431 0.013
1.593 0.010 
1.440 0.016 
1.511 0.008 
1.528 0.010 
1.446 0.011 
1.420 0.013 
1.521 0.014 
1.518 0.010 
1.473 0.014
1.378 0.007 
1.376 0.008 
1.339 0.007  
1.220 0.008 
1.373 0.008 
1.382 0.008 
1.445 0.009 
1.496 0.006
1.228 0.009
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  114.7  1.5
 2  8  9  105.6  1.0
 2 11 10  102.7  1.4
 3  2  8  110.3  2.2
 3  2 11  110.6  2.7
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.4  1.6
 8  9 10  109.7  1.4
 8  2 11  103.2  1.0
 9 10 11  106.1  1.1
 9 10 12  107.8  0.8
11 10 12  114.2  1.6
10 12 13  120.4  1.5
10 12 18  118.2  1.6
12 13 14  123.7  0.6
12 18 17  114.6  0.6
12 18 19  123.1  0.8
13 14 15  118.0  0.6
13 12 18  121.3  0.5
13 14 20  122.9  0.6
14 15 16  124.9  0.7
14 15 17  115.2  0.6
15 17 18  127.2  0.6
15 14 20  119.0  0.6
16 15 17  119.9  0.6
17 18 19  122.3  0.6
  2  8 23  109.5    5.0
  2 11 25  109.5    5.0
  2 11 26  109.5    5.0
  3  2 24  109.5    5.0
  6  7 21  109.5    5.0
  6  7 22  109.5    5.0
  7  8 23  109.5    5.0
  8  2 24  109.5    5.0
  8  7 21  109.5    5.0
  8  7 22  109.5    5.0
  9  8 23  109.5    5.0
  9 10 27  109.5    5.0
 10 11 25  109.5    5.0
 10 11 26  109.5    5.0
 11  2 24  109.5    5.0
 11 10 27  109.5    5.0
 12 10 27  109.5    5.0
 12 13 32  120.0    5.0
 14 13 32  120.0    5.0
 14 20 29  109.5    5.0
 14 20 30  109.5    5.0
 14 20 31  109.5    5.0
 15 17 28  120.0    5.0
 18 17 28  120.0    5.0
 21  7 22  109.5    5.0
 25 11 26  109.5    5.0
 29 20 30  109.5    5.0
 29 20 31  109.5    5.0
 30 20 31  109.5    5.0
#===== 1-3 connections
2   7
7   7   9  10  23  25  26
4   8  11  24
3   5   6
3   4   6
6   4   5   8  21  22
5   1   2   9  23
8   3   6  10  11  21  22  24
7   2   7  11  12  23  27
7   2   8  13  18  25  26
7   3   8   9  12  24  27
8   9  11  14  17  19  27  32
5  10  15  18  20
8  12  16  17  29  30  31  32
5  13  18  20  28
3  14  17
5  12  14  16  19
5  10  13  15  28
3  12  17
3  13  15
4   6   8  22
4   6   8  21
4   2   7   9
4   3   8  11
4   2  10  26
4   2  10  25
4   9  11  12
3  15  18
4  14  30  31
4  14  29  31
4  14  29  30
3  12  14
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
DGUAG  329.210
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  33  35   4   4   1  35   0  22  61   1   0   1  35    3
#===== Atom names
#---++++----++++----++++----++++----++++----++++----++++----
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N9  C8  N7  
 C5  C6  O6  N1  C2  N2  N3  C4  H5''H5' H4' H3' H2''H2'
 H1' H8  H1  H21 H22 
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 22  13 14  14 15  15 22  21 22  19 21  18 19
15 16 16 18  16 17  19 20   7 23   7 24   8 25   2 26  11 27  11 28
10 29 13 30  18 31  20 32  20 33
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 21 22 25 26 27 28 29 30 31 32 33
 10 12   13 14 15 16 17 18 19 20 21 22 30 31 32 33
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
 -9.934  -1.049  -3.953      F
 -7.582  -4.783  -5.249      F
 -7.904  -5.930  -4.475      T-
-11.235  -1.526  -3.440      T-
 -9.096  -0.359  -2.953      T-
 -9.138  -2.341  -4.451      T-
 -9.774  -3.363  -5.201      F
 -8.773  -4.201  -6.013      F
 -8.171  -3.398  -7.009      T-
 -7.040  -4.103  -7.507      F     1
 -6.651  -5.105  -6.416      F
 -5.931  -3.167  -7.804      F     1
 -5.533  -2.077  -7.079      F     1
 -4.528  -1.425  -7.591      T-    1
 -4.193  -2.175  -8.717      F     1
 -3.136  -2.004  -9.662      F     1
 -2.337  -1.074  -9.747      T-    1
 -3.069  -3.041 -10.570      T+ 2  1
 -3.936  -4.095 -10.621      F     1
 -3.724  -5.033 -11.534      T+ 3  1
 -4.977  -4.236  -9.788      T-    1
 -5.033  -3.256  -8.845      F     1
-10.481  -2.922  -5.906      F 
-10.322  -4.011  -4.518      F
 -9.306  -5.016  -6.507      F
 -7.142  -3.994  -4.636      F
 -5.605  -5.016  -6.117      F
 -6.837  -6.114  -6.787      F
 -7.328  -4.644  -8.408      F
 -6.015  -1.777  -6.161      F     1
 -2.312  -3.015 -11.248      F     1
 -2.930  -4.971 -12.162      F     1
 -4.356  -5.827 -11.600      F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.431 0.013
1.593 0.010 
1.440 0.016 
1.511 0.008 
1.528 0.010 
1.446 0.011 
1.420 0.013 
1.521 0.014 
1.518 0.010 
1.459 0.009  
1.374 0.007 
1.375 0.008 
1.305 0.006 
1.388 0.006 
1.379 0.007  
1.350 0.007  
1.323 0.008 
1.373 0.008 
1.419 0.010 
1.391 0.007  
1.237 0.009 
1.341 0.010 
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  114.7  1.5
 2  8  9  105.6  1.0
 2 11 10  102.7  1.4
 3  2  8  110.3  2.2
 3  2 11  110.6  2.7
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.4  1.6
 8  9 10  109.7  1.4
 8  2 11  103.2  1.0
 9 10 11  106.1  1.1
 9 10 12  107.8  0.8
11 10 12  114.2  1.6
10 12 13  127.0  1.3
10 12 22  126.5  1.3
12 13 14  111.1  0.5
12 22 15  105.4  0.4
12 22 21  126.0  0.6
13 14 15  104.3  0.5
13 12 22  106.4  0.4
14 15 16  130.4  0.6
14 15 22  110.8  0.4
15 16 17  128.6  0.6
15 16 18  111.5  0.5
15 22 21  128.6  0.5
16 18 19  125.1  0.6
16 15 22  118.8  0.6
17 16 18  119.9  0.6
18 19 20  116.2  0.9
18 19 21  123.9  0.6
19 21 22  111.9  0.5
20 19 21  119.9  0.7
  2  8 25  109.5    5.0
  2 11 27  109.5    5.0
  2 11 28  109.5    5.0
  3  2 26  109.5    5.0
  6  7 23  109.5    5.0
  6  7 24  109.5    5.0
  7  8 25  109.5    5.0
  8  2 26  109.5    5.0
  8  7 23  109.5    5.0
  8  7 24  109.5    5.0
  9  8 25  109.5    5.0
  9 10 29  109.5    5.0
 10 11 27  109.5    5.0
 10 11 28  109.5    5.0
 11  2 26  109.5    5.0
 11 10 29  109.5    5.0
 12 10 29  109.5    5.0
 12 13 30  120.0    5.0
 14 13 30  120.0    5.0
 16 18 31  120.0    5.0
 19 18 31  120.0    5.0
 19 20 32  120.0    5.0
 19 20 33  120.0    5.0
 23  7 24  109.5    5.0
 27 11 28  109.5    5.0
 32 20 33  120.0    5.0
#===== 1-3 connections
2   7
7   7   9  10  25  27  28
4   8  11  26
3   5   6
3   4   6
6   4   5   8  23  24
5   1   2   9  25
8   3   6  10  11  23  24  26
7   2   7  11  12  25  29
7   2   8  13  22  27  28
7   3   8   9  12  26  29
8   9  11  14  15  21  29  30
4  10  15  22
5  12  16  22  30
6  12  13  17  18  21
5  14  19  22  31
3  15  18
5  15  17  20  21
6  16  22  31  32  33
3  18  21
5  12  15  18  20
6  10  13  14  16  19
4   6   8  24
4   6   8  23
4   2   7   9
4   3   8  11
4   2  10  28
4   2  10  27
4   9  11  12
3  12  14
3  16  19
3  19  33
3  19  32
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
DCYTC  289.180
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  30  31   4   4   1  31   0  19  54   1   0   1  28    3
#===== Atom names
#---++++----++++----++++----++++----++++----++++----++++----
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N1  C6  C5  
 C4  N4  N3  C2  O2  H5''H5' H4' H3' H2''H2' H1' H41 H42 
 H5  H6 
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 18  13 14  18 19  17 18  15 17  14 15  15 16
7 20   7 21   8 22   2 23  11 24  11 25  10 26  16 27  16 28  14 29
13 30
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 22 23 24 25 26 27 28 29 30
 10 12   13 14 15 16 17 18 19 27 28 29 30
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
-6.853  -6.502  -3.404       F
-3.163  -8.690  -4.573       F
-2.932 -10.080  -4.348       T-
-7.524  -7.589  -2.635       T-
-6.396  -5.342  -2.624       T-
-5.591  -7.114  -4.151       T-
-5.647  -8.290  -4.936       F
-4.281  -8.489  -5.602       F
-3.943  -7.335  -6.344       T-
-2.535  -7.274  -6.477       F     1
-1.984  -7.925  -5.193       F
-2.146  -5.839  -6.624       F     1
-2.565  -4.932  -5.683       F     1
-2.264  -3.619  -5.777       F     1
-1.496  -3.212  -6.907       F     1
-1.203  -1.935  -7.069       T+ 3  1
-1.142  -4.071  -7.868       T-    1
-1.463  -5.395  -7.763       F     1
-1.142  -6.155  -8.682       T-    1
-6.418  -8.200  -5.702       F
-5.857  -9.153  -4.306       F
-4.310  -9.342  -6.284       F
-3.433  -8.194  -3.640       F
-1.651  -7.171  -4.480       F
-1.138  -8.574  -5.412       F
-2.249  -7.884  -7.334       F
-1.480  -1.243  -6.384       F     1
-0.716  -1.638  -7.901       F     1
-2.613  -2.931  -5.018       F     1
-3.150  -5.279  -4.843       F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.431 0.013
1.593 0.010 
1.440 0.016  
1.511 0.008 
1.528 0.010 
1.446 0.011 
1.420 0.013 
1.521 0.014 
1.518 0.010 
1.470 0.012   
1.367 0.006  
1.397 0.010 
1.339 0.008  
1.240 0.009  
1.353 0.008  
1.335 0.007 
1.425 0.008 
1.335 0.009
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  114.7  1.5
 2  8  9  105.6  1.0
 2 11 10  102.7  1.4
 3  2  8  110.3  2.2
 3  2 11  110.6  2.7
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.4  1.6
 8  9 10  109.7  1.4
 8  2 11  103.2  1.0
 9 10 11  106.1  1.1
 9 10 12  107.8  0.8
11 10 12  114.2  1.6
10 12 13  120.8  1.2
10 12 18  118.8  1.1
12 13 14  121.0  0.5
12 18 17  119.2  0.7
12 18 19  118.9  0.6
13 14 15  117.4  0.5
13 12 18  120.3  0.4
14 15 16  120.2  0.7
14 15 17  121.9  0.4
15 17 18  119.9  0.5
16 15 17  118.0  0.7
17 18 19  121.9  0.7
  2  8 22  109.5    5.0
  2 11 24  109.5    5.0
  2 11 25  109.5    5.0
  3  2 23  109.5    5.0
  6  7 20  109.5    5.0
  6  7 21  109.5    5.0
  7  8 22  109.5    5.0
  8  2 23  109.5    5.0
  8  7 20  109.5    5.0
  8  7 21  109.5    5.0
  9  8 22  109.5    5.0
  9 10 26  109.5    5.0
 10 11 24  109.5    5.0
 10 11 25  109.5    5.0
 11  2 23  109.5    5.0
 11 10 26  109.5    5.0
 12 10 26  109.5    5.0
 12 13 30  120.0    5.0
 13 14 29  120.0    5.0
 14 13 30  120.0    5.0
 15 14 29  120.0    5.0
 15 16 27  120.0    5.0
 15 16 28  120.0    5.0
 20  7 21  109.5    5.0
 24 11 25  109.5    5.0
 27 16 28  120.0    5.0
#===== 1-3 connections
2   7
7   7   9  10  22  24  25
4   8  11  23
3   5   6
3   4   6
6   4   5   8  20  21
5   1   2   9  22
8   3   6  10  11  20  21  23
7   2   7  11  12  22  26
7   2   8  13  18  24  25
7   3   8   9  12  23  26
8   9  11  14  17  19  26  30
5  10  15  18  29
5  12  16  17  30
6  13  18  27  28  29
3  14  17
5  12  14  16  19
4  10  13  15
3  12  17
4   6   8  21
4   6   8  20
4   2   7   9
4   3   8  11
4   2  10  25
4   2  10  24
4   9  11  12
3  15  28
3  15  27
3  13  15
3  12  14
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
OADEA  329.210
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  33  35   4   4   1  35   0  22  60   1   0   1  35    3
#===== Atom names
#---++++----++++----++++----++++----++++----++++----++++----
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N9  C8  N7
 C5  C6  N6  N1  C2  N3  C4  O2' H5''H5' H4' H3' H2''HO2
 H1' H8  H61 H62 H2 
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 21  13 14  14 15  15 21  20 21  19 20  18 19
15 16 16 18  16 17  22 11   7 23   7 24   8 25   2 26  11 27  22 28
10 29 13 30  17 31  17 32  19 33
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 21 22 25 26 27 28 29 30 31 32 33
 10 12   13 14 15 16 17 18 19 20 21 30 31 32 33
#===== Torsion angles1
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
  -5.424   5.900  -8.545     F
  -3.936   1.475  -6.568     F
  -2.835   0.857  -7.246     T-
  -4.938   6.835  -9.584     T-
  -6.662   5.142  -8.842     T-
  -4.252   4.900  -8.168     T-
  -4.439   3.493  -8.043     F
  -3.917   2.992  -6.692     F
  -4.823   3.390  -5.632     T-
  -4.515   2.592  -4.495     F     1
  -3.954   1.266  -5.042     F
  -5.674   2.482  -3.567     F     1
  -6.997   2.758  -3.771     F     1
  -7.769   2.554  -2.733     T-    1
  -6.912   2.110  -1.762     F     1
  -7.177   1.734  -0.431     F     1
  -8.462   1.799   0.033     T+ 3  1
  -6.079   1.329   0.286     T-    1
  -4.871   1.316  -0.291     F     1
  -4.522   1.652  -1.532     T-    1
  -5.609   2.049  -2.235     F     1
  -2.663   1.045  -4.476     T  4
  -3.857   2.998  -8.851     F
  -5.501   3.179  -8.162     F
  -2.911   3.400  -6.459     F
  -4.909   1.119  -6.984     F
  -4.620   0.416  -4.794     F
  -2.356   0.222  -4.873     F
  -3.709   3.149  -3.974     F
  -7.312   3.117  -4.745     F     1
  -9.174   2.114  -0.599     F     1
  -8.692   1.541   0.972     F     1
  -4.051   0.982   0.345     F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.423 0.014 
1.593 0.010
1.440 0.016
1.510 0.013
1.524 0.011
1.453 0.012
1.414 0.012
1.528 0.010 
1.525 0.011 
1.462 0.010 
1.373 0.008 
1.374 0.006 
1.311 0.007 
1.388 0.006 
1.383 0.007 
1.344 0.006 
1.331 0.009 
1.339 0.009 
1.406 0.009 
1.351 0.007 
1.335 0.008 
1.413 0.013 
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  115.5  1.5
 2  8  9  105.5  1.4
 2 11 10  101.5  0.9
 3  2  8  110.6  2.6
 3  2 11  111.0  2.8
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.2  1.4
 8  9 10  109.6  0.9
 8  2 11  102.7  1.0
 9 10 11  106.4  1.4
 9 10 12  108.2  1.0
11 10 12  113.4  1.6
10 12 13  127.7  1.8
10 12 21  126.3  1.8
12 13 14  113.8  0.5
12 21 15  105.8  0.4
12 21 20  127.4  0.8
13 14 15  103.9  0.5
13 12 21  105.8  0.4
14 15 16  132.3  0.7
14 15 21  110.7  0.5
15 16 17  123.7  0.8
15 16 18  117.7  0.5
15 21 20  126.8  0.7
16 18 19  118.6  0.6
16 15 21  117.0  0.5
17 16 18  118.6  0.6
18 19 20  129.3  0.5
19 20 21  110.6  0.5
22 11 10  110.6  3.0
22 11 2   113.3  2.9
  2  8 25  109.5    5.0
  2 11 27  109.5    5.0
  3  2 26  109.5    5.0
  6  7 23  109.5    5.0
  6  7 24  109.5    5.0
  7  8 25  109.5    5.0
  8  2 26  109.5    5.0
  8  7 23  109.5    5.0
  8  7 24  109.5    5.0
  9  8 25  109.5    5.0
  9 10 29  109.5    5.0
 10 11 27  109.5    5.0
 11  2 26  109.5    5.0
 11 10 29  109.5    5.0
 11 22 28  109.5    5.0
 12 10 29  109.5    5.0
 12 13 30  120.0    5.0
 14 13 30  120.0    5.0
 16 17 31  120.0    5.0
 16 17 32  120.0    5.0
 18 19 33  120.0    5.0
 20 19 33  120.0    5.0
 22 11 27  109.5    5.0
 23  7 24  109.5    5.0
 31 17 32  120.0    5.0
#===== 1-3 connections
2   7
7   7   9  10  22  25  27
4   8  11  26
3   5   6
3   4   6
6   4   5   8  23  24
5   1   2   9  25
8   3   6  10  11  23  24  26
7   2   7  11  12  25  29
7   2   8  13  21  22  27
8   3   8   9  12  26  28  29
8   9  11  14  15  20  29  30
4  10  15  21
5  12  16  21  30
6  12  13  17  18  20
6  14  19  21  31  32
3  15  18
5  15  17  20  33
3  16  21
5  12  15  18  33
6  10  13  14  16  19
4   2  10  27
4   6   8  24
4   6   8  23
4   2   7   9
4   3   8  11
4   2  10  22
2  11
4   9  11  12
3  12  14
3  16  32
3  16  31
3  18  20
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
OURAU  320.190
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  30  31   4   4   1  31   0  20  53   1   0   1  29    3
#===== Atom names
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N1  C6  C5
 C4  O4  N3  C2  O2  O2' H5''H5' H4' H3' H2''HO2 H1' H3 
 H5  H6 
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 18  13 14  18 19  17 18  15 17  14 15  15 16 
11 20  7 21   7 22   8 23   2 24  11 25  20 26  10 27  17 28  14 29
13 30
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 23 24 25 26 27 28 29 30
 10 12   13 14 15 16 17 18 19 28 29 30
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
-16.988  11.091   4.047      F
-14.887  12.328  -0.306      F
-15.698  12.082  -1.457      T-
-18.395  10.668   3.871      T-
-15.945  10.046   3.950      T-
-16.678  12.250   3.000      T-
-16.358  11.831   1.670      F
-15.666  12.940   0.856      F
-14.636  13.638   1.645      T-
-13.391  13.626   0.929      F     1
-13.784  13.345  -0.533      F
-12.427  12.631   1.519      F     1
-12.804  11.603   2.391      F     1
-11.956  10.708   2.925      F     1
-10.545  10.766   2.608      F     1
 -9.702   9.983   3.054      T-    1
-10.185  11.817   1.722      T+ 2  1
-11.055  12.774   1.144      F     1
-10.627  13.648   0.382      T-    1
-14.257  14.525  -1.182      T  4
-17.288  11.529   1.138      F
-15.679  10.946   1.679      F
-16.390  13.708   0.516      F
-14.415  11.386   0.067      F
-12.939  12.926  -1.109      F
-14.651  14.221  -2.005      F
-12.948  14.638   0.996      F
 -9.207  11.890   1.476      F     1
-12.329   9.938   3.595      F     1
-13.858  11.544   2.640      F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.423 0.014 
1.593 0.010
1.440 0.016
1.510 0.013
1.524 0.011
1.453 0.012
1.414 0.012
1.528 0.010 
1.525 0.011 
1.469 0.009
1.375 0.009
1.381 0.009
1.337 0.009
1.219 0.009
1.373 0.007
1.380 0.009
1.431 0.009
1.232 0.008
1.413 0.013
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  115.5  1.5
 2  8  9  105.5  1.4
 2 11 10  101.5  0.9
 3  2  8  110.6  2.6
 3  2 11  111.0  2.8
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.2  1.4
 8  9 10  109.6  0.9
 8  2 11  102.7  1.0
 9 10 11  106.4  1.4
 9 10 12  108.2  1.0
11 10 12  113.4  1.6
10 12 13  121.2  1.4
10 12 18  117.7  1.2
12 13 14  122.7  0.5
12 18 17  114.9  0.6
12 18 19  122.8  0.7
13 14 15  119.7  0.6
13 12 18  121.0  0.6
14 15 16  125.9  0.6
14 15 17  114.6  0.6
15 17 18  127.0  0.6
16 15 17  119.4  0.7
17 18 19  122.2  0.7
20 11 10  110.6  3.0
20 11 2   113.3  2.9
  2  8 23  109.5    5.0
  2 11 25  109.5    5.0
  3  2 24  109.5    5.0
  6  7 21  109.5    5.0
  6  7 22  109.5    5.0
  7  8 23  109.5    5.0
  8  2 24  109.5    5.0
  8  7 21  109.5    5.0
  8  7 22  109.5    5.0
  9  8 23  109.5    5.0
  9 10 27  109.5    5.0
 10 11 25  109.5    5.0
 11  2 24  109.5    5.0
 11 10 27  109.5    5.0
 11 20 26  109.5    5.0
 12 10 27  109.5    5.0
 12 13 30  120.0    5.0
 13 14 29  120.0    5.0
 14 13 30  120.0    5.0
 15 14 29  120.0    5.0
 15 17 28  120.0    5.0
 18 17 28  120.0    5.0
 20 11 25  109.5    5.0
 21  7 22  109.5    5.0
#===== 1-3 connections
2   7
7   7   9  10  20  23  25
4   8  11  24
3   5   6
3   4   6
6   4   5   8  21  22
5   1   2   9  23
8   3   6  10  11  21  22  24
7   2   7  11  12  23  27
7   2   8  13  18  20  25
8   3   8   9  12  24  26  27
8   9  11  14  17  19  27  30
5  10  15  18  29
5  12  16  17  30
5  13  18  28  29
3  14  17
5  12  14  16  19
5  10  13  15  28
3  12  17
4   2  10  25
4   6   8  22
4   6   8  21
4   2   7   9
4   3   8  11
4   2  10  20
2  11
4   9  11  12
3  15  18
3  13  15
3  12  14
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
OGUAG  345.210
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  34  36   4   4   1  36   0  23  62   1   0   1  37    3
#===== Atom names
#---++++----++++----++++----++++----++++----++++----++++----
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N9  C8  N7  
 C5  C6  O6  N1  C2  N2  N3  C4  O2' H5''H5' H4' H3' H2'
'HO2 H1' H8  H1  H21 H22 
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 22  13 14  14 15  15 22  21 22  19 21  18 19
15 16 16 18  16 17  19 20  23 11   7 24   7 25   8 26   2 27  11 28
23 29 10 30  13 31  18 32  20 33  20 34
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 26 27 28 29 30 31 32 33 34
 10 12   13 14 15 16 17 18 19 20 21 22 31 32 33 34
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
 -2.870  -0.714  -7.516      F
 -3.648  -4.002  -3.761      F
 -2.647  -5.016  -3.696      T-
 -1.826  -1.054  -8.508      T-
 -4.260  -1.116  -7.829      T-
 -2.434  -1.317  -6.112      T-
 -2.910  -2.597  -5.706      F
 -3.029  -2.683  -4.176      F
 -3.985  -1.708  -3.689      T-
 -4.559  -2.162  -2.466      F     1
 -4.340  -3.698  -2.422      F
 -5.989  -1.729  -2.415      F     1
 -6.780  -1.217  -3.409      F     1
 -8.005  -0.932  -3.045      T-    1
 -8.056  -1.269  -1.721      F     1
 -9.142  -1.167  -0.836      F     1
-10.265  -0.740  -1.126      T-    1
 -8.790  -1.622   0.465      T+ 2  1
 -7.533  -2.116   0.853      F     1
 -7.473  -2.481   2.159      T+ 3  1
 -6.511  -2.198  -0.038      T-    1
 -6.824  -1.769  -1.286      F     1
 -3.575  -4.064  -1.278      T  4
 -2.193  -3.373  -6.053      F
 -3.906  -2.835  -6.147      F
 -2.060  -2.494  -3.670      F
 -4.446  -4.238  -4.505      F
 -5.304  -4.239  -2.361      F
 -3.488  -5.022  -1.324      F
 -3.984  -1.663  -1.657      F 
 -6.357  -1.083  -4.406      F     1
 -9.527  -1.575   1.154      F     1
 -8.282  -2.383   2.742      F     1
 -6.619  -2.845   2.529      F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.423 0.014 
1.593 0.010
1.440 0.016
1.510 0.013
1.524 0.011
1.453 0.012
1.414 0.012
1.528 0.010 
1.525 0.011 
1.459 0.009  
1.374 0.007 
1.375 0.008 
1.305 0.006 
1.388 0.006 
1.379 0.007  
1.350 0.007  
1.323 0.008 
1.373 0.008 
1.419 0.010 
1.391 0.007  
1.237 0.009 
1.341 0.010 
1.413 0.013
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  115.5  1.5
 2  8  9  105.5  1.4
 2 11 10  101.5  0.9
 3  2  8  110.6  2.6
 3  2 11  111.0  2.8
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.2  1.4
 8  9 10  109.6  0.9
 8  2 11  102.7  1.0
 9 10 11  106.4  1.4
 9 10 12  108.2  1.0
11 10 12  113.4  1.6
10 12 13  127.0  1.3
10 12 22  126.5  1.3
12 13 14  111.1  0.5
12 22 15  105.4  0.4
12 22 21  126.0  0.6
13 14 15  104.3  0.5
13 12 22  106.4  0.4
14 15 16  130.4  0.6
14 15 22  110.8  0.4
15 16 17  128.6  0.6
15 16 18  111.5  0.5
15 22 21  128.6  0.5
16 18 19  125.1  0.6
16 15 22  118.8  0.6
17 16 18  119.9  0.6
18 19 20  116.2  0.9
18 19 21  123.9  0.6
19 21 22  111.9  0.5
20 19 21  119.9  0.7
23 11 10  110.6  3.0
23 11 2   113.3  2.9
  2  8 26  109.5    5.0
  2 11 28  109.5    5.0
  3  2 27  109.5    5.0
  6  7 24  109.5    5.0
  6  7 25  109.5    5.0
  7  8 26  109.5    5.0
  8  2 27  109.5    5.0
  8  7 24  109.5    5.0
  8  7 25  109.5    5.0
  9  8 26  109.5    5.0
  9 10 30  109.5    5.0
 10 11 28  109.5    5.0
 11  2 27  109.5    5.0
 11 10 30  109.5    5.0
 11 23 29  109.5    5.0
 12 10 30  109.5    5.0
 12 13 31  120.0    5.0
 14 13 31  120.0    5.0
 16 18 32  120.0    5.0
 19 18 32  120.0    5.0
 19 20 33  120.0    5.0
 19 20 34  120.0    5.0
 23 11 28  109.5    5.0
 24  7 25  109.5    5.0
 33 20 34  120.0    5.0
#===== 1-3 connections
2   7
7   7   9  10  23  26  28
4   8  11  27
3   5   6
3   4   6
6   4   5   8  24  25
5   1   2   9  26
8   3   6  10  11  24  25  27
7   2   7  11  12  26  30
7   2   8  13  22  23  28
8   3   8   9  12  27  29  30
8   9  11  14  15  21  30  31
4  10  15  22
5  12  16  22  31
6  12  13  17  18  21
5  14  19  22  32
3  15  18
5  15  17  20  21
6  16  22  32  33  34
3  18  21
5  12  15  18  20
6  10  13  14  16  19
4   2  10  28
4   6   8  25
4   6   8  24
4   2   7   9
4   3   8  11
4   2  10  23
2  11
4   9  11  12
3  12  14
3  16  19
3  19  34
3  19  33
#===== Hydrogen bond donor groups
#===== Parameter groups
#============================================================================
*
OCYTC  289.180
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  31  32   4   4   1  32   0  20  55   1   0   1  30    3
#===== Atom names
 P   C3' O3' OP1 OP2 O5' C5' C4' O4' C1' C2' N1  C6  C5  
 C4  N4  N3  C2  O2  O2' H5''H5' H4' H3' H2''HO2 H1' H41
 H42 H5  H6  
#===== Bonds
1 4    1 5    2 3    1 6    6 7    7 8    8 2    8 9    9 10  10 11
2 11  10 12  12 13  12 18  13 14  18 19  17 18  15 17  14 15  15 16
20 11  7 21   7 22   8 23   2 24  11 25  20 26  10 27  16 28  16 29
14 30 13 31
#===== Rotatable bonds
  1  6    2  3  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
  6  7    2  3  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
  7  8    2  3  9 10 11 12 13 14 15 16 17 18 19 20 23 24 25 26 27 28 29 30 31
 10 12   13 14 15 16 17 18 19 28 29 30 31
#===== Torsion angles
4 1 6 7    1 6 7 8    6 7 8 9    9 10 12 13
#===== Coordinates, hydrogen bonding and planarity
-12.650   6.979  15.519      F
-14.027  10.306  11.898      F
-15.452  10.376  11.977      T-
-12.856   6.974  16.985      T-
-13.458   6.032  14.719      T-
-12.874   8.463  14.989      T-
-13.528   8.662  13.737      F
-13.412  10.120  13.273      F
-12.018  10.501  13.070      T-
-11.938  11.499  12.040      F     1
-13.360  11.597  11.420      F
-10.851  11.170  11.056      F     1
-10.716   9.902  10.496      F     1
 -9.747   9.603   9.613      F     1
 -8.818  10.640   9.238      F     1
 -7.827  10.373   8.345      T+ 3  1
 -8.931  11.879   9.773      T-    1
 -9.914  12.205  10.677      F     1
 -9.961  13.368  11.131      T-    1
-14.033  12.759  11.903      T  4 
-14.608   8.411  13.843      F 
-13.103   8.004  12.946      F
-13.851  10.820  14.016      F
-13.691   9.462  11.247      F
-13.325  11.642  10.313      F
-14.956  12.628  11.664      F
-11.703  12.454  12.548      F
 -7.201  11.122   8.107      F     1
 -7.726   9.464   7.937      F     1
 -9.648   8.608   9.178      F     1
-11.449   9.165  10.819      F     1
#===== Standard bonds
1.485 0.017
1.485 0.017
1.423 0.014 
1.593 0.010
1.440 0.016
1.510 0.013
1.524 0.011
1.453 0.012
1.414 0.012
1.528 0.010 
1.525 0.011 
1.470 0.012 
1.367 0.006  
1.397 0.010 
1.339 0.008  
1.240 0.009  
1.353 0.008  
1.335 0.007 
1.425 0.008 
1.335 0.009
1.413 0.013
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
1.000 0.05
#===== Standard angles
 1  6  7  120.9  1.6
 2  8  7  115.5  1.5
 2  8  9  105.5  1.4
 2 11 10  101.5  0.9
 3  2  8  110.6  2.6
 3  2 11  111.0  2.8
 4  1  5  119.6  1.5
 4  1  6  108.1  2.9
 5  1  6  108.3  2.7
 6  7  8  110.2  1.4
 7  8  9  109.2  1.4
 8  9 10  109.6  0.9
 8  2 11  102.7  1.0
 9 10 11  106.4  1.4
 9 10 12  108.2  1.0
11 10 12  113.4  1.6
10 12 13  120.8  1.2
10 12 18  118.8  1.1
12 13 14  121.0  0.5
12 18 17  119.2  0.7
12 18 19  118.9  0.6
13 14 15  117.4  0.5
13 12 18  120.3  0.4
14 15 16  120.2  0.7
14 15 17  121.9  0.4
15 17 18  119.9  0.5
16 15 17  118.0  0.7
17 18 19  121.9  0.7
20 11 10  110.6  3.0
20 11 2   113.3  2.9
  2  8 23  109.5    5.0
  2 11 25  109.5    5.0
  3  2 24  109.5    5.0
  6  7 21  109.5    5.0
  6  7 22  109.5    5.0
  7  8 23  109.5    5.0
  8  2 24  109.5    5.0
  8  7 21  109.5    5.0
  8  7 22  109.5    5.0
  9  8 23  109.5    5.0
  9 10 27  109.5    5.0
 10 11 25  109.5    5.0
 11  2 24  109.5    5.0
 11 10 27  109.5    5.0
 11 20 26  109.5    5.0
 12 10 27  109.5    5.0
 12 13 31  120.0    5.0
 13 14 30  120.0    5.0
 14 13 31  120.0    5.0
 15 14 30  120.0    5.0
 15 16 28  120.0    5.0
 15 16 29  120.0    5.0
 20 11 25  109.5    5.0
 21  7 22  109.5    5.0
 28 16 29  120.0    5.0
#===== 1-3 connections
2   7
7   7   9  10  20  23  25
4   8  11  24
3   5   6
3   4   6
6   4   5   8  21  22
5   1   2   9  23
8   3   6  10  11  21  22  24
7   2   7  11  12  23  27
7   2   8  13  18  20  25
8   3   8   9  12  24  26  27
8   9  11  14  17  19  27  31
5  10  15  18  30
5  12  16  17  31
6  13  18  28  29  30
3  14  17
5  12  14  16  19
4  10  13  15
3  12  17
4   2  10  25
4   6   8  22
4   6   8  21
4   2   7   9
4   3   8  11
4   2  10  20
2  11
4   9  11  12
3  15  29
3  15  28
3  13  15
3  12  14
#===== Hydrogen bond donor groups
#===== Parameter groups
*
TNT 1
  35  36  10  10   1  36   0  23  55   1   1   1  31   6
H1N1N1  H2N1H3N1C7  N2  H2N2H3N2H1N2C4  C5  C6  C3  C2  
C1  O1  CA  CB  CA' O1' C1' C6' C5' C2' C3' C4' C7' N2' 
H2N2H3N2H1N2N1' H2N1H3N1H1N1
   1  2   2  3   2  4   2  5   5  6   5 10   6  7   6  8   6  9  10 11
  10 13  11 12  12 15  13 14  14 15  15 16  16 17  17 18  18 19  19 20
  20 21  21 22  21 24  22 23  23 26  24 25  25 26  26 27  27 28  27 32
  28 29  28 30  28 31  32 33  32 34  32 35
   2   5   1   3   4
   5  10   1   2   3   4   6   7   8   9
  15  16   1   2   3   4   5   6   7   8   9  10  11  12  13  14
  16  17   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
  17  18   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
  18  19   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
  19  20   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
  20  21   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
  26  27   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
  27  32   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  28  29  30  31
   1   2   5  10   2   5  10  13  12  15  16  17  15  16  17  18  16  17  18  19
  17  18  19  20  18  19  20  21  19  20  21  24  23  26  27  32  26  27  32  35
 13.75100 25.21000 64.73499  F  0
 13.73100 25.19000 64.71500  T+ 0
 13.77100 25.23000 64.75500  F  0
 13.79100 25.25000 64.77499  F  0
 13.01700 26.27700 65.10100  F  0
 12.89700 27.33600 64.26300  T+ 0
 12.93700 27.37600 64.30300  F  0
 12.95700 27.39600 64.32300  F  0
 12.91700 27.35600 64.28300  F  0
 12.32700 26.31300 66.34400  F  0
 11.31700 27.26800 66.57600  F  0
 10.56400 27.21100 67.74700  F  0
 12.55800 25.33300 67.32200  F  0
 11.79800 25.26100 68.48700  F  0
 10.78100 26.21700 68.70500  F  0
  9.93100 26.18400 69.76000  T- 0
  9.88600 24.94500 70.48400  F  0
  8.72000 25.01300 71.46700  F  0
  8.95800 24.07000 72.63600  F  0
  7.69100 23.86700 73.32600  T- 0
  7.76000 23.09800 74.47100  F  0
  8.99700 22.60900 74.92900  F  0
  9.06200 21.85500 76.06800  F  0
  6.58700 22.78100 75.21000  F  0
  6.63600 22.01700 76.39500  F  0
  7.94300 21.55000 76.82100  F  0
  8.18200 20.78200 77.98100  F  0
  7.13700 20.05600 78.48900  T+ 0
  7.17700 20.09600 78.52900  F  0
  7.19700 20.11600 78.54900  F  0
  7.15700 20.07600 78.50900  F  0
  9.36600 20.86900 78.62300  T+ 0
  9.40600 20.90900 78.66300  F  0
  9.42600 20.92900 78.68300  F  0
  9.38600 20.88900 78.64300  F  0
1.000 0.050
1.000 0.050
1.000 0.050
1.330 0.020
1.330 0.020
1.390 0.020
1.000 0.050
1.000 0.050
1.000 0.050
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.430 0.020
1.430 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.430 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.330 0.020
1.330 0.020
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
 2  5  6  120.0  2.0
 2  5 10  120.0  2.0
 6  5 10  120.0  2.0
 5 10 11  120.0  2.0
 5 10 13  120.0  2.0
11 10 13  120.0  2.0
10 11 12  120.0  2.0
11 12 15  120.0  2.0
10 13 14  120.0  2.0
13 14 15  120.0  2.0
12 15 14  120.0  2.0
12 15 16  120.0  2.0
14 15 16  120.0  2.0
15 16 17  122.0  2.0
16 17 18  109.5  2.0
17 18 19  111.0  2.0
18 19 20  109.5  2.0
19 20 21  122.0  2.0
20 21 22  120.0  2.0
20 21 24  120.0  2.0
22 21 24  120.0  2.0
21 22 23  120.0  2.0
22 23 26  120.0  2.0
21 24 25  120.0  2.0
24 25 26  120.0  2.0
23 26 25  120.0  2.0
23 26 27  120.0  2.0
25 26 27  120.0  2.0
26 27 28  120.0  2.0
26 27 32  120.0  2.0
28 27 32  120.0  2.0
 1  2  3  109.5  2.0
 1  2  4  109.5  2.0
 3  2  4  109.5  2.0
 1  2  5  109.5  2.0
 3  2  5  109.5  2.0
 4  2  5  109.5  2.0
 5  6  7  109.5  2.0
 5  6  8  109.5  2.0
 7  6  8  109.5  2.0
 5  6  9  109.5  2.0
 7  6  9  109.5  2.0
 8  6  9  109.5  2.0
27 28 29  109.5  2.0
27 28 30  109.5  2.0
29 28 30  109.5  2.0
27 28 31  109.5  2.0
29 28 31  109.5  2.0
30 28 31  109.5  2.0
27 32 33  109.5  2.0
27 32 34  109.5  2.0
33 32 34  109.5  2.0
27 32 35  109.5  2.0
33 32 35  109.5  2.0
34 32 35  109.5  2.0
   4   3   4   5
   3   6  10
   4   1   4   5
   4   1   3   5
   9   1   3   4   7   8   9  11  13
   3   2  10
   4   5   8   9
   4   5   7   9
   4   5   7   8
   5   2   6  12  14
   4   5  13  15
   4  10  14  16
   4   5  11  15
   4  10  12  16
   4  11  13  17
   4  12  14  18
   3  15  19
   3  16  20
   3  17  21
   4  18  22  24
   4  19  23  25
   4  20  24  26
   4  21  25  27
   4  20  22  26
   4  21  23  27
   5  22  24  28  32
   9  23  25  29  30  31  33  34  35
   3  26  32
   4  27  30  31
   4  27  29  31
   4  27  29  30
   3  26  28
   4  27  34  35
   4  27  33  35
   4  27  33  34
0.000 0.000 0.000 0.000 0.127 0.129 0.248 0.248 0.248 0.000
0.000 0.000 0.000 0.000 0.180 -.360 0.180 0.000 0.180 -.360
0.180 0.000 0.000 0.000 0.000 0.000 0.127 0.129 0.248 0.248
0.248 0.000 0.000 0.000 0.000
#
# A tyrosine-like attached group. This is the protonated form
# Added by Jens
*
HC4  2
  12  12   3   3   1  12   0  11  14   1   1   1  13   6
C2' C3' C4' O4' HA  C5' C6' C1' C3  C2  C1  O1  
   1   2   1   8   2   3   3   4   3   6   4   5   6   7   7   8   8   9   9  10
  10  11  11  12
   8   9   1   2   3   4   5   6   7
   9  10   1   2   3   4   5   6   7   8
  10  11   1   2   3   4   5   6   7   8   9
   1   8   9  10   8   9  10  11   9  10  11  12
  0.67400 -0.90300 -0.96000  F  0
  1.73100 -1.59700 -0.38000  F  0
  1.93700 -1.49400  0.97700  F  0
  3.09500 -2.20200  1.61400  T  4
  3.59684 -2.71480  0.91744  F  0
  1.10900 -0.72400  1.76700  F  0
  0.06000 -0.04400  1.18500  F  0
 -0.17600 -0.11700 -0.18500  F  0
 -1.23800  0.71100 -0.90200  F  0
 -2.01800  1.50900 -0.22700  F  0
 -2.95000  2.41200 -0.99700  F  0
 -2.77100  2.75400 -2.17200  T- 0
1.390 0.020
1.390 0.020
1.390 0.020
1.360 0.020
1.390 0.020
1.000 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.230 0.020
 2  1  8  120.0  2.0
 1  2  3  120.0  2.0
 2  3  4  120.0  2.0
 2  3  6  120.0  2.0
 4  3  6  120.0  2.0
 3  6  7  120.0  2.0
 6  7  8  120.0  2.0
 1  8  7  120.0  2.0
 1  8  9  120.0  2.0
 7  8  9  120.0  2.0
 8  9 10  120.0  2.0
 9 10 11  120.0  2.0
10 11 12  121.0  2.0
 3  4  5  109.5  2.0
   4   3   7   9
   4   8   4   6
   4   1   5   7
   3   2   6
   2   3
   4   2   4   8
   4   3   1   9
   4   2   6  10
   4   1   7  11
   3   8  12
   2   9
   2  10
0.000 0.000 0.150 -.548 0.398 0.000 0.000 0.000 0.000 0.000
0.380 -.380
#
# Here is the deprotonated form
# Jens 29/11-2000
*
HCD 4
  10  10   3   3   1  10   0  10  11   1   1   1  11   6
C2' C3' C4' C5' C6' C1' C3  C2  C1  O1  
   1   2   1   6   2   3   3   4   4   5   5   6   6   7   7   8   8   9   9  10
   6   7   1   2   3   4   5
   7   8   1   2   3   4   5   6
   8   9   1   2   3   4   5   6   7
   1   6   7   8   6   7   8   9   7   8   9  10
 15.72000 -0.40300-19.66000  F  0
 16.77700 -1.09700-19.08000  F  0
 16.98300 -0.99400-17.72300  F  0
 16.15500 -0.22400-16.93300  F  0
 15.10600  0.45600-17.51500  F  0
 14.87000  0.38300-18.88500  F  0
 13.83400  1.19000-19.58400  F  0
 13.05400  1.98800-18.90900  F  0
 12.12200  2.89100-19.67900  F  0
 12.30100  3.23300-20.85400  T- 0
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.390 0.020
1.230 0.020
 2  1  6  120.0  2.0
 1  2  3  120.0  2.0
 2  3  4  120.0  2.0
 3  4  5  120.0  2.0
 4  5  6  120.0  2.0
 1  6  5  120.0  2.0
 1  6  7  120.0  2.0
 5  6  7  120.0  2.0
 6  7  8  120.0  2.0
 7  8  9  120.0  2.0
 8  9 10  121.0  2.0
   4   3   5   7
   3   6   4
   3   1   5
   3   2   6
   4   3   1   7
   4   2   4   8
   4   1   5   9
   3   6  10
   2   7
   2   8
0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.380 -.380
# test molecule (was Gly-Asp)
*
DRG 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
   12  11  2   2   1   11  1   8  17   1   2   1  8    6
#===== Atom names
#...||||....||||....||||....||||....||||....||||....||||....||||
 N   CA  C   O   CB  CG  OD1 OD2 H   HA 2HB 3HB 
#===== Bonds
1 2  2 3  3 4  2 5  5 6  6 7  6 8  1 9  2 10  5 11
5 12 
#===== Rotatable bonds
  2  5    6  7  8 11 12
  5  6    7  8
#===== Torsion angles
1 2 5 6   2 5 6 7
#===== Coordinates, hydrogen bonding and planarity
  1.20131  0.84661  0.00001  T+ 1
  0.01001  0.00001  0.00001  F
 -1.25031  0.88111  0.00001  F
 -2.18531  0.66031 -0.78411  T-
  0.04341 -1.01681  1.12341  F     1
  1.22931 -1.96551  1.05271  F     1
  2.05771 -1.86411  0.12721  T  7 11
  1.31911 -2.81041  1.95881  T  7 11
  1.20131  1.84661  0.00001  F         1
  0.00002 -0.51093 -0.85961  F
  0.08630 -0.52385  1.99241  F
 -0.79593 -1.55935  1.08923  F
#===== Standard bonds
1.458 0.019
1.525 0.021
1.231 0.020
1.530 0.020
1.516 0.025
1.249 0.019 
1.249 0.019
1.000 0.050
1.000 0.050
1.000 0.050
1.000 0.050
#===== Standard angles
 1  2  3  111.2  2.8
 2  3  4  120.8  1.7
 1  2  5  110.5  1.7
 3  2  5  110.1  1.9
 2  5  6  112.6  1.0
 5  6  8  118.4  2.3
 5  6  7  118.4  2.3
 8  6  7  122.9  2.4
  1  2 10  109.5    5.0
  2  1  9  120.0    5.0
  2  5 11  109.5    5.0
  2  5 12  109.5    5.0
  3  2 10  109.5    5.0
  5  2 10  109.5    5.0
  6  5 11  109.5    5.0
  6  5 12  109.5    5.0
 11  5 12  109.5    5.0
#===== 1-3 connections
   4   3   5  10
   6   9   4   6  11  12
   4   1   5  10
   2   2
   6   1   3  10   7   8
   4   2  11  12
   3   5   8
   3   5   7
   2   2
   4   1   3   5
   4   2   6  12
   4   2   6  11
#===== Hydrogen bond donor groups
    2 0.00 0.75 9.99
#===== Parameter groups
 0.1 0.0 0.38 -0.48 0.0 0.27 -0.635 -0.635
 0.1 0.0 0.38 -0.48 0.0 0.53 -0.18  -0.18 
*
ANP 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  39  41  11  11   1  41   0  31  65   1   1   1  52   6
#===== Atom names
O3G PG  O1G O2G N3B HAC PB  O1B O2B O3A PA  O1A O2A O5' 
C5' C4' O4' C3' O3' HAA C2' O2' HAB C1' N9  C8  N7  HAD 
C5  C4  N3  HAH C2  N1  HAG C6  N6  HAF HAE 
#===== Bonds
   1   2   2   3   2   4   2   5   5   6   5   7   7   8   7   9   7  10  10  11
  11  12  11  13  11  14  14  15  15  16  16  17  16  18  17  24  18  19  18  21
  19  20  21  22  21  24  22  23  24  25  25  26  25  30  26  27  27  28  27  29
  29  30  29  36  30  31  31  32  31  33  33  34  34  35  34  36  36  37  37  38
  37  39
#===== Rotatable bonds
   2   5   1   3   4
   5   7   1   2   3   4   6
   7  10   1   2   3   4   5   6   8   9
  10  11   1   2   3   4   5   6   7   8   9
  11  14   1   2   3   4   5   6   7   8   9  10  12  13
  14  15   1   2   3   4   5   6   7   8   9  10  11  12  13
  15  16   1   2   3   4   5   6   7   8   9  10  11  12  13  14
  18  19  20
  21  22  23
  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39
  36  37  38  39
#===== Torsion angles
   1   2   5   7   2   5   7  10   5   7  10  11   7  10  11  14  10  11  14  15
  11  14  15  16  14  15  16  18  16  18  19  20  18  21  22  23  17  24  25  30
  29  36  37  39
#===== Coordinates, hydrogen bonding and planarity
-23.63000 28.14000 14.12000  T-  
-23.26000 27.89000 12.57000  F   
-21.94000 27.21000 12.51000  T-  
-23.36000 29.21000 11.92000  T-  
-24.41000 26.87000 12.09000  T+ 5
-24.77000 26.40000 12.90000  F   
-25.67000 27.46000 11.30000  F   
-26.62000 26.33000 11.18000  T-  
-26.18000 28.69000 11.94000  T-  
-25.00000 27.80000  9.88000  T-  
-25.88000 28.04000  8.55000  F   
-26.70000 26.82000  8.37000  T-  
-26.60000 29.32000  8.70000  T-  
-24.73000 28.15000  7.43000  T-  
-23.88000 27.03000  7.17000  F   
-22.61000 27.42000  6.41000  F   
-22.83000 27.78000  5.01000  T-  
-21.90000 28.64000  7.02000  F   
-20.76000 28.20000  7.76000  T  4
-21.00000 28.14000  8.73000  F   
-21.48000 29.56000  5.85000  F   
-20.08000 29.50000  5.60000  T  4
-19.58000 29.71000  6.44000  F   
-22.21000 29.03000  4.59000  F   
-23.22000 30.02000  4.15000  F   
-24.37000 30.31000  4.83000  F   
-24.94000 31.28000  4.03000  T+ 2
-25.82000 31.70000  4.26000  F   
-24.20000 31.59000  2.92000  F   
-23.05000 30.75000  3.00000  F   
-22.07000 30.80000  2.00000  T+ 2
-21.27000 30.20000  2.08000  F   
-22.20000 31.68000  0.92000  F   
-23.33000 32.50000  0.85000  T+ 2
-23.44000 33.13000  0.09000  F   
-24.32000 32.45000  1.84000  F   
-25.36000 33.27000  1.70000  T+ 3
-25.43000 33.89000  0.92000  F   
-26.08000 33.25000  2.40000  F   
#===== Standard bond lengths
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.000 0.020
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.430 0.020
1.530 0.020
1.435 0.020
1.530 0.020
1.435 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.480 0.020
1.330 0.020
1.330 0.020
1.330 0.020
1.000 0.020
1.330 0.020
1.390 0.020
1.390 0.020
1.400 0.020
1.000 0.020
1.400 0.020
1.400 0.020
1.000 0.020
1.400 0.020
1.330 0.020
1.000 0.020
1.000 0.020
#===== Standard bond angles
 1  2  3  109.6  2.0
 1  2  4  109.6  2.0
 3  2  4  120.0  2.0
 1  2  5  103.0  2.0
 3  2  5  109.6  2.0
 4  2  5  109.6  2.0
 2  5  7  116.0  2.0
 5  7  8  109.6  2.0
 5  7  9  109.6  2.0
 8  7  9  120.0  2.0
 5  7 10  103.0  2.0
 8  7 10  109.6  2.0
 9  7 10  109.6  2.0
 7 10 11  120.0  2.0
10 11 12  109.6  2.0
10 11 13  109.6  2.0
12 11 13  120.0  2.0
10 11 14  103.0  2.0
12 11 14  109.6  2.0
13 11 14  109.6  2.0
11 14 15  120.0  2.0
14 15 16  111.0  2.0
15 16 17  109.5  2.0
15 16 18  109.5  2.0
17 16 18  104.0  2.0
16 17 24  104.0  2.0
16 18 19  109.5  2.0
16 18 21  104.0  2.0
19 18 21  109.5  2.0
18 21 22  109.5  2.0
18 21 24  104.0  2.0
22 21 24  109.5  2.0
17 24 21  104.0  2.0
17 24 25  109.5  2.0
21 24 25  109.5  2.0
24 25 26  126.0  2.0
24 25 30  126.0  2.0
26 25 30  108.0  2.0
25 26 27  108.0  2.0
26 27 29  108.0  2.0
27 29 30  108.0  2.0
27 29 36  132.0  2.0
30 29 36  120.0  2.0
25 30 29  108.0  2.0
25 30 31  132.0  2.0
29 30 31  120.0  2.0
30 31 33  120.0  2.0
31 33 34  120.0  2.0
33 34 36  120.0  2.0
29 36 34  120.0  2.0
29 36 37  120.0  2.0
34 36 37  120.0  2.0
 2  5  6  115.0  2.0
 6  5  7  115.0  2.0
18 19 20  109.5  2.0
21 22 23  109.5  2.0
26 27 28  126.0  2.0
28 27 29  126.0  2.0
30 31 32  120.0  2.0
32 31 33  120.0  2.0
33 34 35  120.0  2.0
35 34 36  120.0  2.0
36 37 38  120.0  2.0
36 37 39  120.0  2.0
38 37 39  120.0  2.0
   4   3   4   5
   3   6   7
   4   1   4   5
   4   1   3   5
   7   1   3   4   8   9  10
   3   2   7
   4   2   6  11
   4   5   9  10
   4   5   8  10
   7   5   8   9  12  13  14
   3   7  15
   4  10  13  14
   4  10  12  14
   5  10  12  13  16
   4  11  17  18
   5  14  24  19  21
   5  15  18  21  25
   6  15  17  20  22  24
   3  16  21
   2  18
   6  16  19  23  17  25
   3  18  24
   2  21
   6  16  18  22  26  30
   6  17  21  27  29  31
   5  24  30  28  29
   4  25  30  36
   3  26  29
   7  26  28  25  31  34  37
   7  24  26  27  36  32  33
   4  25  29  34
   3  30  33
   5  30  32  35  36
   4  31  29  37
   3  33  36
   7  27  30  33  35  38  39
   3  29  34
   3  36  39
   3  36  38
0.000 0.000 0.000 0.000 -.280 0.280 0.000 0.000 0.000 -.360
0.990 -.635 -.635 -.360 0.000 0.180 -.360 0.180 0.000 0.000
0.150 -.548 0.398 0.000 0.000 0.000 -.280 0.280 0.000 0.000
-.280 0.280 0.000 -.280 0.280 0.000 -.830 0.415 0.415
*
ATP 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  38  40  11  11   1  40   0  31  63   1   1   1  52   6
#===== Atom names
O3G P3  O1G O2G O3B P2  O1B O2B O3A P1  O1A O2A O5' C5' 
C4' O4' C3' O3' HAA C2' O2' HAB C1' N9  C8  N7  HAC C5  
C4  N3  HAG C2  N1  HAF C6  N6  HAE HAD 
#===== Bonds
   1   2   2   3   2   4   2   5   5   6   6   7   6   8   6   9   9  10  10  11
  10  12  10  13  13  14  14  15  15  16  15  17  16  23  17  18  17  20  18  19
  20  21  20  23  21  22  23  24  24  25  24  29  25  26  26  27  26  28  28  29
  28  35  29  30  30  31  30  32  32  33  33  34  33  35  35  36  36  37  36  38
#===== Rotatable bonds
   2   5   1   3   4
   5   6   1   2   3   4
   6   9   1   2   3   4   5   7   8
   9  10   1   2   3   4   5   6   7   8
  10  13   1   2   3   4   5   6   7   8   9  11  12
  13  14   1   2   3   4   5   6   7   8   9  10  11  12
  14  15   1   2   3   4   5   6   7   8   9  10  11  12  13
  17  18  19
  20  21  22
  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38
  35  36  37  38
#===== Torsion angles
   1   2   5   6   2   5   6   9   5   6   9  10   6   9  10  13   9  10  13  14
  10  13  14  15  13  14  15  17  15  17  18  19  17  20  21  22  16  23  24  29
  28  35  36  38
#===== Coordinates, hydrogen bonding and planarity
 16.48000  6.93000 -0.04000  T-  
 16.20000  8.40000 -0.66000  F   
 16.46000  8.37000 -2.14000  T-  
 17.10000  9.41000 -0.01000  T-  
 14.65000  8.81000 -0.40000  T-  
 14.37000 10.27000 -1.02000  F   
 15.52000 11.22000 -0.98000  T-  
 13.88000  9.93000 -2.37000  T-  
 13.13000 10.73000 -0.28000  T-  
 11.86000  9.90000 -0.17000  F   
 10.81000 10.09000 -1.20000  T-  
 12.19000  8.52000  0.24000  T-  
 11.17000 10.49000  1.14000  T-  
 11.93000 10.67000  2.34000  F   
 11.16000 10.19000  3.55000  F   
  9.95000 10.98000  3.71000  T-  
 10.66000  8.76000  3.53000  F   
 11.69000  7.85000  3.91000  T  4
 12.41000  7.85000  3.22000  F   
  9.53000  8.80000  4.55000  F   
  9.97000  8.78000  5.90000  T  4
 10.37000  7.88000  6.10000  F   
  8.92000 10.18000  4.27000  F   
  7.80000 10.15000  3.33000  F   
  7.85000 10.10000  1.96000  F   
  6.67000 10.10000  1.39000  T+ 2
  6.46000 10.07000  0.41000  F   
  5.78000 10.16000  2.45000  F   
  6.46000 10.18000  3.66000  F   
  5.92000 10.22000  4.89000  T+ 2
  6.45000 10.24000  5.74000  F   
  4.59000 10.22000  4.81000  F   
  3.80000 10.22000  3.73000  T+ 2
  2.80000 10.23000  3.82000  F   
  4.38000 10.19000  2.51000  F   
  3.59000 10.20000  1.43000  T+ 3
  2.60000 10.23000  1.53000  F   
  4.00000 10.18000  0.51000  F   
#===== Standard bond lengths
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.430 0.020
1.530 0.020
1.435 0.020
1.530 0.020
1.435 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.480 0.020
1.330 0.020
1.330 0.020
1.330 0.020
1.000 0.020
1.330 0.020
1.390 0.020
1.390 0.020
1.400 0.020
1.000 0.020
1.400 0.020
1.400 0.020
1.000 0.020
1.400 0.020
1.330 0.020
1.000 0.020
1.000 0.020
#===== Standard bond angles
 1  2  3  109.6  2.0
 1  2  4  109.6  2.0
 3  2  4  120.0  2.0
 1  2  5  103.0  2.0
 3  2  5  109.6  2.0
 4  2  5  109.6  2.0
 2  5  6  120.0  2.0
 5  6  7  109.6  2.0
 5  6  8  109.6  2.0
 7  6  8  120.0  2.0
 5  6  9  103.0  2.0
 7  6  9  109.6  2.0
 8  6  9  109.6  2.0
 6  9 10  120.0  2.0
 9 10 11  109.6  2.0
 9 10 12  109.6  2.0
11 10 12  120.0  2.0
 9 10 13  103.0  2.0
11 10 13  109.6  2.0
12 10 13  109.6  2.0
10 13 14  120.0  2.0
13 14 15  111.0  2.0
14 15 16  109.5  2.0
14 15 17  109.5  2.0
16 15 17  104.0  2.0
15 16 23  104.0  2.0
15 17 18  109.5  2.0
15 17 20  104.0  2.0
18 17 20  109.5  2.0
17 20 21  109.5  2.0
17 20 23  104.0  2.0
21 20 23  109.5  2.0
16 23 20  104.0  2.0
16 23 24  109.5  2.0
20 23 24  109.5  2.0
23 24 25  126.0  2.0
23 24 29  126.0  2.0
25 24 29  108.0  2.0
24 25 26  108.0  2.0
25 26 28  108.0  2.0
26 28 29  108.0  2.0
26 28 35  132.0  2.0
29 28 35  120.0  2.0
24 29 28  108.0  2.0
24 29 30  132.0  2.0
28 29 30  120.0  2.0
29 30 32  120.0  2.0
30 32 33  120.0  2.0
32 33 35  120.0  2.0
28 35 33  120.0  2.0
28 35 36  120.0  2.0
33 35 36  120.0  2.0
17 18 19  109.5  2.0
20 21 22  109.5  2.0
25 26 27  126.0  2.0
27 26 28  126.0  2.0
29 30 31  120.0  2.0
31 30 32  120.0  2.0
32 33 34  120.0  2.0
34 33 35  120.0  2.0
35 36 37  120.0  2.0
35 36 38  120.0  2.0
37 36 38  120.0  2.0
   4   3   4   5
   2   6
   4   1   4   5
   4   1   3   5
   7   1   3   4   7   8   9
   3   2  10
   4   5   8   9
   4   5   7   9
   7   5   7   8  11  12  13
   3   6  14
   4   9  12  13
   4   9  11  13
   5   9  11  12  15
   4  10  16  17
   5  13  23  18  20
   5  14  17  20  24
   6  14  16  19  21  23
   3  15  20
   2  17
   6  15  18  22  16  24
   3  17  23
   2  20
   6  15  17  21  25  29
   6  16  20  26  28  30
   5  23  29  27  28
   4  24  29  35
   3  25  28
   7  25  27  24  30  33  36
   7  23  25  26  35  31  32
   4  24  28  33
   3  29  32
   5  29  31  34  35
   4  30  28  36
   3  32  35
   7  26  29  32  34  37  38
   3  28  33
   3  35  38
   3  35  37
-.360 0.990 -.635 -.635 -.360 0.000 0.000 0.000 -.360 0.990
-.635 -.635 -.360 0.000 0.180 -.360 0.180 0.000 0.000 0.150
-.548 0.398 0.000 0.000 0.000 -.280 0.280 0.000 0.000 -.280
0.280 0.000 -.280 0.280 0.000 -.830 0.415 0.415
*
PO3 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
   4   3   0   0   1   3   0   4   3   1   1   1   3   6
#===== Atom names
OP3 P   OP2 OP1 
#===== Bonds
   1   2   2   3   2   4
#===== Rotatable bonds
#===== Torsion angles
#===== Coordinates, hydrogen bonding and planarity
-42.50000 53.37000 14.12000  T-  
-43.54000 53.15000 15.34000  F   
-44.61000 54.16000 15.22000  T-  
-42.72000 53.17000 16.57000  T-  
#===== Standard bond lengths
1.610 0.020
1.480 0.020
1.480 0.020
#===== Standard bond angles
 1  2  3  109.6  2.0
 1  2  4  109.6  2.0
 3  2  4  120.0  2.0
   3   3   4
   1
   3   1   4
   3   1   3
0.000 0.000 0.000 0.000
*
MPD 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  10   9   4   4   1   9   0   8  12   1   1   1  10   6
#===== Atom names
HAA OAH CAD CAC CAB CAE CAF CAG OAI HAB 
#===== Bonds
   1   2   2   3   3   4   3   5   3   6   6   7   7   8   7   9   9  10
#===== Rotatable bonds
   2   3   1
   3   6   7   8   9  10
   6   7   8   9  10
   7   9  10
#===== Torsion angles
   1   2   3   6   2   3   6   7   3   6   7   9   6   7   9  10
#===== Coordinates, hydrogen bonding and planarity
 -4.32000 -9.67000 -9.53000  F   
 -5.30000 -9.69000 -9.32000  T  4
 -5.64000 -8.56000 -8.64000  F   
 -6.71000 -7.79000 -9.41000  F   
 -6.20000 -8.98000 -7.28000  F   
 -4.43000 -7.66000 -8.39000  F   
 -3.31000 -7.69000 -9.44000  F   
 -1.86000 -7.68000 -8.98000  F   
 -3.60000 -7.74000-10.77000  T  7
 -4.59000 -7.78000-10.91000  F   
#===== Standard bond lengths
1.000 0.020
1.360 0.020
1.530 0.020
1.530 0.020
1.530 0.020
1.530 0.020
1.530 0.020
1.360 0.020
1.000 0.020
#===== Standard bond angles
 2  3  4  109.5  2.0
 2  3  5  109.5  2.0
 4  3  5  109.5  2.0
 2  3  6  109.5  2.0
 4  3  6  109.5  2.0
 5  3  6  109.5  2.0
 3  6  7  111.0  2.0
 6  7  8  120.0  2.0
 6  7  9  120.0  2.0
 8  7  9  120.0  2.0
 1  2  3  109.5  2.0
 7  9 10  109.5  2.0
   2   3
   4   4   5   6
   3   1   7
   4   2   5   6
   4   2   4   6
   6   2   4   5   8   9
   3   3  10
   3   6   9
   3   6   8
   2   7
0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.150 -.548 0.398
*
AT2 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  38  40  11  11   1  40   0  31  63   1   1   1  52   6
#===== Atom names
O2G PG  O3G O1G O3B PB  O1B O2B O3A PA  O1A O2A O5' C5' 
C4' O4' C3' O3' HAA C2' O2' HAB C1' N9  C8  N7  HAC C5  
C4  N3  HAD C2  N1  HAE C6  N6  HAG HAF 
#===== Bonds
   1   2   2   3   2   4   2   5   5   6   6   7   6   8   6   9   9  10  10  11
  10  12  10  13  13  14  14  15  15  16  15  17  16  23  17  18  17  20  18  19
  20  21  20  23  21  22  23  24  24  25  24  29  25  26  26  27  26  28  28  29
  28  35  29  30  30  31  30  32  32  33  33  34  33  35  35  36  36  37  36  38
#===== Rotatable bonds
   2   5   1   3   4
   5   6   1   2   3   4
   6   9   1   2   3   4   5   7   8
   9  10   1   2   3   4   5   6   7   8
  10  13   1   2   3   4   5   6   7   8   9  11  12
  13  14   1   2   3   4   5   6   7   8   9  10  11  12
  14  15   1   2   3   4   5   6   7   8   9  10  11  12  13
  17  18  19
  20  21  22
  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38
  35  36  37  38
#===== Torsion angles
   1   2   5   6   2   5   6   9   5   6   9  10   6   9  10  13   9  10  13  14
  10  13  14  15  13  14  15  17  15  17  18  19  17  20  21  22  16  23  24  29
  28  35  36  38
#===== Coordinates, hydrogen bonding and planarity
-23.36000 29.26000 12.10000  T-  
-23.15000 27.94000 12.76000  F   
-23.46000 27.93000 14.23000  T-  
-21.88000 27.27000 12.42000  T-  
-24.32000 26.99000 12.02000  T-  
-25.72000 27.43000 11.20000  F   
-26.67000 26.29000 11.31000  T-  
-26.22000 28.73000 11.66000  T-  
-25.16000 27.53000  9.73000  T-  
-25.95000 27.94000  8.40000  F   
-26.78000 26.83000  7.94000  T-  
-26.58000 29.29000  8.53000  T-  
-24.68000 28.10000  7.42000  T-  
-23.67000 27.08000  7.26000  F   
-22.49000 27.60000  6.47000  F   
-22.93000 27.82000  5.13000  T-  
-21.96000 28.97000  6.93000  F   
-20.95000 28.86000  7.99000  T  4
-21.37000 28.50000  8.82000  F   
-21.30000 29.51000  5.70000  F   
-19.96000 28.95000  5.50000  T  4
-19.35000 29.30000  6.21000  F   
-22.22000 28.96000  4.60000  F   
-23.18000 30.02000  4.17000  F   
-24.33000 30.42000  4.78000  F   
-24.98000 31.34000  4.11000  T+ 2
-25.84000 31.78000  4.36000  F   
-24.21000 31.57000  3.00000  F   
-23.10000 30.76000  3.03000  F   
-22.11000 30.70000  2.10000  T+ 2
-21.31000 30.10000  2.16000  F   
-22.36000 31.55000  1.13000  F   
-23.39000 32.39000  0.97000  T+ 2
-23.43000 32.98000  0.16000  F   
-24.35000 32.44000  1.89000  F   
-25.36000 33.32000  1.73000  T+ 3
-25.37000 33.92000  0.93000  F   
-26.09000 33.38000  2.41000  F   
#===== Standard bond lengths
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.610 0.020
1.480 0.020
1.480 0.020
1.610 0.020
1.430 0.020
1.530 0.020
1.435 0.020
1.530 0.020
1.435 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.480 0.020
1.330 0.020
1.330 0.020
1.330 0.020
1.000 0.020
1.330 0.020
1.390 0.020
1.390 0.020
1.400 0.020
1.000 0.020
1.400 0.020
1.400 0.020
1.000 0.020
1.400 0.020
1.330 0.020
1.000 0.020
1.000 0.020
#===== Standard bond angles
 1  2  3  109.6  2.0
 1  2  4  109.6  2.0
 3  2  4  120.0  2.0
 1  2  5  103.0  2.0
 3  2  5  109.6  2.0
 4  2  5  109.6  2.0
 2  5  6  120.0  2.0
 5  6  7  109.6  2.0
 5  6  8  109.6  2.0
 7  6  8  120.0  2.0
 5  6  9  103.0  2.0
 7  6  9  109.6  2.0
 8  6  9  109.6  2.0
 6  9 10  120.0  2.0
 9 10 11  109.6  2.0
 9 10 12  109.6  2.0
11 10 12  120.0  2.0
 9 10 13  103.0  2.0
11 10 13  109.6  2.0
12 10 13  109.6  2.0
10 13 14  120.0  2.0
13 14 15  111.0  2.0
14 15 16  109.5  2.0
14 15 17  109.5  2.0
16 15 17  104.0  2.0
15 16 23  104.0  2.0
15 17 18  109.5  2.0
15 17 20  104.0  2.0
18 17 20  109.5  2.0
17 20 21  109.5  2.0
17 20 23  104.0  2.0
21 20 23  109.5  2.0
16 23 20  104.0  2.0
16 23 24  109.5  2.0
20 23 24  109.5  2.0
23 24 25  126.0  2.0
23 24 29  126.0  2.0
25 24 29  108.0  2.0
24 25 26  108.0  2.0
25 26 28  108.0  2.0
26 28 29  108.0  2.0
26 28 35  132.0  2.0
29 28 35  120.0  2.0
24 29 28  108.0  2.0
24 29 30  132.0  2.0
28 29 30  120.0  2.0
29 30 32  120.0  2.0
30 32 33  120.0  2.0
32 33 35  120.0  2.0
28 35 33  120.0  2.0
28 35 36  120.0  2.0
33 35 36  120.0  2.0
17 18 19  109.5  2.0
20 21 22  109.5  2.0
25 26 27  126.0  2.0
27 26 28  126.0  2.0
29 30 31  120.0  2.0
31 30 32  120.0  2.0
32 33 34  120.0  2.0
34 33 35  120.0  2.0
35 36 37  120.0  2.0
35 36 38  120.0  2.0
37 36 38  120.0  2.0
   4   3   4   5
   2   6
   4   1   4   5
   4   1   3   5
   7   1   3   4   7   8   9
   3   2  10
   4   5   8   9
   4   5   7   9
   7   5   7   8  11  12  13
   3   6  14
   4   9  12  13
   4   9  11  13
   5   9  11  12  15
   4  10  16  17
   5  13  23  18  20
   5  14  17  20  24
   6  14  16  19  21  23
   3  15  20
   2  17
   6  15  18  22  16  24
   3  17  23
   2  20
   6  15  17  21  25  29
   6  16  20  26  28  30
   5  23  29  27  28
   4  24  29  35
   3  25  28
   7  25  27  24  30  33  36
   7  23  25  26  35  31  32
   4  24  28  33
   3  29  32
   5  29  31  34  35
   4  30  28  36
   3  32  35
   7  26  29  32  34  37  38
   3  28  33
   3  35  38
   3  35  37
-.360 0.990 -.635 -.635 -.360 0.000 0.000 0.000 -.360 0.990
-.635 -.635 -.360 0.000 0.180 -.360 0.180 0.000 0.000 0.150
-.548 0.398 0.000 0.000 0.000 -.280 0.280 0.000 0.000 -.280
0.280 0.000 -.280 0.280 0.000 -.830 0.415 0.415
*
NPO 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  11  11   2   2   1  11   0  10  14   1   1   1  13   6
#===== Atom names
C2  C3  C4  OH  HAA C5  C6  C1  N1  O3  O2  
#===== Bonds
   1   2   1   8   2   3   3   4   3   6   4   5   6   7   7   8   8   9   9  10
   9  11
#===== Rotatable bonds
   3   4   5
   8   9  10  11
#===== Torsion angles
   2   3   4   5   1   8   9  11
#===== Coordinates, hydrogen bonding and planarity
 10.09000 40.12000 61.39000  F   
 11.17000 39.26000 61.25000  F   
 11.03000 37.94000 61.67000  F   
 12.09000 37.10000 61.54000  T  4
 11.81000 36.17000 61.78000  F   
  9.85000 37.47000 62.22000  F   
  8.77000 38.35000 62.35000  F   
  8.89000 39.67000 61.94000  F   
  7.85000 40.50000 62.06000  F   
  6.78000 40.15000 62.54000  T-  
  7.90000 41.68000 61.70000  T-  
#===== Standard bond lengths
1.390 0.020
1.390 0.020
1.390 0.020
1.360 0.020
1.390 0.020
1.000 0.020
1.390 0.020
1.390 0.020
1.330 0.020
1.230 0.020
1.230 0.020
#===== Standard bond angles
 2  1  8  120.0  2.0
 1  2  3  120.0  2.0
 2  3  4  120.0  2.0
 2  3  6  120.0  2.0
 4  3  6  120.0  2.0
 3  6  7  120.0  2.0
 6  7  8  120.0  2.0
 1  8  7  120.0  2.0
 1  8  9  120.0  2.0
 7  8  9  120.0  2.0
 8  9 10  120.0  2.0
 8  9 11  120.0  2.0
10  9 11  120.0  2.0
 3  4  5  109.5  2.0
   4   3   7   9
   4   8   4   6
   4   1   5   7
   3   2   6
   2   3
   4   2   4   8
   4   3   1   9
   5   2   6  10  11
   3   1   7
   3   8  11
   3   8  10
0.000 0.000 0.150 -.548 0.398 0.000 0.000 0.000 0.720 -.360
-.360
*
CMN 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  19  19   6   6   1  19   0  12  27   1   1   1  17   6
#===== Atom names
HAF O4  C4  C3  O3  HAE C2  O2  HAD C1  N1  HAB HAC HAA 
C7  C5  C6  O6  HAG 
#===== Bonds
   1   2   2   3   3   4   3  16   4   5   4   7   5   6   7   8   7  10   8   9
  10  11  10  15  11  12  11  13  11  14  15  16  16  17  17  18  18  19
#===== Rotatable bonds
   2   3   1
   4   5   6
   7   8   9
  10  11  12  13  14
  16  17  18  19
  17  18  19
#===== Torsion angles
   1   2   3  16   3   4   5   6   4   7   8   9   7  10  11  14   3  16  17  18
  16  17  18  19
#===== Coordinates, hydrogen bonding and planarity
  3.54000-12.14000 -1.56000  F   
  2.73000-12.47000 -2.07000  T  4
  2.76000-11.93000 -3.40000  F   
  2.72000-13.09000 -4.43000  F   
  2.98000-12.60000 -5.76000  T  4
  3.97000-12.47000 -5.91000  F   
  1.38000-13.86000 -4.45000  F   
  1.46000-15.03000 -3.62000  T  4
  2.10000-15.69000 -4.03000  F   
  0.12000-13.06000 -4.02000  F   
 -0.11000-13.31000 -2.58000  T+ 5
  0.65000-12.90000 -2.07000  F   
 -0.09000-14.30000 -2.43000  F   
 -0.97000-12.91000 -2.28000  F   
  0.25000-11.55000 -4.26000  F   
  1.52000-11.04000 -3.55000  F   
  1.50000 -9.64000 -2.94000  F   
  0.72000 -9.68000 -1.74000  T  4
  0.69000 -8.76000 -1.34000  F   
#===== Standard bond lengths
1.000 0.020
1.430 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.470 0.020
1.530 0.020
1.000 0.020
1.000 0.020
1.000 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.000 0.020
#===== Standard bond angles
 2  3  4  109.5  2.0
 2  3 16  109.5  2.0
 4  3 16  109.5  2.0
 3  4  5  109.5  2.0
 3  4  7  111.0  2.0
 5  4  7  109.5  2.0
 4  7  8  109.5  2.0
 4  7 10  111.0  2.0
 8  7 10  109.5  2.0
 7 10 11  109.5  2.0
 7 10 15  111.0  2.0
11 10 15  109.5  2.0
10 15 16  111.0  2.0
 3 16 15  120.0  2.0
 3 16 17  120.0  2.0
15 16 17  120.0  2.0
16 17 18  111.0  2.0
 1  2  3  109.5  2.0
 4  5  6  109.5  2.0
 7  8  9  109.5  2.0
10 11 12  109.5  2.0
10 11 13  109.5  2.0
12 11 13  109.5  2.0
10 11 14  109.5  2.0
12 11 14  109.5  2.0
13 11 14  109.5  2.0
17 18 19  109.5  2.0
   2   3
   3   4  16
   6   1   5   7  15  17
   6   2  16   6   8  10
   3   3   7
   2   4
   6   3   5   9  11  15
   3   4  10
   2   7
   7   4   8  12  13  14  16
   3   7  15
   4  10  13  14
   4  10  12  14
   4  10  12  13
   5   7  11   3  17
   5   2   4  10  18
   4   3  15  19
   2  16
   2  17
0.398 -.548 0.150 0.150 -.548 0.398 0.150 -.548 0.398 0.127
0.129 0.248 0.248 0.248 0.000 0.000 0.150 -.548 0.398
*
G6D 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  11  11   2   2   1  11   0   9  14   1   1   1  12   6
#===== Atom names
HAB O3  C3  C4  C5  C6  O5  C1  C2  O2  HAA 
#===== Bonds
   1   2   2   3   3   4   3   9   4   5   5   6   5   7   7   8   8   9   9  10
  10  11
#===== Rotatable bonds
   2   3   1
   9  10  11
#===== Torsion angles
   1   2   3   9   3   9  10  11
#===== Coordinates, hydrogen bonding and planarity
 27.54049 25.25936 34.26093  F   
 27.73700 26.12100 33.79300  T  4
 28.79100 25.92100 32.84700  F   
 29.09700 27.27000 32.11600  F   
 30.23400 26.90000 31.10000  F   
 30.70900 28.04000 30.19300  F   
 31.35500 26.56900 31.92700  T-  
 31.26300 25.38700 32.64500  F   
 30.04500 25.40000 33.56400  F   
 29.71200 24.09500 34.02800  T  4
 30.49620 23.70094 34.50732  F   
#===== Standard bond lengths
1.000 0.020
1.430 0.020
1.530 0.020
1.530 0.020
1.530 0.020
1.530 0.020
1.435 0.020
1.435 0.020
1.530 0.020
1.430 0.020
1.000 0.020
#===== Standard bond angles
 2  3  4  109.5  2.0
 2  3  9  109.5  2.0
 4  3  9  111.0  2.0
 3  4  5  109.5  2.0
 4  5  6  109.5  2.0
 4  5  7  109.5  2.0
 6  5  7  109.5  2.0
 5  7  8  109.5  2.0
 7  8  9  109.5  2.0
 3  9  8  111.0  2.0
 3  9 10  109.5  2.0
 8  9 10  109.5  2.0
 1  2  3  109.5  2.0
 9 10 11  109.5  2.0
   2   3
   3   4   9
   5   1   5   8  10
   5   2   9   6   7
   3   3   8
   3   4   7
   4   4   6   9
   4   5   3  10
   5   2   4   7  11
   3   3   8
   2   9
0.398 -.548 0.150 0.000 0.000 0.180 -.360 0.180 0.150 -.548
0.398
*
CMN 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  19  19   6   6   1  19   0  12  27   1   1   1  17   6
#===== Atom names
HAF O4  C4  C3  O3  HAE C2  O2  HAD C1  N1  HAB HAC HAA 
C7  C5  C6  O6  HAG 
#===== Bonds
   1   2   2   3   3   4   3  16   4   5   4   7   5   6   7   8   7  10   8   9
  10  11  10  15  11  12  11  13  11  14  15  16  16  17  17  18  18  19
#===== Rotatable bonds
   2   3   1
   4   5   6
   7   8   9
  10  11  12  13  14
  16  17  18  19
  17  18  19
#===== Torsion angles
   1   2   3  16   3   4   5   6   4   7   8   9   7  10  11  14   3  16  17  18
  16  17  18  19
#===== Coordinates, hydrogen bonding and planarity
 38.52729 38.08137 39.88862  F   
 38.61800 38.32100 38.92200  T  4
 39.14100 37.19800 38.19700  F   
 38.24100 35.94400 38.39200  F   
 38.28400 35.51000 39.74700  T  4
 37.70197 34.70446 39.85813  F   
 38.67900 34.84500 37.43400  F   
 37.98400 33.65600 37.73600  T  4
 38.27065 32.93296 37.10746  F   
 38.58900 35.27700 35.96100  F   
 37.18900 35.17900 35.46800  T+ 5
 36.60947 35.81046 35.98318  F   
 37.15964 35.41636 34.49702  F   
 36.85595 34.24434 35.59240  F   
 39.08800 36.69100 35.77100  F   
 39.26900 37.57800 36.75200  F   
 39.58500 39.03800 36.49800  F   
 39.38300 39.41700 35.15900  T  4
 39.60571 40.38566 35.04902  F   
#===== Standard bond lengths
1.000 0.020
1.430 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.470 0.020
1.530 0.020
1.000 0.020
1.000 0.020
1.000 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.000 0.020
#===== Standard bond angles
 2  3  4  109.5  2.0
 2  3 16  109.5  2.0
 4  3 16  109.5  2.0
 3  4  5  109.5  2.0
 3  4  7  111.0  2.0
 5  4  7  109.5  2.0
 4  7  8  109.5  2.0
 4  7 10  111.0  2.0
 8  7 10  109.5  2.0
 7 10 11  109.5  2.0
 7 10 15  111.0  2.0
11 10 15  109.5  2.0
10 15 16  111.0  2.0
 3 16 15  120.0  2.0
 3 16 17  120.0  2.0
15 16 17  120.0  2.0
16 17 18  111.0  2.0
 1  2  3  109.5  2.0
 4  5  6  109.5  2.0
 7  8  9  109.5  2.0
10 11 12  109.5  2.0
10 11 13  109.5  2.0
12 11 13  109.5  2.0
10 11 14  109.5  2.0
12 11 14  109.5  2.0
13 11 14  109.5  2.0
17 18 19  109.5  2.0
   2   3
   3   4  16
   6   1   5   7  15  17
   6   2  16   6   8  10
   3   3   7
   2   4
   6   3   5   9  11  15
   3   4  10
   2   7
   7   4   8  12  13  14  16
   3   7  15
   4  10  13  14
   4  10  12  14
   4  10  12  13
   5   7  11   3  17
   5   2   4  10  18
   4   3  15  19
   2  16
   2  17
0.398 -.548 0.150 0.150 -.548 0.398 0.150 -.548 0.398 0.127
0.129 0.248 0.248 0.248 0.000 0.000 0.150 -.548 0.398
*
GLC 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  15  15   5   5   1  15   0  11  19   1   1   1  15   6
#===== Atom names
HAA O4  C4  C3  O3  HAC C2  O2  HAB C1  O5  C5  C6  O6  
HAD 
#===== Bonds
   1   2   2   3   3   4   3  12   4   5   4   7   5   6   7   8   7  10   8   9
  10  11  11  12  12  13  13  14  14  15
#===== Rotatable bonds
   2   3   1
   4   5   6
   7   8   9
  12  13  14  15
  13  14  15
#===== Torsion angles
   1   2   3  12   3   4   5   6   4   7   8   9   3  12  13  14  12  13  14  15
#===== Coordinates, hydrogen bonding and planarity
 30.74003 24.18137 32.69724  F   
 31.25100 24.24500 31.84000  T  4
 32.39500 23.35800 31.89200  F   
 31.91500 21.90400 32.10900  F   
 31.18100 21.83800 33.31600  T  4
 30.42478 22.49149 33.28297  F   
 33.09800 20.91600 32.10000  F   
 32.53400 19.64600 32.36300  T  4
 33.25670 18.95488 32.37051  F   
 33.84500 21.06100 30.73600  F   
 34.22500 22.41800 30.64300  T-  
 33.21200 23.41500 30.60700  F   
 33.92900 24.73900 30.37500  F   
 34.87100 24.99300 31.39700  T  4
 35.32666 25.86571 31.22161  F   
#===== Standard bond lengths
1.000 0.020
1.430 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.435 0.020
1.435 0.020
1.530 0.020
1.430 0.020
1.000 0.020
#===== Standard bond angles
 2  3  4  109.5  2.0
 2  3 12  109.5  2.0
 4  3 12  111.0  2.0
 3  4  5  109.5  2.0
 3  4  7  111.0  2.0
 5  4  7  109.5  2.0
 4  7  8  109.5  2.0
 4  7 10  111.0  2.0
 8  7 10  109.5  2.0
 7 10 11  109.5  2.0
10 11 12  109.5  2.0
 3 12 11  109.5  2.0
 3 12 13  109.5  2.0
11 12 13  109.5  2.0
12 13 14  109.5  2.0
 1  2  3  109.5  2.0
 4  5  6  109.5  2.0
 7  8  9  109.5  2.0
13 14 15  109.5  2.0
   2   3
   3   4  12
   6   1   5   7  11  13
   6   2  12   6   8  10
   3   3   7
   2   4
   5   3   5   9  11
   3   4  10
   2   7
   4   4   8  12
   4   7   3  13
   5   2   4  10  14
   4   3  11  15
   2  12
   2  13
0.398 -.548 0.150 0.150 -.548 0.398 0.150 -.548 0.398 0.180
-.360 0.180 0.150 -.548 0.398
*
FH7 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  24  25   7   7   1  25   0  18  34   1   1   1  27   6
#===== Atom names
 OAF HAA CAB CAC OAG HAB CAD OAH HAC CAE OAI CAA OAP CAO 
 CAQ NAR HAF CAJ CAK OAL HAD CAM OAN HAE 
#===== Bonds
   1   2   1   3   3   4   3  12   4   5   4   7   5   6   7   8   7  10   8   9
  10  11  11  12  12  13  13  14  14  15  14  22  15  16  16  17  16  18  18  19
  19  20  19  22  20  21  22  23  23  24
#===== Rotatable bonds
   1   3   2
   4   5   6
   7   8   9
  12  13  14  15  16  17  18  19  20  21  22  23  24
  13  14  15  16  17  18  19  20  21  22  23  24
  19  20  21
  22  23  24
#===== Torsion angles
   2   1   3  12   3   4   5   6   4   7   8   9   3  12  13  14  12  13  14  22
  18  19  20  21  14  22  23  24
#===== Coordinates, hydrogen bonding and planarity
 44.97950 42.12161 46.08649  T  4  
 45.92678 42.44153 46.10405  F     
 44.10150 43.22761 46.25849  F     
 42.64650 42.71361 46.22849  F     
 42.38750 41.94861 47.40049  T  4  
 43.01579 41.17152 47.43750  F     
 41.68850 43.89261 46.17849  F     
 40.34550 43.41661 46.09849  T  4  
 40.14549 42.84622 46.89515  F     
 42.03450 44.77361 44.97149  F     
 43.36850 45.26861 45.13749  T-    
 44.35750 44.23661 45.12249  F     
 45.66450 44.81361 45.31249  T-    
 45.93950 45.75761 44.28649  F     
 46.53850 45.01761 43.05750  F     
 46.90950 45.96661 42.02049  T+ 5  
 47.32335 45.45637 41.26658  F     
 47.79550 46.93761 42.41449  F     
 47.44850 47.66861 43.63549  F     
 48.39050 48.59761 44.03749  T  4  
 48.67225 49.14865 43.25203  F     
 46.96750 46.77761 44.81149  F     
 46.34850 47.56561 45.81749  T  4  
 47.00436 48.23375 46.16884  F     
#===== Standard bond lengths
1.000 0.020
1.430 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.435 0.020
1.435 0.020
1.435 0.020
1.435 0.020
1.530 0.020
1.530 0.020
1.470 0.020
1.000 0.020
1.470 0.020
1.530 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.000 0.020
#===== Standard bond angles
  1   3   4  109.5  2.0
  1   3  12  109.5  2.0
  4   3  12  111.0  2.0
  3   4   5  109.5  2.0
  3   4   7  111.0  2.0
  5   4   7  109.5  2.0
  4   7   8  109.5  2.0
  4   7  10  111.0  2.0
  8   7  10  109.5  2.0
  7  10  11  109.5  2.0
 10  11  12  109.5  2.0
  3  12  11  109.5  2.0
  3  12  13  109.5  2.0
 11  12  13  109.5  2.0
 12  13  14  109.5  2.0
 13  14  15  109.5  2.0
 13  14  22  109.5  2.0
 15  14  22  109.5  2.0
 14  15  16  109.5  2.0
 15  16  18  109.5  2.0
 16  18  19  109.5  2.0
 18  19  20  109.5  2.0
 18  19  22  111.0  2.0
 20  19  22  109.5  2.0
 14  22  19  111.0  2.0
 14  22  23  109.5  2.0
 19  22  23  109.5  2.0
  2   1   3  109.5  2.0
  4   5   6  109.5  2.0
  7   8   9  109.5  2.0
 15  16  17  109.5  2.0
 17  16  18  109.5  2.0
 19  20  21  109.5  2.0
 22  23  24  109.5  2.0
   3   4  12
   2   3
   6   2   5   7  11  13
   6   1  12   6   8  10
   3   3   7
   2   4
   5   3   5   9  11
   3   4  10
   2   7
   4   4   8  12
   4   7   3  13
   5   1   4  10  14
   5   3  11  15  22
   5  12  16  19  23
   5  13  22  17  18
   3  14  19
   3  15  18
   5  15  17  20  22
   5  16  21  14  23
   3  18  22
   2  19
   6  13  15  18  20  24
   3  14  19
   2  22
#===== Charges
-0.192  0.045  0.124  0.124 -0.192  0.045  0.124 -0.192  0.045  0.214
-0.183  0.210 -0.183  0.210  0.061 -0.275  0.000  0.061  0.124 -0.192
 0.045  0.124 -0.192  0.045
*
FH9 1
# AT  BND CHI TRS CRD #BL HGR HAT ANG 1-3 PAR PC HANG TYP
  28  29   9   9   1  29   0  20  40   1   1   1  30   6
#===== Atom names
 OAF HAA CAB CAC OAG HAB CAD OAH HAC CAE OAI CAA OAT CAM 
 CAN NAO HAG CAJ NAQ HAH OAP HAD CAK OAR HAE CAL OAS HAF 
#===== Bonds
   1   2   1   3   3   4   3  12   4   5   4   7   5   6   7   8   7  10   8   9
  10  11  11  12  12  13  13  14  14  15  14  26  15  16  16  17  16  18  18  19
  18  23  19  20  19  21  21  22  23  24  23  26  24  25  26  27  27  28
#===== Rotatable bonds
   1   3   2
   4   5   6
   7   8   9
  12  13   1   2   3   4   5   6   7   8   9  10  11
  13  14   1   2   3   4   5   6   7   8   9  10  11  12
  18  19  20  21  22
  19  21  22
  23  24  25
  26  27  28
#===== Torsion angles
   2   1   3  12   3   4   5   6   4   7   8   9   3  12  13  14  12  13  14  26
  16  18  19  21  18  19  21  22  18  23  24  25  14  26  27  28
#===== Coordinates, hydrogen bonding and planarity
 44.37726 41.94876 46.85027  T  4  
 45.32842 42.25643 46.87560  F     
 43.51326 43.06576 47.00527  F     
 42.05626 42.60077 46.96727  F     
 41.78826 41.79776 48.10727  T  4  
 42.40688 41.01211 48.11507  F     
 41.10826 43.80277 46.95727  F     
 39.77126 43.35376 46.77127  T  4  
 39.52515 42.72384 47.50790  F     
 41.49326 44.76876 45.83227  F     
 42.87226 45.17676 45.98027  T-    
 43.78726 44.07476 45.88827  F     
 45.12926 44.58177 46.04427  T-    
 45.53126 45.32177 44.89827  F     
 46.19526 44.32676 43.90127  F     
 46.69526 45.14476 42.75327  T+ 5  
 47.51938 44.67055 42.44353  F     
 47.04726 46.51277 42.75327  F     
 47.53726 47.02677 41.63227  T+ 5  
 48.43975 47.41211 41.82467  F     
 47.65926 46.11377 40.65427  T  4  
 48.03102 46.55002 39.83483  F     
 46.87126 47.36576 44.07927  F     
 47.97526 48.16577 44.36827  T  4  
 48.17160 48.75851 43.58717  F     
 46.53926 46.41477 45.30527  F     
 45.99126 47.16977 46.39527  T  4  
 46.63856 47.88203 46.66672  F     
#===== Standard bond lengths
1.000 0.020
1.430 0.020
1.530 0.020
1.530 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.435 0.020
1.435 0.020
1.435 0.020
1.435 0.020
1.530 0.020
1.530 0.020
1.470 0.020
1.000 0.020
1.470 0.020
1.470 0.020
1.530 0.020
1.000 0.020
1.460 0.020
1.000 0.020
1.430 0.020
1.530 0.020
1.000 0.020
1.430 0.020
1.000 0.020
#===== Standard bond angles
  1   3   4  109.5  2.0
  1   3  12  109.5  2.0
  4   3  12  111.0  2.0
  3   4   5  109.5  2.0
  3   4   7  111.0  2.0
  5   4   7  109.5  2.0
  4   7   8  109.5  2.0
  4   7  10  111.0  2.0
  8   7  10  109.5  2.0
  7  10  11  109.5  2.0
 10  11  12  109.5  2.0
  3  12  11  109.5  2.0
  3  12  13  109.5  2.0
 11  12  13  109.5  2.0
 12  13  14  109.5  2.0
 13  14  15  109.5  2.0
 13  14  26  109.5  2.0
 15  14  26  109.5  2.0
 14  15  16  109.5  2.0
 15  16  18  109.5  2.0
 16  18  19  109.5  2.0
 16  18  23  109.5  2.0
 19  18  23  109.5  2.0
 18  19  21  109.5  2.0
 18  23  24  109.5  2.0
 18  23  26  111.0  2.0
 24  23  26  109.5  2.0
 14  26  23  111.0  2.0
 14  26  27  109.5  2.0
 23  26  27  109.5  2.0
  2   1   3  109.5  2.0
  4   5   6  109.5  2.0
  7   8   9  109.5  2.0
 15  16  17  109.5  2.0
 17  16  18  109.5  2.0
 18  19  20  109.5  2.0
 20  19  21  109.5  2.0
 19  21  22  109.5  2.0
 23  24  25  109.5  2.0
 26  27  28  109.5  2.0
   3   4  12
   2   3
   6   2   5   7  11  13
   6   1  12   6   8  10
   3   3   7
   2   4
   5   3   5   9  11
   3   4  10
   2   7
   4   4   8  12
   4   7   3  13
   5   1   4  10  14
   5   3  11  15  26
   5  12  16  23  27
   5  13  26  17  18
   4  14  19  23
   3  15  18
   7  15  17  20  21  24  26
   4  16  23  22
   3  18  21
   3  18  20
   2  19
   6  16  19  25  14  27
   3  18  26
   2  23
   6  13  15  18  24  28
   3  14  23
   2  26
#===== Charges
-0.179  0.059  0.136  0.136 -0.179  0.058  0.136 -0.179  0.058  0.227
-0.170  0.223 -0.170  0.223  0.073 -0.262  0.013  0.136 -0.262  0.014
-0.179  0.058  0.136 -0.179  0.058  0.136 -0.179  0.058
*END

