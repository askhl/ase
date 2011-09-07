# atomization energies in kcal / mol (= 43.364 meV).
# All values evaluated with self-consistent orbitals and densities.
# Geometry optimization to within 10meV/Ang.
# Data from reference [3].
atomization_vasp = {
# Molecule        Expt. PBE_VASP PBE_G03 PBE0_VASP PBE0_G03
'LiH'       : (     58,   53.5,   53.5,   52.6,   52.9,),
'BeH'       : (     48,   55.5,   55.6,   55.4,   56.0,),
'CH'        : (     84,   84.7,   84.8,   83.3,   83.0,),
'CH2_s3B1d' : (    189,  194.4,  194.6,  193.9,  193.8,),
'CH2_s1A1d' : (    182,  178.8,  179.1,  176.3,  176.5,),
'CH3'       : (    306,  309.7,  310.1,  308.3,  308.6,),
'CH4'       : (    420,  419.6,  420.2,  417.2,  417.9,),
'NH'        : (     82,   88.6,   88.6,   85.2,   85.3,),
'NH2'       : (    182,  188.7,  188.9,  183.4,  183.3,),
'NH3'       : (    297,  301.7,  302.3,  294.7,  295.3,),
'OH'        : (    107,  109.7,  110.1,  105.4,  105.8,),
'H2O'       : (    233,  233.7,  234.5,  226.4,  227.3,),
'HF'        : (    142,  141.5,  142.2,  136.2,  137.0,),
'SiH2_s1A1d': (    154,  147.9,  148.0,  147.2,  147.4,),
'SiH2_s3B1d': (    131,  131.3,  131.8,  132.2,  132.5,),
'SiH3'      : (    226,  222.2,  222.6,  223.5,  223.8,),
'SiH4'      : (    324,  313.3,  313.7,  315.0,  315.7,),
'PH2'       : (    153,  154.5,  154.6,  153.0,  153.2,),
'PH3'       : (    241,  239.0,  239.3,  237.3,  237.5,),
'SH2'       : (    182,  182.0,  182.2,  179.8,  180.0,),
'HCl'       : (    107,  106.3,  106.5,  104.4,  105.0,),
'Li2'       : (     26,   19.9,   20.1,   19.3,   19.3,),
'LiF'       : (    139,  138.4,  139.0,  131.0,  131.9,),
'C2H2'      : (    404,  414.5,  415.1,  404.5,  404.7,),
'C2H4'      : (    562,  571.0,  571.9,  563.8,  564.2,),
'C2H6'      : (    711,  716.0,  717.1,  711.6,  712.4,),
'CN'        : (    179,  197.5,  197.7,  179.1,  179.1,),
'HCN'       : (    313,  326.3,  326.5,  311.1,  311.5,),
'CO'        : (    261,  268.6,  269.1,  255.3,  255.8,),
'HCO'       : (    279,  294.9,  295.5,  280.5,  280.9,),
'H2CO'      : (    376,  385.5,  386.3,  371.9,  372.8,),
'CH3OH'     : (    513,  519.3,  520.4,  509.0,  510.3,),
'N2'        : (    227,  243.7,  243.9,  225.3,  225.9,),
'N2H4'      : (    437,  452.7,  453.7,  437.9,  438.8,),
'NO'        : (    153,  172.0,  172.5,  153.3,  153.8,),
'O2'        : (    118,  143.3,  144.0,  124.1,  124.9,),
'H2O2'      : (    268,  281.6,  282.6,  262.7,  263.8,),
'F2'        : (     38,   52.6,   53.0,   35.2,   35.3,),
'CO2'       : (    392,  415.4,  416.5,  390.8,  392.0,),
'Na2'       : (     19,   17.7,   18.1,   15.6,   15.9,),
'Si2'       : (     74,   81.3,   81.4,   76.5,   77.3,),
'P2'        : (    116,  121.5,  121.7,  111.8,  111.7,),
'S2'        : (     98,  115.4,  115.2,  107.3,  107.0,),
'Cl2'       : (     57,   65.8,   65.8,   60.1,   59.9,),
'NaCl'      : (     99,   93.6,   94.5,   92.1,   93.6,),
'SiO'       : (    191,  195.6,  196.6,  182.2,  183.3,),
'CS'        : (    172,  179.5,  179.6,  168.0,  168.2,),
'SO'        : (    122,  141.5,  141.3,  127.9,  127.3,),
'ClO'       : (     62,   81.6,   81.5,   67.4,   67.6,),
'ClF'       : (     62,   72.3,   72.5,   61.3,   61.3,),
'Si2H6'     : (    533,  519.5,  520.4,  522.2,  523.3,),
'CH3Cl'     : (    395,  399.4,  400.2,  395.0,  395.7,),
'CH3SH'     : (    473,  477.8,  478.6,  472.7,  473.5,),
'HOCl'      : (    165,  175.2,  175.7,  162.9,  163.3,),
'SO2'       : (    253,  281.1,  280.7,  254.1,  253.5,),
}

# Experimental, and calculated bindinglengths of 16 diatomic molecules
# of the G2-1 test set. In Angstroms.
# Data from reference [3].
diatomic = {
#System Expt. PBEVASP PBEG03 PBE0VASP PBE0G03
'BeH': (1.343, 1.354, 1.353, 1.350, 1.348),
'CH' : (1.120, 1.136, 1.136, 1.124, 1.124),
'Cl2': (1.988, 1.999, 2.004, 1.973, 1.978),
'ClF': (1.628, 1.648, 1.650, 1.614, 1.617),
'ClO': (1.570, 1.576, 1.577, 1.554, 1.555),
'CN' : (1.172, 1.173, 1.174, 1.159, 1.159),
'CO' : (1.128, 1.136, 1.135, 1.122, 1.122),
'F2' : (1.412, 1.414, 1.413, 1.377, 1.376),
'FH' : (0.917, 0.932, 0.930, 0.919, 0.918),
'HCl': (1.275, 1.287, 1.288, 1.276, 1.278),
'Li2': (2.673, 2.728, 2.728, 2.727, 2.727),
'LiF': (1.564, 1.583, 1.575, 1.571, 1.561),
'LiH': (1.595, 1.604, 1.604, 1.602, 1.597),
'N2' : (1.098, 1.103, 1.102, 1.089, 1.089),
'O2' : (1.208, 1.218, 1.218, 1.193, 1.192),
'Na2': (3.079, 3.087, 3.076, 3.086, 3.086),
}

def convert(input, column):
    keys = input.keys()
    keys.sort()
    data = {}
    for k in keys:
        data[k] = input[k][column]
    return data

info = {}

info['atomization energy'] = {}

info['atomization energy'].update({'reference': convert(atomization_vasp, 0)})
info['atomization energy'].update({'PBE': convert(atomization_vasp, 1)})
info['atomization energy'].update({'PBE0': convert(atomization_vasp, 3)})

info['bondlength'] = {}

info['bondlength'].update({'reference': convert(diatomic, 0)})
info['bondlength'].update({'PBE': convert(diatomic, 1)})
info['bondlength'].update({'PBE0': convert(diatomic, 3)})

## Note the difference between the experimental data from Blaha [1] and
## from Kresse [3]:
## Mol  BLAHA KRESSE
## LiH :  57.8  58.0
## CH4 : 419.3 420.0
## NH3 : 297.4 297.0
## OH  : 106.4 107.0
## H2O : 232.2 233.0
## HF  : 140.8 142.0
## Li2 :  24.4  26.0
## LiF : 138.9 139.0
## C2H2: 405.4 404.0
## C2H4: 562.6 562.0
## HCN : 311.9 313.0
## CO  : 259.3 261.0
## N2  : 228.5 227.0
## NO  : 152.9 153.0
## O2  : 120.5 118.0
## F2  :  38.5  38.0
## P2  : 117.3 116.0
## Cl2 :  58.0  57.0
## -----------------
## MAE :   0.0   1.0

## and the difference between the PBE results of Blaha [1] and those from
## VASP and GAUSSIAN-03 [3]
## Mol  BLAHA  VASP   G03
## LiH :  53.5  53.5  53.5
## CH4 : 419.8 419.6 420.2
## NH3 : 301.7 301.7 302.3
## OH  : 109.8 109.7 110.1
## H2O : 234.2 233.7 234.5
## HF  : 142.0 141.5 142.2
## Li2 :  19.9  19.9  20.1
## LiF : 138.6 138.4 139.0
## C2H2: 414.9 414.5 415.1
## C2H4: 571.5 571.0 571.9
## HCN : 326.1 326.3 326.5
## CO  : 268.8 268.6 269.1
## N2  : 243.2 243.7 243.9
## NO  : 171.9 172.0 172.5
## O2  : 143.7 143.3 144.0
## F2  :  53.4  52.6  53.0
## P2  : 121.1 121.5 121.7
## Cl2 :  65.1  65.8  65.8
## -----------------------
## MAE :   0.0   0.3   0.4
##
## Where in the last two, geometry optimization has been performed, but not in
## the first.

# References:
# [1]:
# Kurth, Perdew, and Blaha
# Molecular and Solid-State Tests of Density Functional Approximations
# International Journal of Quantum Chemistry, Vol. 85, 889-909 (1999)
# [2]:
# Krieger, Chen, Iafrate, and Savin
# In Electron Correlations and Materials Properties
# Gonis and Kioussis; eds.
# Plenum: New York, 1999.
# [3]:
# Paier, Hirschl, Marsman, and Kresse
# The Perdew-Burke-Ernzerhof exchange-correlation functional applied to the
# G2-1 test set using a plane wave basis set
# The Journal of Chemical Physics, Vol 122, 234102 (2005)
