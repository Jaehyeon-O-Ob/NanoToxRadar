# -*- coding: utf-8 -*-
"""
Created on Mon May 18 14:23:58 2020

@author: Gil

*Reference: http://www.wiredchemist.com/chemistry/data/metallic-radii
A.F. Wells, "Structural Inorganic Chemistry," 5th ed., Clarendon Press, Oxford,
1984, p. 1288 (metallic radii for 12-coordination);
Huheey, pp. 292 (covalent radii for nonmetals);
R.D. Shannon, Acta Crystallogr., Sect. A: Found. Crystallogr., 32, 751 (1976) (ionic radii for 6-coordination).

Unit of radius is pm.
*Metallic radii for 12-coordination are given for all metals.
*Ionic radii are for six-coordination.

Metallic Radius Reference (Empirical size) - https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
*Ta - https://en.wikipedia.org/wiki/Tantalum
*Gd - https://en.wikipedia.org/wiki/Gadolinium
*Mo - https://en.wikipedia.org/wiki/Molybdenum
*Pr - https://en.wikipedia.org/wiki/Praseodymium
*Sm - https://en.wikipedia.org/wiki/Samarium
*B - https://en.wikipedia.org/wiki/Boron
"""

metallic_radii = {'Ag': 144,'Al': 143, 'Au':144, 'B': 85, 'Ba':224, 'Be':112,
                  'Bi':182, 'Br':114, 'Ca':197, 'Cd':152, 'Ce':182,
                  'Co':125, 'Cr':129, 'Cs':272, 'Cu':128, 'Fe':126,
                  'Ga':153, 'Ge':139, 'Hf':159, 'Hg': 155, 'In': 167,
                  'Ir':136, 'K':235, 'La':188, 'Li':157, 'Lu':172,
                  'Mg':160, 'Mn':137, 'Na':191, 'Ni':125, 'Os':135,
                  'Pb':175, 'Pd':137, 'Po':153, 'Pt':139, 'Rb':250, 'Rh':134,
                  'Ru':134, 'Sb':161, 'Sc':164, 'Sn':158, 'Sr':215,
                  'Th':180, 'Ti':147, 'Tl':171, 'U':156, 'V':135,
                  'W':141, 'Y':182, 'Zn':137, 'Zr':160, 'Ta':146, 'Gd':180,
                 'Mo':145, 'Pr':185, 'Sm':185}

ionic_radii = {'Ag+': 115, 'Al+3':53, 'As+3':58, 'Au+':137, 'Ba+2':135, 'Be+2':45,
               'Bi+3':103, 'Br-':196, 'Ca+2':100, 'Cd+2':95, 'Ce+3':102, 'Cl-':184,
               'Co+2':70, 'Co+3':60, 'Cr+3':62, 'Cs+':167, 'Cu+':77, 'Cu+2':73,
               'F-':133, 'Fe+2':70, 'Fe+3':60, 'Ga+3':62, 'Hg+2':102, 'I-':220,
               'In+3':80, 'Ir+3':82, 'K+':138, 'La+3':103, 'Li+':76, 'Lu+3':86,
               'Mg+2':72, 'Mn+2':70, 'Na+':102, 'Ni+2':70, 'O-2':140, 'Pb+2':119, 'Pd+2':86,
               'Pt+2':80, 'Pt+4':77, 'Rb+':152, 'S-2':184, 'Sb+3':76, 'Sc+3':75,
               'Se-2':198, 'Sn+2':118, 'Sr+2':118, 'Te-2':221, 'Ti+3':67, 'Tl+':150,
               'Tl+3':88, 'V+2':79, 'Y+3':90, 'Zn+2':74}

"""
*Reference: https://en.wikipedia.org/wiki/Ionic_radius#cite_note-Identification-7
Shannon radii
Crystal ionic radii is ionic radii within ionic crystal system. (O-2: 126 pm)
Effective ionic radii is ionic radii consistent with Pauling's radii. (O-2: 140 pm)
"""
crystal_ionic_radii = {'H+1':-4,'Li+1':90,'Be+2':59, 'B+3':41, 'C+4':30, 'N-3':132,
                       'N+3':30, 'N+5':27,'O-2':126,'F-1':119,'F+7':22, 'Na+1':116,
                       'Mg+2':86,'Al+3':67.5,'Si+4':54,'P+3':58,'P+5':52,'S-2':170,
                       'S+4':51, 'S+6':43,'Cl-1':167,' Cl+5':26,' Cl+7':41,'K+1':152,
                       'Ca+2':114, 'Sc+3':88.5, 'Ti+2':100, 'Ti+3':81, 'Ti+4':74.5,
                       'V+2':93, 'V+3':78, 'V+4':72, 'V+5':68, 'Cr+2':87, 'Cr+3':75.5 ,
                       'Cr+4':69, 'Cr+5':63,'Cr+6':58,'Mn+2':81, 'Mn+3':72, 'Mn+4':67,
                       'Mn+5':47, 'Mn+6':39.5, 'Mn+7':60, 'Fe+2':75, 'Fe+3':69, 'Fe+4':72.5,
                       'Fe+6':39,'Co+2':79, 'Co+3':68.5, 'Ni+2':83,' Ni+3':70, 'Ni+4':62,
                       'Cu+1':91, 'Cu+2':87, 'Cu+3':68,'Zn+2':88, 'Ga+3':76, 'Ge+2':87,
                       'Ge+4':67, 'As+3':72, 'As+5':60, 'Se-2':184, 'Se+4':64, 'Se+6':56,
                       'Br-1':182, 'Br+3':73, 'Br+5':45, 'Br+7':53, 'Rb+1':166, 'Sr+2':132,
                       'Y+3':104, 'Zr+4':86, 'Nb+3':86, 'Nb+4':82, 'Nb+5':78, 'Mo+3':83,
                       'Mo+4':79, 'Mo+5':75, 'Mo+6':73, 'Tc+4':78.5, 'Tc+5':74, 'Tc+7':70,
                       'Ru+3':82, 'Ru+4':76, 'Ru+5':70.5, 'Ru+7':52, 'Ru+8':50,
                       'Rh+3':80.5, 'Rh+4':74, 'Rh+5':69, 'Pd+1':73, 'Pd+2':100, 'Pd+3':90, 'Pd+4':75.5,
                       'Ag+1':129, 'Ag+2':108, 'Ag+3':89, 'Cd+2':109, 'In+3':94, 'Sn+4':83,
                       'Sb+3':90, 'Sb+5':74, 'Te-2':207, 'Te+4':111, 'Te+6':70, 'I-1':206, 'I+5':109,
                       'I+7':67, 'Xe+8':62, 'Cs+1':181, 'Ba+2':149, 'La+3':117.2, 'Ce+3':115, 'Ce+4':101,
                       'Pr+3':113, 'Pr+4':99, 'Nd+2':143, 'Nd+3':112.3, 'Pm+3':111,
                       'Sm+2':136, 'Sm+3':109.8, 'Eu+2':131, 'Eu+3':108.7, 'Gd+3':107.8,
                       'Tb+3':106.3, 'Tb+4':90, 'Dy+2':121, 'Dy+3':105.2, 'Ho+3':104.1,
                       'Er+3':103, 'Tm+2':117, 'Tm+3':102, 'Yb+2':116, 'Yb+3':100.8,
                       'Lu+3':100.1, 'Hf+4':85, 'Ta+3':86, 'Ta+4':82, 'Ta+5':78,
                       'W+4':80, 'W+5':76, 'W+6':74, 'Re+4':77, 'Re+5':72, 'Re+6':69, 'Re+7':67,
                       'Os+4':77, 'Os+5':71.5, 'Os+6':68.5, 'Os+7':66.5, 'Os+8':53,
                       'Ir+3':82, 'Ir+4':76.5, 'Ir+5':71, 'Pt+2':94, 'Pt+4':76.5, 'Pt+5':71,
                       'Au+1':151, 'Au+3':99, 'Au+5':71, 'Hg+1':133, 'Hg+2':116,
                       'Tl+1':164, 'Tl+3':102.5, 'Pb+2':133, 'Pb+4':91.5,
                       'Bi+3':117, 'Bi+5':90, 'Po+4':108, 'Po+6':81, 'At+7':76,
                       'Fr+1':194, 'Ra+2':162,' Ac+3':126, 'Th+4':108, 'Pa+3':116,
                       'Pa+4':104, 'Pa+5':92, 'U+3':116.5, 'U+4':103, 'U+5':90, 'U+6':87,
                       'Np+2':124, 'Np+3':115, 'Np+4':101, 'Np+5':89, 'Np+6':86, 'Np+7':85,
                       'Pu+3':114, 'Pu+4':100, 'Pu+5':88, 'Pu+6':85, 'Am+2':140, 'Am+3':111.5,
                       'Am+4':99, 'Cm+3':111, 'Cm+4':99, 'Bk+3':110, 'Bk+4':97,
                       'Cf+3':109, 'Cf+4':96.1, 'Es+3':92.8}

effective_ionic_radii = {'H+1':-18,'Li+1':76,'Be+2':45, 'B+3':27, 'C+4':16, 'N-3':146,
                         'N+3':16, 'N+5':13,'O-2':140,'F-1':133,'F+7':8, 'Na+1':102,
                         'Mg+2':72,'Al+3':53.5,'Si+4':40,'P+3':44,'P+5':38,'S-2':184,
                         'S+4':37, 'S+6':29,'Cl-1':181,'Cl+5':12,'Cl+7':27,'K+1':138,
                         'Ca+2':100, 'Sc+3':74.5, 'Ti+2':86, 'Ti+3':67, 'Ti+4':60.5,
                         'V+2':79, 'V+3':64, 'V+4':58, 'V+5':54, 'Cr+2':73, 'Cr+3':61.5,
                         'Cr+4':55, 'Cr+5':49,'Cr+6':44,'Mn+2':67, 'Mn+3':58, 'Mn+4':53,
                         'Mn+5':33, 'Mn+6':25.5, 'Mn+7':46, 'Fe+2':61, 'Fe+3':55, 'Fe+4':58.5,
                         'Fe+6':25,'Co+2':65, 'Co+3':54.5, 'Ni+2':69,'Ni+3':56, 'Ni+4':48,
                         'Cu+1':77, 'Cu+2':73, 'Cu+3':54,'Zn+2':74, 'Ga+3':62, 'Ge+2':73,
                         'Ge+4':53, 'As+3':58, 'As+5':46, 'Se-2':198, 'Se+4':50, 'Se+6':42,
                         'Br-1':196, 'Br+3':59, 'Br+5':31, 'Br+7':39, 'Rb+1':152, 'Sr+2':118,
                         'Y+3':90, 'Zr+4':72, 'Nb+3':72, 'Nb+4':68, 'Nb+5':64, 'Mo+3':69,
                         'Mo+4':65, 'Mo+5':61, 'Mo+6':59, 'Tc+4':64.5, 'Tc+5':60, 'Tc+7':56,
                         'Ru+3':68, 'Ru+4':62, 'Ru+5':56.5, 'Ru+7':38, 'Ru+8':36,
                         'Rh+3':66.5, 'Rh+4':60, 'Rh+5':55, 'Pd+1':59, 'Pd+2':86, 'Pd+3':76, 'Pd+4':61.5,
                         'Ag+1':115, 'Ag+2':94, 'Ag+3':75, 'Cd+2':95, 'In+3':80, 'Sn+4':69,
                         'Sb+3':76, 'Sb+5':60, 'Te-2':221, 'Te+4':97, 'Te+6':56, 'I-1':220, 'I+5':95,
                         'I+7':53, 'Xe+8':48, 'Cs+1':167, 'Ba+2':135, 'La+3':103.2, 'Ce+3':101, 'Ce+4':87,
                         'Pr+3':99, 'Pr+4':85, 'Nd+2':129, 'Nd+3':98.3, 'Pm+3':97,
                         'Sm+2':122, 'Sm+3':95.8, 'Eu+2':117, 'Eu+3':94.7, 'Gd+3':93.5,
                         'Tb+3':92.3, 'Tb+4':76, 'Dy+2':107, 'Dy+3':91.2, 'Ho+3':90.1,
                         'Er+3':89, 'Tm+2':103, 'Tm+3':88, 'Yb+2':102, 'Yb+3':86.8,
                         'Lu+3':86.1, 'Hf+4':71, 'Ta+3':72, 'Ta+4':68, 'Ta+5':64,
                         'W+4':66, 'W+5':62, 'W+6':60, 'Re+4':63, 'Re+5':58, 'Re+6':55, 'Re+7':53,
                         'Os+4':63, 'Os+5':57.5, 'Os+6':54.5, 'Os+7':52.5, 'Os+8':39,
                         'Ir+3':68, 'Ir+4':62.5, 'Ir+5':57, 'Pt+2':80, 'Pt+4':62.5, 'Pt+5':57,
                         'Au+1':137, 'Au+3':85, 'Au+5':57, 'Hg+1':119, 'Hg+2':102,
                         'Tl+1':150, 'Tl+3':88.5, 'Pb+2':119, 'Pb+4':77.5,
                         'Bi+3':103, 'Bi+5':76, 'Po+4':94, 'Po+6':67, 'At+7':62,
                         'Fr+1':180, 'Ra+2':148,' Ac+3':112, 'Th+4':94, 'Pa+3':104,
                         'Pa+4':90, 'Pa+5':78, 'U+3':102.5, 'U+4':89, 'U+5':76, 'U+6':73,
                         'Np+2':110, 'Np+3':101, 'Np+4':87, 'Np+5':75, 'Np+6':72, 'Np+7':71,
                         'Pu+3':100, 'Pu+4':86, 'Pu+5':74, 'Pu+6':71, 'Am+2':126, 'Am+3':97.5,
                         'Am+4':85, 'Cm+3':97, 'Cm+4':85, 'Bk+3':96, 'Bk+4':83,
                         'Cf+3':95, 'Cf+4':82.1, 'Es+3':83.5}

"""
RDKit stores radii of atoms (unit: Angstrom, 1pm == 1Ang)
from rdkit import Chem
periodic_table = Chem.GetPeriodicTable()
atoms = [periodic_table.GetElementSymbol(i) for i in range(1,121)]
vdw = {a:periodic_table.GetRvdw(a) for a in atoms}
covalent = {a:periodic_table.GetRcovalent(a) for a in atoms}

Volume of particular conformer is calculated in RDKit.
from rdkit.Chem import AllChem
mol = Chem.AddHs(Chem.MolFromSmiles('C'))
AllChem.EmbedMolecule(mol)
AllChem.ComputeMolVolume(mol)

*Ref: https://en.wikipedia.org/wiki/Carbon%E2%80%93carbon_bond
Carbon-carbon single bond length: 154 pm

*Ref: https://en.wikipedia.org/wiki/Carbon%E2%80%93hydrogen_bond
Carbon-hydrogen single bond length: 109 pm

*Ref: https://en.wikipedia.org/wiki/Benzene
Carbon-carbon aromtatic bond length (Benzene): 139 pm
Carbon-hydrogen aromtatic bond length (Benzene): 109 pm

*Ref: https://en.wikipedia.org/wiki/Carbon%E2%80%93nitrogen_bond
Carbon-nitrogen (simple amine): 147.5 pm

*Ref: https://en.wikipedia.org/wiki/Carbon%E2%80%93oxygen_bond
Carbon-Oxygen with H: 143 pm (OH1R)
Carbon-Oxygen (-1): 136 pm (Ocarboxyl)
Carbon=Oxygen (double bond): 123 pm (O2C)

*Ref: https://en.wikipedia.org/wiki/Organosulfur_compounds
Carbon-Surfur: 183pm

*Ref: https://en.wikipedia.org/wiki/Silicon_dioxide
Silicon-Oxygen-Silicon: 164.6 pm
"""

covalent_radii = {'C1.5C':69.5, 'C':77, 'H1.5C':39.5, 'H':32,'N':70.5,
                  'OH1R':66, 'O2C':46, 'Ocarboxyl':59, 'S':106, 'Si':98.6}


"""
Volume of coating materials is needed covalent_radii and neutral radii.
Using covalent and neutral radii is more reasonable in this study to simplify the structure of particles.\
*Ref
F - Radius of Fluorine - https://en.wikipedia.org/wiki/Covalent_radius_of_fluorine#:~:text=The%20covalent%20radius%20of%20fluorine,radius%20is%20difficult%20to%20evaluate.
"""
neutral_radii = {'C':77, 'H':32, 'N':70.5, 'S': 106, 'Si':98.5, 'Zn':137,
                 'Na':191, 'Ga':153, 'Br':114, 'P':100, 'O':60, 'Cl':100, 'Ag': 144,
                 'Pd':137,'Cu':128, 'Pt':139, 'Gd':180, 'Au':144, 'B': 85, 'F':71}