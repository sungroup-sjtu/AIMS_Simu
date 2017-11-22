#!/usr/bin/env python3

import sys

sys.path.append('..')

import numpy as np
from app.models import *
from app.models_cv import *
from app.models_yaws import *

content = ''
for f in sys.argv[2:]:
    with open(f) as fin:
        content += fin.read()

molecules = []

for line in content.splitlines():
    line = line.strip()
    if line == '' or line.startswith('#'):
        continue

    mol = YawsMolecule.query.filter(YawsMolecule.cas == line.split()[2]).first()
    molecules.append(mol)

dev_density_list = []
dev_expansion_list = []
dev_hvap_list = []
dev_cp_list = []
dev_compressibility_list = []
n_sim = 0
n_exp = 0

print('#ID formula smiles density ref simu Cp ref simu Hvap ref simu expan ref simu')
for mol in molecules:
    task = Task.query.filter(Task.smiles_list == '["%s"]' % mol.smiles).first()
    if task is None:
        continue
    if task.post_result is None:
        continue

    n_sim += 1

    if sys.argv[1] == '298':
        T = 298
    elif sys.argv[1] == 'Tvap':
        T = mol.Tvap
        if T == None:
            continue
    elif sys.argv[1] == 'Tc':
        Tc = mol.Tc
        if Tc == None:
            continue
        else:
            T = Tc * 0.8

    Pvap = mol.get_Pvap(T)
    if Pvap is None:
        continue

    post_result = task.get_post_result(T, Pvap)
    if post_result['score-density'] < 0.999 or post_result['score-e_inter'] < 0.999:
        continue

    cv = Cv.query.filter(Cv.cas == mol.cas).first()
    if cv is None:
        cp_sim = None
    else:
        cv_intra = cv.get_post_result(T)
        cp_sim = cv_intra + post_result['Cp_inter'] + post_result['Cp_PV']

    density = mol.get_density(T)
    expansion = mol.get_expansion(T)
    hvap = mol.get_Hvap(T)
    cp = mol.get_Cp(T)
    # print(mol.id, mol.formula, mol.smiles,
    #       density, density, post_result['density'],
    #       cp, cp, cp_sim,
    #       hvap, hvap, post_result['Hvap'],
    #       expansion, expansion, post_result['expansion']
    #       )

    # density_sim = post_result['density']
    # hvap_sim = post_result['Hvap']
    # expansion_sim = post_result['expansion']
    #
    # if density != None:
    #     dev = (density_sim - density) / density * 100
    #     if density_sim > 0.1 and dev < 100:
    #         dev_density_list.append(dev)
    # if cp != None and cp_sim != None:
    #     dev = (cp_sim - cp) / cp * 100
    #     if cp_sim > 10 and dev < 100:
    #         dev_cp_list.append(dev)
    # if expansion != None:
    #     dev = (expansion_sim - expansion) / expansion * 100
    #     if expansion_sim > 0.0002 and dev < 100:
    #         dev_expansion_list.append(dev)
    # if hvap != None:
    #     dev = (hvap_sim - hvap) / hvap * 100
    #     if hvap_sim > 10 and dev < 100:
    #         dev_hvap_list.append(dev)

    T, compressibility = mol.get_compressibility_T()
    if T is None:
        continue
    Pvap = mol.get_Pvap(T)
    if Pvap is None:
        continue
    post_result = task.get_post_result(T, Pvap)
    if post_result['score-density'] < 0.999 or post_result['score-e_inter'] < 0.999:
        continue

    compressibility_sim = post_result['compressibility']
    print(mol.id, mol.formula, mol.smiles, compressibility, compressibility,
          compressibility_sim, mol.compressibility_is_exp)
    if compressibility != None and compressibility_sim > 0 and mol.compressibility_is_exp:
        dev = (compressibility_sim - compressibility) / compressibility * 100
        if compressibility_sim > 0 and dev < 100:
            dev_compressibility_list.append(dev)

    n_exp += 1

print(n_sim, n_exp)

# print('%.1f/%.1f/%.1f\t%.1f/%.1f/%.1f\t%.1f/%.1f/%.1f\t%.1f/%.1f/%.1f' % (
#     np.mean(dev_density_list), np.mean(np.abs(dev_density_list)), np.std(dev_density_list),
#     np.mean(dev_cp_list), np.mean(np.abs(dev_cp_list)), np.std(dev_cp_list),
#     np.mean(dev_hvap_list), np.mean(np.abs(dev_hvap_list)), np.std(dev_hvap_list),
#     np.mean(dev_expansion_list), np.mean(np.abs(dev_expansion_list)), np.std(dev_expansion_list)
# ))
print(np.mean(dev_compressibility_list), np.std(dev_compressibility_list))
