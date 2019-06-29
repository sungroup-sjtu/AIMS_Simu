#!/usr/bin/env python3
# coding=utf-8

import os, sys, time
import pybel

sys.path.append('..')
from app import create_app, db
from app.models import Task, Compute, PbsJob

app = create_app(sys.argv[1])
app.app_context().push()

smarts_bad = {
    # Not covered by TEAM FF
    'radicalC'                     : '[#6;v0,v1,v2,v3]',
    '*=*=*'                        : '*=*=*',
    '*#*~*#*'                      : '*#*~*#*',
    '[F,Cl,Br]~[!#6]'              : '[F,Cl,Br]~[!#6]',
    '*#*~[!#6]'                    : '*#*~[!#6]',
    '[NX2,NX4]'                    : '[NX2,NX4]',
    'O~N(~[!$([OX1])])~[!$([OX1])]': 'O~N(~[!$([OX1])])~[!$([OX1])]',
    'peroxide'                     : 'O~O',
    'N~N'                          : 'N~N',
    '[O,N]*[O,N;H1,H2]'            : '[O,N]*[O,N;H1,H2]',
    'C=C~[O,N;H1,H2]'              : 'C=C~[O,N;H1,H2]',
    'beta-dicarbonyl'              : 'O=C~*~C=O',
    'a=*'                          : 'a=*',
    'o'                            : 'o',
    '[n;r5]'                       : '[n;r5]',
    'pyridine-N-oxide'             : '[nX3;r6]',
    'triazine(zole)'               : '[$(nnn),$(nnan),$(nanan)]',
    '[R3]'                         : '[R3]',
    '[r3,r4;R2]'                   : '[r3,r4;R2]',
    '[r3,r4;#6X3]'                 : '[r3,r4;#6X3]',
    '[r3,r4]~[!#6]'                : '[r3,r4]~[!#6]',
    'nitrate'                      : 'O[NX3](~[OX1])~[OX1]',
    'amide'                        : 'O=C[NX3]',
    'acyl-halide'                  : 'O=C[F,Cl,Br]',
    'polybenzene'                  : 'c1ccc2c(c1)cccc2',
    # Covered by TEAM FF but the results are not good
    '[r5;#6X3]'                    : '[r5;#6X3]',
    '[r5]~[!#6]'                   : '[r5]~[!#6]',
    'cyclo-ester'                  : '[C;R](=O)O',
    'C=C~[O,N;H0]'                 : 'C=C~[O,N;H0]',
    'C=C-X'                        : 'C=C[F,Cl,Br]',
    '[F,Cl,Br][#6][F,Cl,Br]'       : '[F,Cl,Br][#6][F,Cl,Br]',
    'alkyne'                       : '[CX2]#[CX2]',
    'acid'                         : 'C(=O)[OH]',
    'nitrile'                      : '[NX1]#[CX2][C,c]',
    'nitro'                        : '[C,c][NX3](~[OX1])~[OX1]',
    'N-aromatic'                   : 'n',
    'halogen'                      : '[F,Cl,Br]',
}

smarts_list = [pybel.Smarts(smarts) for smarts in smarts_bad.values()]


def main():
    tasks = Task.query
    tasks.update({'remark': None})
    n_total = tasks.count()
    for i, task in enumerate(tasks):
        m = task.get_mol_list()[0]
        for s in smarts_list:
            if s.findall(m) != []:
                print(f'{i} / {n_total}', task)
                task.remark = 'bad'
                break
    db.session.commit()


if __name__ == '__main__':
    main()
