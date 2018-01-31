from flask import request
from ..models import *
from ..models_cv import *
from ..models_yaws import *


class StatAction():
    def update_mol_task_list(self, smiles_list):
        self.yaws_list = []
        self.task_list = []

        for smiles in smiles_list:
            yaws: YawsMolecule = YawsMolecule.query.filter(YawsMolecule.smiles == smiles).first()
            if yaws == None:
                continue
            task = Task.query.filter(Task.smiles_list == '["%s"]' % smiles).first()
            if task == None or task.post_result == None:
                continue

            self.yaws_list.append(yaws)
            self.task_list.append(task)

    def get_T_P_from_T(self, yaws, T):
        if T == 'Tvap':
            T = yaws.get_Tvap()
        elif T == 'Tc':
            T = yaws.get_Tc()
            if T != None:
                T *= 0.8
                T = min(T, 600)
        if T == None:
            return None, None

        Pvap = yaws.get_Pvap(T)
        return T, Pvap

    def get_density(self, smiles_list, T=298):
        dens_list = []
        for i, yaws in enumerate(self.yaws_list):
            task = self.task_list[i]
            _T, _P = self.get_T_P_from_T(yaws, T)

            if None in (_T, _P):
                continue

            dens_exp = yaws.get_density(_T)
            dens_sim = task.get_post_result(_T, _P)['density']
            print(yaws, _T, _P, dens_exp, dens_sim)

            if dens_exp != None and dens_sim != None:
                dens_list.append([dens_exp, dens_sim])

        return dens_list

    def get_Cp(self, smiles_list, T=298):
        Cp_list = []
        for i, yaws in enumerate(self.yaws_list):
            task = self.task_list[i]
            _T, _P = self.get_T_P_from_T(yaws, T)

            if None in (_T, _P):
                continue

            cv = Cv.query.filter(Cv.cas == yaws.cas).first()
            if cv == None:
                continue

            Cp_exp = yaws.get_Cp(_T)
            post_result = task.get_post_result(_T, _P)
            Cp_inter = post_result['Cp_inter']
            Cp_PV = post_result['Cp_PV']
            Cp_sim = Cp_inter + Cp_PV + cv.get_post_result(_T)

            print(yaws, _T, _P, Cp_exp, Cp_sim)

            if Cp_exp != None and Cp_sim != None:
                Cp_list.append([Cp_exp, Cp_sim])

        return Cp_list

    def get_Hvap(self, smiles_list, T):
        Hvap_list = []
        for i, yaws in enumerate(self.yaws_list):
            task = self.task_list[i]
            _T, _P = self.get_T_P_from_T(yaws, T)

            if None in (_T, _P):
                continue

            Hvap_exp = yaws.get_Hvap(_T)
            Hvap_sim = task.get_post_result(_T, _P)['Hvap']

            print(yaws, _T, _P, Hvap_exp, Hvap_sim)

            if Hvap_exp != None and Hvap_sim != None:
                Hvap_list.append([Hvap_exp, Hvap_sim])

        return Hvap_list

    def get_expansion(self, smiles_list, T):
        expan_list = []
        for i, yaws in enumerate(self.yaws_list):
            task = self.task_list[i]
            _T, _P = self.get_T_P_from_T(yaws, T)

            if None in (_T, _P):
                continue

            expan_exp = yaws.get_expansion(_T)
            expan_sim = task.get_post_result(_T, _P)['expansion']

            print(_T, _P, expan_exp, expan_sim)

            if expan_exp != None and expan_sim != None:
                expan_list.append([expan_exp, expan_sim])

        return expan_list
