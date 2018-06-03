from flask import request
import json
from ..models import *
from ..models_cv import *
from ..models_yaws import *
from ..models_nist import *


class StatAction():
    def update_mol_task_list(self, smiles_list):
        self.yaws_list = []
        self.nist_list = []
        self.task_list = []
        self.slab_list = []
        self.TP_list = {298: [], 'Tm25': [], 'Tvap': [], 'Tcx8': []}

        for smiles in smiles_list:
            yaws = YawsMolecule.query.filter(YawsMolecule.isomeric_smiles == smiles).first()
            nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
            if yaws == None and nist == None:
                continue

            # ignore molecules without expt data
            # if not yaws.density_is_exp and not yaws.Cp_is_exp and not yaws.Hvap_is_exp and not yaws.expansion_is_exp:
            #     continue

            task = Task.query.filter(Task.smiles_list == json.dumps([smiles])).filter_by(procedure='npt').first()
            slab = Task.query.filter(Task.smiles_list == json.dumps([smiles])).filter_by(procedure='nvt-slab').first()
            if (task == None or task.post_result == None) and (slab == None or slab.post_result == None):
                continue

            self.yaws_list.append(yaws)
            self.nist_list.append(nist)
            self.task_list.append(task)
            self.slab_list.append(slab)

            T, P = self.get_T_P_from_T(yaws, 298)
            self.TP_list[298].append((T, P))
            T, P = self.get_T_P_from_T(yaws, 'Tm25')
            self.TP_list['Tm25'].append((T, P))
            T, P = self.get_T_P_from_T(yaws, 'Tvap')
            self.TP_list['Tvap'].append((T, P))
            T, P = self.get_T_P_from_T(yaws, 'Tcx8')
            self.TP_list['Tcx8'].append((T, P))

    def get_T_P_from_T(self, yaws, T):
        if T == 'Tm25':
            T = yaws.get_Tfus()
            if T != None:
                T += 25
        if T == 'Tvap':
            T = yaws.get_Tvap()
        elif T == 'Tcx8':
            T = yaws.get_Tc()
            if T != None:
                T *= 0.8
                T = min(T, 600)
        if T == None:
            return None, None

        Pvap = yaws.get_Pvap(T)
        return T, Pvap

    def get_density(self, T=298):
        dens_list = []
        for i, yaws in enumerate(self.yaws_list):
            if yaws == None:
                continue
            if not yaws.density_is_exp:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            print(yaws, _T, _P)
            task = self.task_list[i]
            dens_exp = yaws.get_density(_T)
            dens_sim = task.get_post_data(_T, _P)['density']
            if dens_exp != None:
                dens_list.append([dens_exp, dens_sim])

        return dens_list

    def get_Cp(self, T=298):
        Cp_list = []
        for i, yaws in enumerate(self.yaws_list):
            task = self.task_list[i]
            if not yaws.Cp_is_exp:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            cv = Cv.query.filter(Cv.cas == yaws.cas).first()
            if cv == None:
                continue

            Cp_exp = yaws.get_Cp(_T)
            post_result = task.get_post_data(_T, _P)
            Cp_inter = post_result['cp_inter']
            Cp_PV = post_result['cp_pv']
            Cp_sim = Cp_inter + Cp_PV + cv.get_post_data(_T)

            print(yaws, _T, _P, Cp_exp, Cp_sim)

            if Cp_exp != None and Cp_sim != None:
                Cp_list.append([Cp_exp, Cp_sim])

        return Cp_list

    def get_density_nist(self, T=298):
        dens_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(_T)
            if value == None:
                continue

            task = self.task_list[i]
            dens_sim = task.get_post_data(_T, _P)['density']

            dens_list.append([value / 1000, dens_sim, uncertainty / 1000])

        return dens_list

    def get_cp_nist(self, T=298):
        cp_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i%s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'cp-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(_T)
            if value == None:
                continue

            cv = Cv.query.filter(Cv.smiles == nist.smiles).first()
            if cv == None:
                continue

            task = self.task_list[i]
            post_result = task.get_post_data(_T, _P)
            cp_inter = post_result['cp_inter']
            cp_pv = post_result['cp_pv']
            cp_sim = cp_inter + cp_pv + cv.get_post_data(_T)

            cp_list.append([value, cp_sim + 8.314 * nist.n_heavy / 15, uncertainty])

        return cp_list

    def get_sound_nist(self, T=298):
        from scipy.interpolate import interp2d
        sound_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i%s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'sound-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(_T)
            if value == None:
                continue

            cv = Cv.query.filter(Cv.smiles == nist.smiles).first()
            if cv == None:
                continue

            task = self.task_list[i]
            post_result = task.get_post_data(_T, _P)
            cp_inter = post_result['cp_inter']
            cp_pv = post_result['cp_pv']
            cp_sim = cp_inter + cp_pv + cv.get_post_data(_T)

            dens = post_result['density']  # g/mL
            expan = post_result['expansion']  # /K

            t_list = []
            p_list = []
            c_list = []
            for job in task.jobs:
                if not job.converged:
                    continue
                t_list.append(job.t)
                p_list.append(job.p)
                c_list.append(json.loads(job.result)['compressibility'][0])
            interp_c = interp2d(t_list, p_list, c_list)
            compr = interp_c(_T, _P)[0]  # /bar
            # compr = post_result['compressibility']  # /bar

            m = pybel.readstring('smi', nist.smiles)
            cv_sim = cp_sim - _T * m.molwt / dens * expan ** 2 / compr * 0.1  # J/mol.K
            sound_sim = (cp_sim / cv_sim / compr / dens) ** 0.5 * 10  # m/s

            print(nist.smiles, _T, _P, cp_sim, cv_sim, sound_sim, value)

            sound_list.append([value, sound_sim, uncertainty])

        return sound_list

    def get_hvap(self, T):
        hvap_list = []
        for i, yaws in enumerate(self.yaws_list):
            sys.stdout.write('\r\t%5i %s' % (i, yaws))
            sys.stdout.flush()

            if yaws == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            if not yaws.Hvap_is_exp:
                continue
            hvap_exp = yaws.get_Hvap(_T)

            task = self.task_list[i]
            hvap_sim = task.get_post_data(_T, _P)['hvap']

            if hvap_exp != None and hvap_sim != None:
                hvap_list.append([hvap_exp, hvap_sim, 0])

        return hvap_list

    def get_hvap_nist(self, T=298):
        hvap_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'hvap-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(_T)
            if value == None:
                continue

            if math.isnan(value):
                continue

            task = self.task_list[i]
            if _P < 2:
                hvap_sim = task.get_post_data(_T, _P)['hvap']
            else:
                ei_sim = task.get_post_data(_T, _P)['e_inter']
                dens_sim = task.get_post_data(_T, _P)['density']  # g/mL

                spline_dgas = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
                if None in [spline_dgas]:
                    continue

                v_dgas, _ = spline_dgas.get_data(_T)  # kg/m^3
                if None in [v_dgas]:
                    continue

                pVg_l = _P * nist.weight * (1 / v_dgas - 0.001 / dens_sim) / 10  # kJ/mol
                hvap_sim = pVg_l - ei_sim

            hvap_list.append([value, hvap_sim - nist.n_heavy / 15 * 8.314 * _T / 1000, uncertainty])

        return hvap_list

    def get_ei_nist(self, T=298):
        ei_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            spline_hvap = nist.splines.join(NistProperty).filter(NistProperty.name == 'hvap-lg').first()
            spline_dliq = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-lg').first()
            spline_dgas = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
            if None in [spline_hvap, spline_dliq, spline_dgas]:
                continue
            v_hvap, uncertainty = spline_hvap.get_data(_T)
            if v_hvap == None or math.isnan(v_hvap):
                continue

            if _P < 2:
                pVg_l = 8.314 * _T / 1000  # ideal gas approximation
            else:
                v_dliq, _ = spline_dliq.get_data(_T)
                if v_dliq == None:
                    continue
                v_dgas, _ = spline_dgas.get_data(_T)
                if v_dgas == None:
                    continue
                pVg_l = _P * nist.weight * (1 / v_dgas - 1 / v_dliq) / 10  # kJ/mol

            task = self.task_list[i]
            ei_sim = task.get_post_data(_T, _P)['e_inter']

            ei_list.append([pVg_l - v_hvap, ei_sim + nist.n_heavy / 15 * 8.314 * _T / 1000, -uncertainty])

        return ei_list

    def get_pv_nist(self, T=298):
        pv_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            spline_dliq = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-lg').first()
            spline_dgas = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
            if None in [spline_dliq, spline_dgas]:
                continue
            v_dliq, _ = spline_dliq.get_data(_T)
            if v_dliq == None:
                continue
            v_dgas, _ = spline_dgas.get_data(_T)
            if v_dgas == None:
                continue
            pVg_l = _P * nist.weight * (1 / v_dgas - 1 / v_dliq) / 10  # kJ/mol
            RT = 8.314 * _T / 1000

            pv_list.append([pVg_l, RT, 0])

        return pv_list

    def get_tc_nist(self):
        tc_list = []
        for i, nist in enumerate(self.nist_list):
            if nist == None:
                continue
            slab = self.slab_list[i]
            if slab == None:
                continue

            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            value, uncertainty = nist.tc, nist.tc_u
            if value == None:
                continue

            tc = slab.get_post_data(100)['tc']

            tc_list.append([value, tc, uncertainty])

        return tc_list

    def get_dc_nist(self):
        dc_list = []
        for i, nist in enumerate(self.nist_list):
            if nist == None:
                continue
            slab = self.slab_list[i]
            if slab == None:
                continue

            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            value, uncertainty = nist.dc, nist.dc_u
            if value == None:
                continue

            dc = slab.get_post_data(100)['dc']

            dc_list.append([value / 1000, dc, uncertainty / 1000])

        return dc_list

    def get_dgas_nist(self, T=298):
        dgas_list = []
        for i, nist in enumerate(self.nist_list):
            if nist == None:
                continue
            slab = self.slab_list[i]
            if slab == None:
                continue

            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            _T, _P = self.TP_list[T][i]
            if None in (_T,):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(_T)
            if value == None:
                continue

            dgas_sim = slab.get_post_data(_T)['density_gas']

            dgas_list.append([value / 1000, dgas_sim, uncertainty / 1000])

        return dgas_list

    def get_st_nist(self, T=298):
        st_list = []
        for i, nist in enumerate(self.nist_list):
            if nist == None:
                continue
            slab = self.slab_list[i]
            if slab == None:
                continue

            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            _T, _P = self.TP_list[T][i]
            if None in (_T,):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'st-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(_T)
            if value == None:
                continue

            slab = self.slab_list[i]
            st_sim = slab.get_post_data(_T)['surface_tension']

            dev = (st_sim / 1000 - value) / value * 100
            if dev > 50:
                print(nist, slab)

            st_list.append([value * 1000, st_sim, uncertainty * 1000])

        return st_list

    def get_expansion(self, T):
        expan_list = []
        for i, yaws in enumerate(self.yaws_list):
            sys.stdout.write('\r\t%5i %s' % (i, yaws))
            sys.stdout.flush()

            if yaws == None:
                continue

            _T, _P = self.TP_list[T][i]
            if None in (_T, _P):
                continue

            if not yaws.expansion_is_exp:
                continue
            expan_exp = yaws.get_expansion(_T)

            task = self.task_list[i]
            expan_sim = task.get_post_data(_T, _P)['expansion']

            if expan_exp != None and expan_sim != None:
                expan_list.append([expan_exp, expan_sim, 0])

        return expan_list
