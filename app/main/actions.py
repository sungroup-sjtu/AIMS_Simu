from flask import request
import json
from ..models import *
from ..models_cv import *
from ..models_yaws import *
from ..models_nist import *


class StatAction():
    def update_mol_task_list(self, smiles_list):
        self.nist_list = []
        self.task_list = []
        self.TP_list = {
            # 298: [],
            'Tm25': [],
            'Tvap': [],
            'Tcx8': [],
        }

        for smiles in smiles_list:
            nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
            if nist == None:
                continue

            task = Task.query.filter(Task.smiles_list == json.dumps([smiles])).filter_by(procedure='npt').first()
            if task != None:
                post_result = task.get_post_result()
                if post_result == None \
                        or post_result['density-poly4-score'] < 0.999 \
                        or post_result['e_inter-poly4-score'] < 0.999:
                    task = None

            if task == None:
                continue

            self.nist_list.append(nist)
            self.task_list.append(task)

            for _T in self.TP_list.keys():
                T = self.get_T_nist(nist, _T)

                if T == None:
                    P = None
                else:
                    P = self.get_P_nist(nist, T)

                self.TP_list[_T].append((T, P))

    def update_mol_slab_list(self, smiles_list):
        self.nist_list = []
        self.slab_list = []
        self.TP_list = {
            # 298: [],
            'Tm25': [],
            'Tvap': [],
            'Tcx8': [],
        }

        for smiles in smiles_list:
            nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
            if nist == None:
                continue

            slab = Task.query.filter(Task.smiles_list == json.dumps([smiles])).filter_by(procedure='nvt-slab').first()
            if slab == None or slab.post_result == None:
                continue

            self.nist_list.append(nist)
            self.slab_list.append(slab)

            for _T in self.TP_list.keys():
                T = self.get_T_nist(nist, _T)

                if T == None:
                    P = None
                else:
                    P = self.get_P_nist(nist, T)

                self.TP_list[_T].append((T, P))

    def get_T_yaws(self, yaws, _T):
        if type(_T) == str and yaws == None:
            return None

        if _T == 'Tm25':
            T = yaws.get_Tfus()
            if T != None:
                T += 25
        elif _T == 'Tvap':
            T = yaws.get_Tvap()
        elif _T == 'Tcx8':
            T = yaws.get_Tc()
            if T != None:
                T *= 0.8
                T = min(T, 600)
        else:
            T = _T

        if T != None:
            T = int(round(T))

        return T

    def get_T_nist(self, nist, _T):
        if type(_T) == str and nist == None:
            return None

        if _T == 'Tm25':
            T = nist.tt
            if T != None:
                T += 25
        elif _T == 'Tvap':
            T = nist.tb
        elif _T == 'Tcx8':
            T = nist.tc
            if T != None:
                T *= 0.8
                T = min(T, 600)
        else:
            T = _T

        if T != None:
            T = int(round(T))

        return T

    def get_P_nist(self, nist, T):
        if nist == None:
            return None

        spline_pvap = nist.splines.join(NistProperty).filter(NistProperty.name == 'pvap-lg').first()
        if spline_pvap == None:
            P = None
        else:
            P = spline_pvap.get_data(T)[0]  # kPa
            if P != None:
                P /= 100  # bar
        return P

    def get_density_nist(self, _T=298):
        dens_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(T)
            if value == None:
                continue

            task = self.task_list[i]
            if task == None or task.post_result == None:
                continue
            dens_sim = task.get_post_data(T, P)['density']

            dens_list.append([value / 1000, dens_sim, uncertainty / 1000])

        return dens_list

    def get_cp_nist(self, _T=298):
        cp_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i%s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'cp-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(T)
            if value == None:
                continue

            cv = Cv.query.filter(Cv.smiles == nist.smiles).first()
            if cv == None:
                continue

            task = self.task_list[i]
            if task == None or task.post_result == None:
                continue
            post_data = task.get_post_data(T, P)
            cp_inter = post_data['cp_inter']
            cp_pv = post_data['cp_pv']
            cp_sim = cp_inter + cp_pv + cv.get_post_data(T)

            cp_list.append([value, cp_sim, uncertainty])

        return cp_list

    def get_sound_nist(self, _T=298):
        from scipy.interpolate import interp2d
        sound_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i%s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'sound-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(T)
            if value == None:
                continue

            cv = Cv.query.filter(Cv.smiles == nist.smiles).first()
            if cv == None:
                continue

            task = self.task_list[i]
            if task == None or task.post_result == None:
                continue
            post_data = task.get_post_data(T, P)
            cp_inter = post_data['cp_inter']
            cp_pv = post_data['cp_pv']
            cp_sim = cp_inter + cp_pv + cv.get_post_data(T)

            dens = post_data['density']  # g/mL
            expan = post_data['expansion']  # /K

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
            compr = interp_c(T, P)[0]  # /bar
            # compr = post_data['compressibility']  # /bar

            m = pybel.readstring('smi', nist.smiles)
            cv_sim = cp_sim - T * m.molwt / dens * expan ** 2 / compr * 0.1  # J/mol.K
            sound_sim = (cp_sim / cv_sim / compr / dens) ** 0.5 * 10  # m/s

            print(nist.smiles, T, P, cp_sim, cv_sim, sound_sim, value)

            sound_list.append([value, sound_sim, uncertainty])

        return sound_list

    def get_hvap_nist(self, _T=298):
        hvap_list = []
        for i, nist in enumerate(self.nist_list):
            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            if nist == None:
                continue

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'hvap-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(T)
            if value == None:
                continue

            task = self.task_list[i]
            if task == None or task.post_result == None:
                continue
            hvap_sim = task.get_post_data(T, P)['hvap-spline']
            # if _P < 2:
            #     hvap_sim = task.get_post_data(_T, _P)['hvap-spline']
            # else:
            #     ei_sim = task.get_post_data(_T, _P)['e_inter']
            #     dens_sim = task.get_post_data(_T, _P)['density']  # g/mL
            #
            #     spline_dgas = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
            #     if None in [spline_dgas]:
            #         continue
            #
            #     v_dgas, _ = spline_dgas.get_data(_T)  # kg/m^3
            #     if None in [v_dgas]:
            #         continue
            #
            #     pVg_l = _P * nist.weight * (1 / v_dgas - 0.001 / dens_sim) / 10  # kJ/mol
            #     hvap_sim = pVg_l - ei_sim

            dev = abs(hvap_sim - value) / value * 100
            if dev > 50:
                print(nist, task)

            # hvap_list.append([value, hvap_sim, uncertainty])
            # hvap_list.append([value, hvap_sim - 1.4 * (8.314 * T) / 1000 + 1.2 * math.log(nist.n_heavy), uncertainty])
            hvap_list.append([value, hvap_sim - nist.n_CON / 15 * (8.314 * T) / 1000, uncertainty])

        return hvap_list

    def get_tc_nist(self):
        tc_list = []
        for i, nist in enumerate(self.nist_list):
            if nist == None:
                continue
            slab = self.slab_list[i]
            if slab == None or slab.post_result == None:
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
            if slab == None or slab.post_result == None:
                continue

            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            value, uncertainty = nist.dc, nist.dc_u
            if value == None:
                continue

            dc = slab.get_post_data(100)['dc']

            dc_list.append([value / 1000, dc, uncertainty / 1000])

        return dc_list

    def get_dgas_nist(self, _T=298):
        dgas_list = []
        for i, nist in enumerate(self.nist_list):
            if nist == None:
                continue
            slab = self.slab_list[i]
            if slab == None or slab.post_result == None:
                continue

            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            T, _ = self.TP_list[_T][i]
            if None in (T,):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(T)
            if value == None:
                continue

            dgas_sim = slab.get_post_data(T)['density_gas']

            dgas_list.append([value / 1000, dgas_sim, uncertainty / 1000])

        return dgas_list

    def get_st_nist(self, _T=298):
        st_list = []
        for i, nist in enumerate(self.nist_list):
            if nist == None:
                continue
            slab = self.slab_list[i]
            if slab == None or slab.post_result == None:
                continue

            sys.stdout.write('\r\t%5i %s' % (i, nist))
            sys.stdout.flush()

            T, _ = self.TP_list[_T][i]
            if None in (T,):
                continue

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'st-lg').first()
            if spline == None:
                continue
            value, uncertainty = spline.get_data(T)
            if value == None:
                continue

            st_sim = slab.get_post_data(T)['surface_tension']

            dev = (st_sim / 1000 - value) / value * 100
            if dev > 50:
                print(nist, slab)

            st_list.append([value * 1000, st_sim, uncertainty * 1000])

        return st_list


class VleAction():
    def get_vle(self, smiles) -> bool:
        import numpy as np

        self.nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
        slab = Task.query.filter(Task.smiles_list == json.dumps([smiles])).filter_by(procedure='nvt-slab').first()
        if self.nist == None or slab == None or slab.post_result == None:
            return False

        spline_dl = self.nist.splines.join(NistProperty).filter(NistProperty.name == 'density-lg').first()
        spline_dg = self.nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
        spline_pvap = self.nist.splines.join(NistProperty).filter(NistProperty.name == 'pvap-lg').first()
        if spline_dl == None and spline_dg == None:
            return False

        self.t_list = np.linspace(slab.t_min, slab.t_max / 0.85 * 0.9, 20)
        self.dl_list = []
        self.dg_list = []
        self.dl_u_list = []
        self.dg_u_list = []
        self.log10pvap_list = []

        for t in self.t_list:
            if spline_dl == None:
                dl, dl_u = None, None
            else:
                dl, dl_u = spline_dl.get_data(t)
            if dl != None:
                dl /= 1000
                dl_u /= 1000

            if spline_dg == None:
                dg, dg_u = None, None
            else:
                dg, dg_u = spline_dg.get_data(t)
            if dg != None:
                dg /= 1000
                dg_u /= 1000

            if spline_pvap == None:
                pvap, pvap_u = None, None
            else:
                pvap, pvap_u = spline_pvap.get_data(t)
            if pvap != None:
                pvap /= 100
                pvap = np.log10(pvap)

            self.dl_list.append(dl)
            self.dl_u_list.append(dl_u)
            self.dg_list.append(dg)
            self.dg_u_list.append(dg_u)
            self.log10pvap_list.append(pvap)

        self.dc = self.nist.dc
        if self.dc != None:
            self.dc /= 1000
        self.tc = self.nist.tc

        self.ts_list = []
        self.dls_list = []
        self.ts_good_list = []
        self.dgs_good_list = []
        self.ts_bad_list = []
        self.dgs_bad_list = []
        self.log10pvaps_good_list = []
        self.log10pvaps_u_good_list = []
        for job in slab.jobs:
            result = job.get_result()
            if result == None or result.get('failed') == True:
                continue
            self.ts_list.append(job.t)
            if result['density_gas'][0] >= 0.01:  # density of gas phase smaller than 0.01 g/mL are not reliable
                self.ts_good_list.append(job.t)
                log10pvaps = np.log10(result['pzz'][0])
                log10pvaps_min = np.log10(result['pzz'][0] - result['pzz'][1])
                self.log10pvaps_good_list.append(log10pvaps)
                self.log10pvaps_u_good_list.append(log10pvaps - log10pvaps_min)
            else:
                self.ts_bad_list.append(job.t)

        for i, t in enumerate(self.ts_list):
            post_data = slab.get_post_data(t)
            self.dls_list.append(post_data['density_liq'])
            if t in self.ts_good_list:
                self.dgs_good_list.append(post_data['density_gas'])
            else:
                self.dgs_bad_list.append(post_data['density_gas'])
            self.tcs = post_data['tc']
            self.dcs = post_data['dc']

        return True
