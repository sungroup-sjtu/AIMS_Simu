import os
from collections import OrderedDict
from typing import Dict


class PPF():
    def __init__(self, ppf_file):
        with open(ppf_file) as f:
            self.terms = f.read().splitlines()

    @property
    def adj_lj_paras(self) -> OrderedDict:
        terms = OrderedDict()
        for term in self.terms:
            if not (term.startswith('N12_6') or term.startswith('BINC')):
                continue

            words = term.split(':')
            words = [w.strip() for w in words]
            if term.startswith('N12_6'):
                a_type = words[1]
                paras = words[2]

                words = paras.split(',')
                words = [w.strip() for w in words]
                r0 = words[0]
                e0 = words[1]
                if not r0.endswith('*'):
                    terms['%s_r0' % a_type] = float(r0)
                if not e0.endswith('*'):
                    terms['%s_e0' % a_type] = float(e0)
            elif term.startswith('BINC'):
                a_types = words[1]
                para = words[2]

                words = a_types.split(',')
                words = [w.strip() for w in words]
                a1_type = words[0]
                a2_type = words[1]
                if not para.endswith('*'):
                    terms['%s_%s_bi' % (a1_type, a2_type)] = float(para)
        return terms

    def set_lj_para(self, new_paras: Dict):
        terms = [''] * len(self.terms)
        replaced = False
        for i, term in enumerate(self.terms):
            if not (term.startswith('N12_6') or term.startswith('BINC')):
                terms[i] = term
                continue

            if term.startswith('N12_6'):
                words = term.split(':')
                words = [w.strip() for w in words]
                a_type = words[1]
                paras = words[2]
                words = paras.split(',')
                words = [w.strip() for w in words]
                r0 = words[0]
                e0 = words[1]

                r0_key = '%s_r0' % a_type
                if r0_key in new_paras.keys():
                    r0 = new_paras[r0_key]
                    replaced = True
                e0_key = '%s_e0' % a_type
                if e0_key in new_paras.keys():
                    e0 = new_paras[e0_key]
                    replaced = True

                new_term = 'N12_6: %s: %s, %s:' % (a_type, str(r0), str(e0))
                terms[i] = new_term
            elif term.startswith('BINC'):
                words = term.split(':')
                words = [w.strip() for w in words]
                a_types = words[1]
                bi = words[2]

                words = a_types.split(',')
                words = [w.strip() for w in words]
                a1_type = words[0]
                a2_type = words[1]

                bi_key = '%s_%s_bi' % (a1_type, a2_type)
                if bi_key in new_paras.keys():
                    bi = new_paras[bi_key]
                    replaced = True

                new_term = 'BINC: %s, %s: %s:' % (a1_type, a2_type, str(bi))
                terms[i] = new_term

        self.terms = terms
        return replaced

    def write(self, ppf_out):
        with open(ppf_out, 'w') as f:
            for term in self.terms:
                f.write(term + '\n')
