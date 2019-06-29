from mstools.utils import get_T_list_from_range, random_string

def get_t_min_max(Tvap=None, Tfus=None, Tc=None):
    if Tvap is None:
        return None, None

    if Tfus is None:
        t_min = int(round(Tvap * 0.4 + 100))
    else:
        t_min = int(round(Tfus + 25))

    if Tc is None:
        t_max = int(round(Tvap * 1.2))
    else:
        t_max = int(round(Tc * 0.85))
    t_max = min(t_max, 650)

    if t_min >= t_max:
        return None, None
    return t_min, t_max

def get_sim_t_list(Tvap=None, Tm=None, Tc=None):
    t_min, t_max = get_t_min_max(Tvap=Tvap, Tfus=Tm, Tc=Tc)
    if t_min is None or t_max is None:
        return None
    return get_T_list_from_range(t_min, t_max)
