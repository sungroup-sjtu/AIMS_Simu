class Procedure:
    NPT = 'npt'
    NPT_SOLVATION = 'npt-solvation'
    NPT_HOV_IONIC_LIQUID = 'npt-hov-ionic-liquid'
    NVT = 'nvt'
    NVT_MSD = 'nvt-msd'
    NVT_VISCOSITY = 'nvt-viscosity'
    NVT_SLAB = 'nvt-slab'
    NPT_BINARY_SLAB = 'npt-binary-slab'
    choices = [NPT, NPT_SOLVATION, NPT_HOV_IONIC_LIQUID, NVT, NVT_MSD, NVT_VISCOSITY, NVT_SLAB, NPT_BINARY_SLAB]
    T_RELEVANT = choices
    P_RELEVANT = [NPT, NPT_SOLVATION, NVT, NVT_MSD, NVT_VISCOSITY, NPT_BINARY_SLAB]

