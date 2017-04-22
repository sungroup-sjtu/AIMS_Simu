class Procedure:
    NPT = 'npt'
    NPT_CP = 'npt-cp'
    NPT_HVAP = 'npt-hvap'
    NVT_SINGLE = 'nvt-single'
    NPT_SOLVATION = 'npt-solvation'
    NVT_MSD = 'nvt-msd'
    NVT_VISCOSITY = 'nvt-viscosity'
    NVT_SLAB = 'nvt-slab'
    NPT_BINARY_SLAB = 'npt-binary-slab'
    choices = [NPT, NPT_CP, NPT_HVAP, NPT_SOLVATION, NVT_SINGLE, NVT_MSD, NVT_VISCOSITY, NVT_SLAB, NPT_BINARY_SLAB]
    T_RELEVANT = choices
    P_RELEVANT = [NPT, NPT_CP, NPT_SOLVATION, NVT_MSD, NVT_VISCOSITY, NPT_BINARY_SLAB]
