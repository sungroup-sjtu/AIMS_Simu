import sys
sys.path.append('..')
from app import create_app
from app.models_nist import *
from app.selection import ClassificationAtomType
from rdkit.Chem import AllChem as Chem
from mstools.smiles.smiles import has_stereo_isomer_from_inchi
procedure = 'npt'
app = create_app(procedure)
app.app_context().push()
molecules = NistMolecule.query.filter(NistMolecule.id < 85)
classifier = ClassificationAtomType(AllowedAtomicNum=[1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53])
for i, mol in enumerate(molecules):
    sys.stdout.write('\r%i / %i. %s\t\t\t\t' % (i, molecules.count(), mol.smiles))
    rdk_mol = Chem.MolFromInchi(mol.inchi)
    if '.' in mol.smiles:
        mol.remark = 'mixture'
    elif rdk_mol is None:
        mol.remark = 'unrecognizable'
    elif classifier.classify(mol.inchi):
        if has_stereo_isomer_from_inchi(mol.inchi):
            mol.remark = 'has undefined stereo'
        else:
            mol.remark = 'selected'
    db.session.commit()
