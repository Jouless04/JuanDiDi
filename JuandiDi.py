
from __future__ import print_function
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit import rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import MolsToImage
from IPython.display import display
from pubchempy import get_compounds as compounds
from PIL import Image
#Librerias adicionales 
import pandas as pd
import numpy as np
from itertools import chain
from IPython.display import Image
from 
IPythonConsole.ipython_useSVG=False

comp = ['alpha-Maltotriose','Amylopectin','Hydrochloric Acid']


#print(comp)
smiles = []
for i in comp: 
    try:
        compi = compounds(i,'name')[0]
        smiles.append(compi.canonical_smiles)
       
        #print(smiles)
    except IndexError:
        smiles.append('None')

mol_list=[]    
for sm in smiles:
    mol = Chem.MolFromSmiles(sm)
    mol_list.append(mol)
img = MolsToImage(mol_list)


# Crea un objeto MolDraw2DCairo para dibujar la estructura molecular
drawer = MolDraw2DCairo.MolDraw2DCairo(300, 300)
drawer.drawMolecule(mol)
drawer.FinishDrawing()

# Crea un objeto de imagen PNG utilizando valores RGB y alfa
image = Image.frombytes('RGBA', (300, 300), drawer.GetDrawingText())

# Mostrar la imagen en una ventana emergente
image.show()

display(img)
reactivo1 = smiles[1]
reactivo2 = smiles[2]
#reactivo_HCl = smiles[3]

#rxn = AllChem.ReactionFromSmarts
#reaccion = 
print(smiles)
