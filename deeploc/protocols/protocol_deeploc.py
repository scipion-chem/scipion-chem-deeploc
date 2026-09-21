# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import glob
import os
import pandas as pd

from pwem.protocols import EMProtocol
from pyworkflow.object import String
from pyworkflow.protocol import params
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder

from pwem.objects.data import AtomStruct, SetOfAtomStructs

from pwchem.__init__ import Plugin as pwchemPlugin

from deeploc import DEEPLOC_DIC


class ProtDeepLoc(EMProtocol):
  """"""
  _label = 'deeploc protein subcellular localization prediction'


  def _defineParams(self, form):
    form.addSection(label='Input')

    form.addParam('inputType', params.EnumParam, label='Input format: ',
                    default=0, choices=['AtomStruct', 'SetOfAtomStructs'],
                    help="Input format.")

    form.addParam('inputStruct', params.PointerParam, allowsNull=True,
                  pointerClass='AtomStruct', condition='inputType==0',
                  label="Input structure: ",
                  help='Select the reference structure.')
    form.addParam('inputSet', params.PointerParam, allowsNull=True,
                  pointerClass='SetOfAtomStructs', condition='inputType==1',
                  label="Input structures: ",
                  help='Select the reference structures.')

  def _insertAllSteps(self):
    self._insertFunctionStep(self.createInputFileStep)
    self._insertFunctionStep(self.runDeepLocStep)
    self._insertFunctionStep(self.createOutputStep)

  def createInputFileStep(self):
    fastaFile = self._getExtraPath("input.fasta")

    with open(fastaFile, "w") as f:
      if self.inputType.get() == 0:
        self.writeStructureFasta(self.inputStruct.get(), f)
      else:
        for atomStruct in self.inputSet.get():
          self.writeStructureFasta(atomStruct, f)


  def runDeepLocStep(self):
    fastaFile = self._getExtraPath("input.fasta")
    outFile = self._getPath('outputs')

    args = [
        '-f', os.path.abspath(fastaFile),
        '-o', os.path.abspath(outFile),
        '-p'
    ]

    pwchemPlugin.runCondaCommand(
        self,
        args=" ".join(args),
        condaDic=DEEPLOC_DIC,
        program="deeploc2",
        cwd=os.path.abspath(pwchemPlugin.getVar(DEEPLOC_DIC['home']))
    )

  def createOutputStep(self):
    outPath = self._getPath('outputs')
    outCsv = glob.glob(os.path.join(outPath, 'results*.csv'))
    df = pd.read_csv(outCsv[0])

    predictions = {}

    for _, row in df.iterrows():
        predictions[str(row['Protein_ID'])] = {
            'localizations': row['Localizations'],
            'signals': row['Signals'],
            'membraneTypes': row['Membrane types']
        }
    if self.inputType.get() == 0:
        model = self.inputStruct.get().clone()
        proteinId = os.path.splitext(
            os.path.basename(model.getFileName())
        )[0]
        prediction = predictions.get(proteinId)

        if prediction is None:
            raise ValueError(
                f"No DeepLoc prediction found for Protein_ID '{proteinId}'"
            )

        model._localizations = String(str(prediction['localizations']))
        model._signals = String(str(prediction['signals']))
        model._membraneTypes = String(str(prediction['membraneTypes']))

        model._localizationPerc = String(str(outCsv[0]))

        self._defineOutputs(outputAtomStruct=model)
    else:
        outputSet = SetOfAtomStructs.create(self._getPath())
        for atomStruct in self.inputSet.get():
            model = atomStruct.clone()
            proteinId = os.path.splitext(os.path.basename(model.getFileName()))[0]
            prediction = predictions.get(proteinId)

            if prediction is None:
                self.warning(
                    f"No DeepLoc prediction found for "
                    f"Protein_ID '{proteinId}'"
                )
                continue

            model._localizations = String(str(prediction['localizations']))
            model._signals = String(str(prediction['signals']))
            model._membraneTypes = String(str(prediction['membraneTypes']))

            outputSet._localizationPerc = String(str(outCsv[0]))

            outputSet.append(model)

        self._defineOutputs(outputAtomStructs=outputSet)

  ##################### UTILS #####################

  def _summary(self):
      summary = []
      if self.inputType.get() == 0:
          model = getattr(self, 'outputAtomStruct', None)
          if model is not None:
              if hasattr(model, '_localizations'):
                  summary.append(
                      f"*Localizations*: {model._localizations.get()}"
                  )
              if hasattr(model, '_signals'):
                  summary.append(
                      f"*Signals*: {model._signals.get()}"
                  )
              if hasattr(model, '_membraneTypes'):
                  summary.append(
                      f"*Membrane types*: {model._membraneTypes.get()}"
                  )
      else:
          outputSet = getattr(self, 'outputAtomStructs', None)
          if outputSet is not None:
              for i, model in enumerate(outputSet):
                  proteinId = os.path.splitext(
                      os.path.basename(model.getFileName())
                  )[0]
                  localizations = (
                      model._localizations.get()
                      if hasattr(model, '_localizations')
                      else 'N/A'
                  )
                  signals = (
                      model._signals.get()
                      if hasattr(model, '_signals')
                      else 'N/A'
                  )
                  membraneTypes = (
                      model._membraneTypes.get()
                      if hasattr(model, '_membraneTypes')
                      else 'N/A'
                  )
                  summary.append(
                      f"*{proteinId}*:\n"
                      f"    *Localizations*: {localizations}\n"
                      f"    *Signals*: {signals}\n"
                      f"    *Membrane types*: {membraneTypes}\n"
                  )
      return summary

  def writeStructureFasta(self, atomStruct, outHandle):
    fileName = atomStruct.getFileName()

    if fileName.endswith(".cif"):
        structure = MMCIFParser(QUIET=True).get_structure("protein", fileName)
    else:
        structure = PDBParser(QUIET=True).get_structure("protein", fileName)

    ppb = PPBuilder()
    seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(structure))
    if not seq:
        self.warning(f"Could not extract sequence from {fileName}")
        return
    seqId = os.path.splitext(os.path.basename(fileName))[0]
    outHandle.write(f">{seqId}\n")
    outHandle.write(f"{seq}\n")
