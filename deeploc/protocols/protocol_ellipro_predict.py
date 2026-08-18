# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, json
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.convert import AtomicStructHandler

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs, StructROI, SetOfStructROIs
from pwchem.utils import runOpenBabel, createPocketFile, pdbFromASFile, getBaseName, cifFromASFile

from .. import Plugin as iedbPlugin
from ..constants import ELLI_DIC

class ProtElliProPrediction(EMProtocol):
  """Run a prediction using ElliPro to extract B-cell structural epitopes"""
  _label = 'ellipro prediction'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                    label='Input protein structure: ',
                    help="Protein structure to perform the prediction on")

    iGroup.addParam("rchains", params.BooleanParam,
                    label='Select specific chains: ', default=False, important=True,
                    help='Keep only the chains selected')
    iGroup.addParam("chain_name", params.StringParam,
                    label="Keep chains: ", important=True, condition="rchains",
                    help="Select the chain(s) you want to use in the analysis")

    pGroup = form.addGroup('Parameters')
    pGroup.addParam('minScore', params.FloatParam, label='Minimum score: ', default=0.5,
                    help="Minimum score for a predicted epitope to be considered")
    pGroup.addParam('maxDist', params.FloatParam, label='Maximum distance (A): ', default=6.0,
                    help="Maximum distance between the epitope points")


  def _insertAllSteps(self):
    self._insertFunctionStep(self.convertInputStep)
    self._insertFunctionStep(self.elliProStep)
    self._insertFunctionStep(self.createOutputStep)

  def convertInputStep(self):
    inpFile = self.getInputPath()
    pdbFromASFile(inpFile, self._getPdbFile(), atomStruct=self.inputAtomStruct.get())

  def elliProStep(self):
    chains = self.getInputChains()
    elliArgs = f' -f {self._getPdbFile()} -s {self.minScore.get()} -d {self.maxDist.get()} '
    if chains:
      elliArgs += f'-c {",".join(chains)}'
    iedbPlugin.runElliPro(self, elliArgs, cwd=os.path.abspath(self._getExtraPath()))

  def createOutputStep(self):
    epiDic = self.parseResults()
    handler = AtomicStructHandler(self._getPdbFile())
    structModel = PDBParser().get_structure('inputAS', self._getPdbFile())[0]

    outSeqROIOuts = {}
    for chainId, chainDic in epiDic['Linear'].items():
      mapDic = self.getResidueIdMap(structModel, chainId)

      outSeqROIOuts[chainId] = SetOfSequenceROIs(filename=self._getPath(f'sequenceROIs_{chainId}.sqlite'))
      chainSeq = str(handler.getSequenceFromChain(modelID=0, chainID=chainId))
      outSeq = Sequence(name=f'Chain {chainId}', sequence=chainSeq, description=f'Chain {chainId}')

      for roiId, roiDic in chainDic.items():
        idxs = roiDic['Index']
        idxs = [mapDic[i] for i in idxs]
        roiSeq = Sequence(sequence=roiDic['Peptide'], name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                          description='Linear ElliPro epitope')
        seqROI = SequenceROI(sequence=outSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
        seqROI._epitopeType = params.String('B')
        seqROI._source = params.String('ElliPro')
        seqROI._type = 'ElliPro'
        setattr(seqROI, 'ElliPro', params.Float(roiDic['Score']))
        outSeqROIOuts[chainId].append(seqROI)

    for chainId, outROIs in outSeqROIOuts.items():
      self._defineOutputs(**{f'outSequenceROIs_{chainId}': outROIs})


    cifProtFile = cifFromASFile(self.getInputPath(), self._getCifFile(), atomStruct=self.inputAtomStruct.get())
    outROIs = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
    for roiId, roiDic in epiDic['Discontinous'].items():
      roiCoords = []
      for res in roiDic['Residues']:
        chainId, resIdx = res.split(':')[0], int(res.split(':')[1][1:])
        residue = structModel[chainId][resIdx]
        atoms = residue.get_atoms()
        for a in atoms:
          roiCoords.append(list(a.get_coord()))

      pocketFile = self._getExtraPath(f'pocketFile_{roiId}.cif')
      createPocketFile(roiCoords, roiId, pocketFile)
      pocket = StructROI(pocketFile, cifProtFile)
      pocket.calculateContacts()
      outROIs.append(pocket)

    if len(outROIs) > 0:
      outROIs.buildPDBhetatmFile()
      self._defineOutputs(outputStructROIs=outROIs)

  ##################### UTILS #####################
  def getResidueIdMap(self, structModel, chainId):
    resIdxs = {}
    for i, res in enumerate(structModel[chainId]):
      resIdxs[res.get_id()[1]] = i+1
    return resIdxs

  def getInputPath(self):
    return self.inputAtomStruct.get().getFileName()

  def getInputFileName(self):
    return self.getInputPath().split('/')[-1]

  def _getPdbFile(self):
    return os.path.abspath(self._getExtraPath(self._getInputName() + '.pdb'))

  def _getCifFile(self):
    return os.path.abspath(self._getExtraPath(self._getInputName() + '.cif'))

  def _getInputName(self):
    return getBaseName(self.getInputPath())

  def getInputChains(self):
    chains = []
    if self.rchains.get() and self.chain_name.get().strip():
      chDic = json.loads(self.chain_name.get())
      if 'residues' in chDic:
        chains = [chDic['chain'].strip()]
      elif 'model-chain' in chDic:
        mcs = chDic['model-chain'].split(',')
        chains = [mc.split('-')[1].strip() for mc in mcs]
    return chains

  def parseResults(self):
    rDic = {'Linear': {}, 'Discontinous': {}}
    with open(self._getExtraPath('output.txt')) as f:
      for line in f:
        line = line.strip()
        if line:
          sl = line.split(',')
          mode = sl[-1]
          if mode == 'Linear':
            chain = sl[2]
            if chain not in rDic['Linear']:
              rDic['Linear'][chain] = {}
            rDic['Linear'][chain][sl[0]] = {'Index': [int(sl[3]), int(sl[4])], 'Peptide': sl[5], 'Score': float(sl[7])}
          elif mode == 'Discontinous':
            rDic['Discontinous'][sl[0]] = {'Residues': [res.strip() for res in sl[2:-3]], 'Score': float(sl[-2])}

    return rDic



