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

import os

from pyworkflow.protocol import params

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs

from .. import Plugin as iedbPlugin
from ..constants import MHCI_alleles_dic
from ..utils import getAllMHCIAlleles
from ..protocols.protocol_mhc_ii_predict import ProtMHCIIPrediction

SEQ, SEQROIS = 0, 1
RANK, IC50, TOPP, NTOP = 0, 1, 2, 3

class ProtMHCIPrediction(ProtMHCIIPrediction):
  """Run a prediction using mhc-i package from IEDB to predict MHC-I epitopes over a sequence
  or to label a set of sequence ROIs with the alleles found on them"""
  _label = 'mhc-i prediction'

  MINLEN, MAXLEN = 8, 14
  selMap = {RANK: 'rank', IC50: 'ic50', TOPP: 'topPerc', NTOP: 'topN'}

  _mhciMethodsDic = {'IEDB recommended': 'netmhcpan', 'Consensus-2.18': 'consensus',
                     'NetMHC_Cons': 'netmhccons', 'ANN-4.0': 'ann', 'SMMPMBEC-1.0': 'smmpmbec', 'SMM-1.0': 'smm',
                     'Combinatorial Library-1.0': 'comblib_sidney2008', 'PickPocket-1.1': 'pìckpocket'}
  _species = ['Chimpanzee', 'Cow', 'Gorilla', 'Human', 'Macaque', 'Mouse', 'Pig']
  _alleleGroups = ['Frequent (>1%)', 'Representative HLA supertypes', 'Most frequent A, B', 'Custom']
  _selTypes = ['Percentile rank', 'IC50', 'Top x%', 'Top x']

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = self._defineInputParams(form)

    pGroup = form.addGroup('Parameters')
    pGroup.addParam('method', params.EnumParam, label='Prediction method: ',
                    default=0, choices=list(self._mhciMethodsDic.keys()),
                    help="Prediction method to use for the MHC-I binding prediction over the protein sequence.")
    pGroup.addParam('predMode', params.EnumParam, label='Prediction mode: ', condition='method==0',
                    default=0, choices=['Binding Affinity', 'Elution'],
                    help="Whether to predict the binding affinity or an elution score.")

    pGroup.addParam('specie', params.EnumParam, label='Host species: ', expertLevel=params.LEVEL_ADVANCED,
                    default=3, choices=self._species,
                    help="Host specie  to predict the MHC-I epitopes on.")
    pGroup.addParam('lengths', params.StringParam, label='Peptide lengths: ', default='9',
                    help="Available lengths to include in the analysis."
                         "You can include several lengths as comma separated  (11, 12,13)"
                         f"Lenght limits are [{self.MINLEN}, {self.MAXLEN}]")

    pGroup.addParam('alleleGroup', params.EnumParam, label='Select allele groups: ', condition='specie==3',
                    default=2, choices=self._alleleGroups,
                    help="Select usual groups of human alleles to perform the analysis on")
    pGroup.addParam('alleleCustom', params.StringParam, label='Select alleles to add: ',
                    condition='specie!=3 or alleleGroup==3',
                    help="Host specie  to predict the MHC-I epitopes on.")

    sGroup = form.addGroup('Selection')
    sGroup.addParam('mergeAlleles', params.BooleanParam, label='Merge alleles: ', default=True,
                    expertLevel=params.LEVEL_ADVANCED, condition='inputSource==0',
                    help="Merges same epitope sequences predicted to interact with different alleles into the same "
                         "sequence ROI")

    sGroup.addParam('selType', params.EnumParam, label='Select output peptides by: ',
                    default=0, choices=self._selTypes,
                    help="Select output peptides in the chosen manner")
    sGroup.addParam('rank', params.FloatParam, label='Percentile rank threshold: ', default=1, condition='selType==0',
                    help="Predicted percentile rank threshold (<=)")
    sGroup.addParam('ic50', params.FloatParam, label='IC50 threshold (nM): ', default=500, condition='selType==1',
                    help="Predicted IC50 threshold (<=)")
    sGroup.addParam('topPerc', params.FloatParam, label='Top percentage threshold (%): ', default=2,
                    condition='selType==2', help="Select top x% peptides based on percentile rank")
    sGroup.addParam('topN', params.IntParam, label='Top threshold: ', default=5, condition='selType==3',
                    help="Select top x peptides based on percentile rank")


  def _insertAllSteps(self):
    self._insertFunctionStep(self.mhcStep)
    self._insertFunctionStep(self.createOutputStep)

  def mhcStep(self):
    inFile = self.writeInputFasta()
    oFile = self.getMHCOutputFile()

    method = self._mhciMethodsDic[self.getEnumText('method')]
    if self.method.get() == 0:
        pMode = 'ba' if self.predMode.get() == 0 else 'el'
        method += f'_{pMode}'
    selAlleles = self.getSelectedAlleles()
    lenList = self.getLenghts()

    alDic = getAllMHCIAlleles(method, specie=self.getEnumText('specie').lower())
    fullAlList, fullLenList = self.filterAlleles(alDic, selAlleles, lenList)
    fullAlStr, fullLenStr = ','.join(fullAlList), ','.join(fullLenList)

    mhcArgs = f'{method} {fullAlStr} {fullLenStr} {inFile} > {oFile} '

    iedbPlugin.runMHC_I(self, mhcArgs)

  def createOutputStep(self):
    epiDic = self.parseResults(self.getMHCOutputFile())

    inpSeq = self.inputSequence.get()
    outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
    method = f"MHCI_{self.getEnumText('method')}"

    if self.inputSource.get() == SEQ:
      epiDic = epiDic['1']
      for (idxI, epitope) in epiDic:
        idxs = [idxI, idxI + len(epitope) - 1]
        roiSeq = Sequence(sequence=epitope, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                          description=f'MHC-I TepiTool epitope')

        if not self.mergeAlleles.get():
          for allele, score in epiDic[(idxI, epitope)].items():
            seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
            seqROI._allelesMHCI = params.String(allele)
            seqROI._epitopeType = params.String('MHC-I')
            seqROI._source = params.String(method)
            setattr(seqROI, method, params.Float(score))

            outROIs.append(seqROI)
        else:
          alleles, scores = list(epiDic[(idxI, epitope)].keys()), list(epiDic[(idxI, epitope)].values())
          allele, score = '/'.join(alleles), min(scores)
          seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
          seqROI._allelesMHCI = params.String(allele)
          seqROI._epitopeType = params.String('MHC-I')
          seqROI._source = params.String(method)
          setattr(seqROI, method, params.Float(score))

          outROIs.append(seqROI)

    else:
      # Each input ROI is labelled with the alleles found inside them
      inROIs = [roi.clone() for roi in self.inputSequenceROIs.get()]
      i, lens = 0, self.getLenghts()
      for curROI in inROIs:
        curAlleles, curScores = [], []
        if len(curROI.getROISequence()) >= min(lens):
          i += 1
          roiId = str(i)
          if roiId in epiDic:
            for (idxI, epitope) in epiDic[roiId]:
              alleles, scores = list(epiDic[roiId][(idxI, epitope)].keys()), list(epiDic[roiId][(idxI, epitope)].values())
              curAlleles += alleles
              curScores += scores

        allele, score = '/'.join(curAlleles), min(curScores) if curScores else 0
        curROI._allelesMHCI = params.String(allele)
        curROI._sourceMHCI = params.String(method)
        setattr(curROI, method, params.Float(score))
        outROIs.append(curROI)

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)

  ##################### UTILS #####################

  def getAvailableAlleles(self):
    methKey = self._mhciMethodsDic[self.getEnumText('method')]
    alleDic = getAllMHCIAlleles(methKey, self.getEnumText('specie'))
    return list(alleDic.keys())

  def getSelectedAlleles(self):
    '''Get list of alleles to perform the analysis on'''
    if self.specie.get() != 3 or self.alleleGroup.get() == 3:
      alList = self.alleleCustom.get().strip().split(', ')
    else:
      alList = MHCI_alleles_dic[self.getEnumText('alleleGroup')]
    return alList

  def filterAlleles(self, alDic, alList, lenList):
    '''Filters the allowed alleles and lengths from alList and lenList according to alDic and return the allele and
    length lists necessary to run predict_binding.py'''
    fAL, fLL = [], []
    for allele in alList:
      for length in lenList:
        if str(length) in alDic[allele]:
          fAL.append(allele), fLL.append(str(length))
    return fAL, fLL

  def getMHCOutputFile(self):
    return os.path.abspath(self._getExtraPath('mhc-I_results.tsv'))

  def getRankIdx(self):
    return 8 if self.selType.get() == IC50 else 9

  def parseResults(self, oFile):
    '''Parse the results in the raw_output.tsv file generated by TepiTools and returns a dictionary
    as {seq_id: {(position, epitopeString): {allele: score}}}
    '''
    resAr = self.getResultsArray(oFile)

    # Build the output from the selection
    epiDic = {}
    for row in resAr:
      allele, seq_id, pos, _, _, peptide, _, _, ic50, rank = row[:10]
      key = (int(pos), peptide)
      if not seq_id in epiDic:
        epiDic[seq_id] = {}
      if not key in epiDic[seq_id]:
        epiDic[seq_id][key] = {}

      rankIdx = self.getRankIdx()
      epiDic[seq_id][key][allele] = float(row[rankIdx])
    return epiDic
