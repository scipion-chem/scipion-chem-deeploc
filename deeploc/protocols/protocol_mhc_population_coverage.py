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

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.utils import runInParallel

from .. import Plugin as iedbPlugin
from ..constants import POP_DIC
from ..utils import translateArea, buildMHCCoverageArgs, parseCoverageResults


class ProtMHCPopulationCoverage(EMProtocol):
  """Calculates the population coverage of a series of MHC epitopes"""
  _label = 'mhc population coverage'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequenceROIs', params.PointerParam, pointerClass="SetOfSequenceROIs",
                    label='Input MHC epitopes: ',
                    help="Set of sequence MHC epitopes to calculate the population coverage on")

    iGroup.addParam('mhc', params.EnumParam, label='MHC type: ', display=params.EnumParam.DISPLAY_HLIST,
                    default=2, choices=['I', 'II', 'combined'],
                    help="Whether to predict the coverage on MHC of type I, II or combined.")
    iGroup.addParam('eachROI', params.BooleanParam, label='Perform analysis for each ROI: ', default=True,
                    help="Perform the analysis for each ROI separately")
    iGroup.addParam('allROI', params.BooleanParam, label='Perform analysis for all ROIs together: ', default=True,
                    help="Perform the analysis for all ROIs together, taking into account all their alleles")

    pGroup = form.addGroup('Population')
    pGroup.addParam('pop', params.EnumParam, label='Population by: ',
                    default=0, choices=['Area', 'Ethnicity'],
                    help="Choose population by geographic area or ethnicity")

    pGroup.addParam('area1', params.StringParam, label='Continent: ', default='All', condition='pop==0',
                    help="Choose a continent to add")
    pGroup.addParam('area2', params.StringParam, label='Country: ', default='All',
                    condition='pop==0 and area1!="All" and len(area1.split(",")) == 1',
                    help="Choose a country to add. Use the wizard to choose among the continent countries")
    pGroup.addParam('area3', params.StringParam, label='Country population: ', default='All',
                    condition='pop==0 and area1!="All" and len(area1.split(",")) == 1 and '
                              'area2!="All" and len(area2.split(",")) == 1',
                    help="Choose a country population to add. Use the wizard to choose among the populations in the "
                         "selected country")

    pGroup.addParam('ethnicity', params.StringParam, label='Ethnicity: ', condition='pop==1', default='All',
                    help="Choose an ethnicity to add")

    pGroup.addParam('addArea', params.LabelParam, label='Add defined population: ',
                    help='Add selected area or ethnicity to be included in the analysis')

    sGroup = form.addGroup('Populations summary')
    sGroup.addParam('inAreas', params.TextParam, width=70, default='',
                    label='Populations summary: ',
                    help='Summary of the areas or ethnicities to be included in the analysis')

    form.addParallelSection(threads=4, mpi=1)

  def _insertAllSteps(self):
    self._insertFunctionStep(self.coverageStep)
    self._insertFunctionStep(self.createOutputStep)

  def coverageStep(self):
    selPops = self.getSelectedPopulations()
    epFile, oDir = self._getExtraPath('inputEpitopes.tsv'), self.getResultsDir()
    os.mkdir(oDir)
    self.performCoverageAnalysis(self.inputSequenceROIs.get(), epFile, selPops, self.getEnumText("mhc"), oDir,
                                 self.eachROI, self.allROI)

  def createOutputStep(self):
    coveDic = self.parseResults()
    outROIs = self.inputSequenceROIs.get().createCopy(self._getPath(), copyInfo=True, copyItems=False)

    for roi in self.inputSequenceROIs.get():
      if self.eachROI:
        roiId = str(roi.getObjId())
        if roiId in coveDic:
            roiDic = coveDic[roiId]
            roi._coverage = params.Float(roiDic['coverage'])
            roi._averageHit = params.Float(roiDic['average_hit'])
            roi._pc90 = params.Float(roiDic['pc90'])
        else:
            roi._coverage = params.Float(0)
            roi._averageHit = params.Float(None)
            roi._pc90 = params.Float(None)
      outROIs.append(roi)

    if self.allROI:
      roiDic = coveDic['All']
      outROIs._coverage = params.Float(roiDic['coverage'])
      outROIs._averageHit = params.Float(roiDic['average_hit'])
      outROIs._pc90 = params.Float(roiDic['pc90'])
      outROIs._coverageFile = params.String(self.getCoverageOutputFile())

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)

  ##################### UTILS #####################
  def getResultsDir(self, path=''):
    return self._getPath('coverage_results', path)

  def parseResults(self):
    resDic = {}
    resDir = self.getResultsDir()
    for resFile in os.listdir(resDir):
      outId = resFile.split('_')[-2]
      resFile = os.path.join(resDir, resFile)
      resDic[outId] = parseCoverageResults(resFile)['average']
    return resDic

  def getCoverageOutputFile(self):
    return os.path.abspath(self.getResultsDir('inputEpitopes_All_results.tsv'))

  def performCoverageAnalysis(self, inputROIs, epFile, populations, mhc, oDir, eachROI, allROI):
    nt = self.numberOfThreads.get()
    if eachROI:
      coveArgs = buildMHCCoverageArgs(inputROIs, epFile, populations, mhc, oDir, separated=True)
      runInParallel(iedbPlugin.runPopulationCoverage, None,
                    paramList=[coveArg for coveArg in coveArgs], jobs=nt)

    if allROI:
      coveArgs = buildMHCCoverageArgs(inputROIs, epFile, populations, mhc, oDir, separated=False)
      runInParallel(iedbPlugin.runPopulationCoverage, None,
                    paramList=[coveArg for coveArg in coveArgs], jobs=nt)

  def getAreaOptions(self, level):
    areas = [self.area1.get(), self.area2.get()]
    prevAreas = areas[:level-1]
    cDic = POP_DIC['Area']
    for pArea in prevAreas:
      cDic = cDic[pArea]

    pops = [area for area in cDic]
    pops.sort()
    return pops

  def getPopulationsElement(self):
    newArea = self.getEnumText('pop')
    if newArea == 'Area':
      areaLevel = 1
      cArea = getattr(self, f'area{areaLevel}').get()
      if cArea != 'All':
        newArea += f' / {cArea}'
      while cArea != 'All' and areaLevel < 3 and len(cArea.split(',')) == 1:
        areaLevel += 1
        cArea = getattr(self, f'area{areaLevel}').get()
        if cArea != 'All':
          newArea += f' / {cArea}'

    elif newArea == 'Ethnicity':
      if self.ethnicity.get() != 'All':
        newArea += f' / {self.ethnicity.get()}'
    return newArea

  def getSelectedPopulations(self):
    pops = []
    for popLine in self.inAreas.get().split('\n'):
      popElements = popLine.split('/ ')[-1].split(',')
      [pops.append(p.strip()) for p in popElements if p.strip() not in ['All', '']]

    pops = translateArea(pops)
    return pops

  def _summary(self):
    summ = []
    if os.path.exists(self.getCoverageOutputFile()):
      coveDic = parseCoverageResults(self.getCoverageOutputFile())['average']
      summ.append(f'Average coverage of the set of sequence ROIs (MHC {self.getEnumText("mhc")} epitopes) for the selected populations: '
                  f'\n  - Coverage: {coveDic["coverage"]}'
                  f'\n  - Average hit: {coveDic["average_hit"]}'
                  f'\n  - PC90: {coveDic["pc90"]}'
                  f'\n\nFor more details, check the output file: {self.getCoverageOutputFile()}')
    return summ


