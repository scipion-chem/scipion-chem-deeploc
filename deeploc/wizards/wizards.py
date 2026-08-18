# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
"""

from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String

from pwchem.wizards import VariableWizard, AddElementWizard, SelectMultiChainWizard

from ..protocols import *

#################### MHC wizards #####################

class SelectMultiAlleleWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, form):
    return form.protocol.getAvailableAlleles()

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    elementList = self.getListOfElements(form)

    finalList = []
    for i in elementList:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Allele for prediction", provider,
                            "Select one or several alleles")

    vals = [v.get() for v in dlg.values]
    form.setVar(outputParam[0], ', '.join(vals))

SelectMultiAlleleWizard().addTarget(protocol=ProtMHCIPrediction,
                                    targets=['alleleCustom'],
                                    inputs=[],
                                    outputs=['alleleCustom'])

class SelectMultiLengthWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, form):
    return form.protocol.getAvailableLengthList()

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    elementList = self.getListOfElements(form)

    finalList = []
    for i in elementList:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Epitope lengths for prediction", provider,
                            "Select one or several lengths")

    vals = [v.get() for v in dlg.values]
    form.setVar(outputParam[0], ', '.join(vals))

SelectMultiLengthWizard().addTarget(protocol=ProtMHCIPrediction,
                                    targets=['lengths'],
                                    inputs=[],
                                    outputs=['lengths'])

SelectMultiLengthWizard().addTarget(protocol=ProtMHCIIPrediction,
                                    targets=['lengths'],
                                    inputs=[],
                                    outputs=['lengths'])

class SelectAreaWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    from ..constants import POP_DIC
    inputParam, outputParam = self.getInputOutput(form)
    if outputParam[0] == 'ethnicity':
      level = 4
      areaOptions = POP_DIC['Ethnicity']
    elif 'area' in outputParam[0]:
      level = int(outputParam[0][-1])
      areaOptions = form.protocol.getAreaOptions(level)

    finalList = [String('All')]
    for i in areaOptions:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Area options", provider,
                            "Select one or several area options")

    vals = [v.get() for v in dlg.values]
    form.setVar(outputParam[0], ', '.join(vals))

    for nLevel in range(level+1, 4):
      form.setVar(f'area{nLevel}', 'All')


SelectAreaWizard().addTarget(protocol=ProtMHCPopulationCoverage,
                             targets=['ethnicity'],
                             inputs=[],
                             outputs=['ethnicity'])

SelectAreaWizard().addTarget(protocol=ProtMHCPopulationCoverage,
                             targets=['area1'],
                             inputs=[],
                             outputs=['area1'])

SelectAreaWizard().addTarget(protocol=ProtMHCPopulationCoverage,
                             targets=['area2'],
                             inputs=[],
                             outputs=['area2'])

SelectAreaWizard().addTarget(protocol=ProtMHCPopulationCoverage,
                             targets=['area3'],
                             inputs=[],
                             outputs=['area3'])

class AddPopElement(AddElementWizard):
  """Add area or ethnicity to be included in the protocol"""
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    protocol = form.protocol

    newAreas = protocol.getPopulationsElement()
    if newAreas and newAreas.strip() != '':
      prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())
      form.setVar(outputParam[0], prevList + '{}\n'.format(newAreas.strip()))


AddPopElement().addTarget(protocol=ProtMHCPopulationCoverage,
                          targets=['addArea'],
                          inputs=['inAreas'],
                          outputs=['inAreas'])


SelectMultiChainWizard().addTarget(protocol=ProtElliProPrediction,
                                targets=['chain_name'],
                                inputs=['inputAtomStruct'],
                                outputs=['chain_name'])

