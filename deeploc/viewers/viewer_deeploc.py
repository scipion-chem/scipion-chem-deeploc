# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.widgets import Button

from pyworkflow.protocol import params
from pwem.wizards.wizard import ColorScaleWizardBase

from pwchem.viewers.viewer_structure_attributes import (
    SASAStructureViewer,
    plotSequenceAttribute
)
from pwchem.viewers.viewers_sequences import SequenceAliView

from ..protocols.protocol_deeploc import ProtDeepLoc
from pwem.viewers import ChimeraAttributeViewer


LOCALIZATION_COLUMNS = [
    'Cytoplasm',
    'Nucleus',
    'Extracellular',
    'Cell membrane',
    'Mitochondrion',
    'Plastid',
    'Endoplasmic reticulum',
    'Lysosome/Vacuole',
    'Golgi apparatus',
    'Peroxisome',
    'Peripheral',
    'Transmembrane',
    'Lipid anchor',
    'Soluble'
]


def plotInteractive(data, histogram=False):
    if not data:
        return

    currentIndex = [0]

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.20)

    def drawPlot():
        ax.clear()

        name, attrName, attrValues = data[currentIndex[0]]
        attrValues = list(map(float, attrValues))

        if histogram:
            ax.hist(attrValues, bins=20)
            ax.set_xlabel(attrName)
            ax.set_ylabel('Frequency')
        else:
            xs = np.arange(1, len(attrValues) + 1)
            ax.bar(xs, attrValues)
            ax.set_xlabel('Residue')
            ax.set_ylabel(attrName)
            ax.set_xlim(
                0.5,
                len(attrValues) + 0.5
            )

        ax.set_title(
            '{} - {}'.format(name, attrName)
        )

        fig.canvas.draw_idle()

    previousAx = plt.axes([0.25, 0.05, 0.20, 0.075])
    nextAx = plt.axes([0.55, 0.05, 0.20, 0.075])

    previousButton = Button(previousAx, 'Previous')
    nextButton = Button(nextAx, 'Next')

    def previous(event):
        if currentIndex[0] > 0:
            currentIndex[0] -= 1
            drawPlot()

    def next(event):
        if currentIndex[0] < len(data) - 1:
            currentIndex[0] += 1
            drawPlot()

    previousButton.on_clicked(previous)
    nextButton.on_clicked(next)

    drawPlot()
    plt.show()


def plotSequenceAttributesInteractive(sequenceData):
    plotInteractive(sequenceData)


def plotAtomStructAttributesInteractive(structureData):
    plotInteractive(structureData, histogram=True)


def plotAtomStructSequenceAttributesInteractive(structureData):
    plotInteractive(structureData)


def plotLocalizationHistogramFromDataFrame(df):
    localizationColumns = [
        column for column in LOCALIZATION_COLUMNS
        if column in df.columns
    ]

    if 'Protein_ID' not in df.columns:
        raise ValueError(
            "The localization file must contain a 'Protein_ID' column."
        )

    if not localizationColumns:
        raise ValueError(
            'No DeepLoc localization probability columns were found.'
        )

    proteinIds = df['Protein_ID'].astype(str).tolist()

    nProteins = len(proteinIds)
    nLocalizations = len(localizationColumns)

    if nProteins == 0:
        return

    x = np.arange(nLocalizations)
    width = 0.8 / nProteins

    fig, ax = plt.subplots(figsize=(14, 7))

    for i, proteinId in enumerate(proteinIds):
        values = df.iloc[i][localizationColumns].astype(float).values

        offset = (
            i - (nProteins - 1) / 2
        ) * width

        ax.bar(
            x + offset,
            values,
            width,
            label=proteinId
        )

    ax.set_xlabel('Localization')
    ax.set_ylabel('Probability')
    ax.set_title('DeepLoc localization probabilities')

    ax.set_xticks(x)
    ax.set_xticklabels(
        localizationColumns,
        rotation=45,
        ha='right'
    )

    ax.set_ylim(0, 1)
    ax.legend(title='Protein ID')

    plt.tight_layout()
    plt.show()


def plotLocalizationHistogram(csvFile):
    df = pd.read_csv(csvFile)
    plotLocalizationHistogramFromDataFrame(df)


class DeepLocStructureViewer(SASAStructureViewer):

    _targets = [ProtDeepLoc]
    _label = 'DeepLoc viewer'

    def _defineParams(self, form):

        # --------------------------------------------------------------
        # Single sequence
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputSequence'):
            form.addSection(
                label='Visualization of DeepLoc residue importance'
            )

            form.addParam(
                'viewLocalization',
                params.LabelParam,
                label='Display localization probabilities: ',
                help='Display the DeepLoc predicted localization probabilities.'
            )

            form.addParam(
                'viewSequence',
                params.LabelParam,
                label='View sequence: ',
                help='View output sequence'
            )

            form.addParam(
                'viewSequenceAttribute',
                params.LabelParam,
                label='Display DeepLoc residue importance: ',
                help='Display the residue attribute over the sequence.'
            )

        # --------------------------------------------------------------
        # Set of sequences
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputSequences'):
            form.addSection(
                label='Visualization of DeepLoc residue importance'
            )

            form.addParam(
                'viewLocalization',
                params.LabelParam,
                label='Display localization probabilities: ',
                help='Display the DeepLoc predicted localization probabilities.'
            )

            form.addParam(
                'viewSequences',
                params.LabelParam,
                label='View all sequences: ',
                help='View all output sequences.'
            )

            form.addParam(
                'viewSequencesAttribute',
                params.LabelParam,
                label='Display residue attribute for all sequences: ',
                help='Display the residue attribute for all output sequences.'
            )

        # --------------------------------------------------------------
        # Single AtomStruct
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputAtomStruct'):
            form.addParam(
                'viewLocalization',
                params.LabelParam,
                label='Display localization probabilities: ',
                help='Display the DeepLoc predicted localization probabilities.'
            )

            ChimeraAttributeViewer._defineParams(self, form)
            self._defineColorScale(form)

        # --------------------------------------------------------------
        # Set of AtomStructs
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputAtomStructs'):
            form.addParam(
                'viewLocalization',
                params.LabelParam,
                label='Display localization probabilities: ',
                help='Display the DeepLoc predicted localization probabilities.'
            )

            self._defineAtomStructSetParams(form)

    def _defineColorScale(self, form):
        group = form.addGroup('Color settings')

        ColorScaleWizardBase.defineColorScaleParams(
            group,
            defaultLowest=0,
            defaultHighest=1,
            defaultIntervals=21,
            defaultColorMap='RdBu'
        )

    def _getVisualizeDict(self):
        visDic = {}

        # --------------------------------------------------------------
        # Localization
        # --------------------------------------------------------------

        if (
            hasattr(self.protocol, 'outputSequence')
            or hasattr(self.protocol, 'outputSequences')
            or hasattr(self.protocol, 'outputAtomStruct')
            or hasattr(self.protocol, 'outputAtomStructs')
        ):
            visDic['viewLocalization'] = self._showLocalization

        # --------------------------------------------------------------
        # Single sequence
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputSequence'):
            visDic.update({
                'viewSequence': self._showSequenceAttrs,
                'viewSequenceAttribute': self._showSequenceAttribute
            })

        # --------------------------------------------------------------
        # Set of sequences
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputSequences'):
            visDic.update({
                'viewSequences': self._showSequences,
                'viewSequencesAttribute':
                    self._showSequencesAttribute
            })

        # --------------------------------------------------------------
        # Single AtomStruct
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputAtomStruct'):
            visDic.update(
                ChimeraAttributeViewer._getVisualizeDict(self)
            )

        # --------------------------------------------------------------
        # Set of AtomStructs
        # --------------------------------------------------------------

        if hasattr(self.protocol, 'outputAtomStructs'):
            visDic.update({
                'viewAtomStruct': self._showAtomStruct,
                'viewAtomStructAttribute':
                    self._showAtomStructAttribute,
                'viewAtomStructSequenceAttribute':
                    self._showAtomStructSequenceAttribute
            })

        return visDic

    # ------------------------------------------------------------------
    # Localization
    # ------------------------------------------------------------------

    def _getLocalizationObject(self):
        if hasattr(self.protocol, 'outputSequence'):
            return self.getOutSequences()

        if hasattr(self.protocol, 'outputSequences'):
            return self.protocol.outputSequences

        if hasattr(self.protocol, 'outputAtomStruct'):
            return self.getAtomStructObject()

        if hasattr(self.protocol, 'outputAtomStructs'):
            return self.protocol.outputAtomStructs

        return self.protocol

    def _getLocalizationFile(self):
        obj = self._getLocalizationObject()

        localizationFile = getattr(
            obj,
            '_localizationPerc',
            None
        )

        if hasattr(localizationFile, 'get'):
            localizationFile = localizationFile.get()

        if not localizationFile:
            raise FileNotFoundError(
                'No DeepLoc localization probability file was found.'
            )

        if not os.path.exists(localizationFile):
            raise FileNotFoundError(
                'DeepLoc localization probability file not found: '
                '{}'.format(localizationFile)
            )

        return localizationFile

    def _showLocalization(self, paramName=None):
        localizationFile = self._getLocalizationFile()

        plotLocalizationHistogram(localizationFile)

    # ------------------------------------------------------------------
    # Single sequence
    # ------------------------------------------------------------------

    def _showSequenceAttrs(self, paramName=None):
        obj = self.protocol.outputSequence

        outPath = os.path.abspath(
            self.protocol._getExtraPath(
                'viewSequences_{}.fasta'.format(obj.getId())
            )
        )

        obj.exportToFile(outPath)

        return [
            SequenceAliView(
                [outPath],
                cwd=self.protocol._getExtraPath()
            )
        ]

    def _showSequenceAttribute(self, paramName=None):
        attrName = self.protocol._ATTRNAME
        attrDic = self.protocol.outputSequence.getAttributesDic()

        plotSequenceAttribute(
            attrDic[attrName],
            attrName=attrName
        )

    # ------------------------------------------------------------------
    # Set of sequences
    # ------------------------------------------------------------------

    def _showSequences(self, paramName=None):
        outPath = os.path.abspath(
            self.protocol._getExtraPath('viewSequences.fasta')
        )

        if os.path.exists(outPath):
            os.remove(outPath)

        self.protocol.outputSequences.exportToFile(outPath)

        return [
            SequenceAliView(
                [outPath],
                cwd=self.protocol._getExtraPath()
            )
        ]

    def _showSequencesAttribute(self, paramName=None):
        attrName = self.protocol._ATTRNAME

        sequenceData = [
            (
                sequence.getSeqName(),
                attrName,
                sequence.getAttributesDic()[attrName]
            )
            for sequence in self.protocol.outputSequences
            if attrName in sequence.getAttributesDic()
        ]

        plotSequenceAttributesInteractive(sequenceData)

    # ------------------------------------------------------------------
    # Set of AtomStructs
    # ------------------------------------------------------------------

    def _defineAtomStructSetParams(self, form):
        form.addSection(label='Visualization of structure set')

        self._structureNames = [
            os.path.splitext(
                os.path.basename(atomStruct.getFileName())
            )[0]
            for atomStruct in self.protocol.outputAtomStructs
        ]

        form.addParam(
            'atomStruct',
            params.EnumParam,
            choices=self._structureNames,
            default=0,
            label='Structure: ',
            help='Select the structure to visualize.'
        )

        form.addParam(
            'viewAtomStruct',
            params.LabelParam,
            label='View structure: ',
            help='View the selected structure in ChimeraX.'
        )

        form.addParam(
            'viewAtomStructAttribute',
            params.LabelParam,
            label='Display attribute histogram: ',
            help='Display the residue attribute distribution '
                 'for the selected structure.'
        )

        form.addParam(
            'chain_name',
            params.StringParam,
            default='A',
            allowsNull=True,
            label='Chain of interest: ',
            help='Specify the chain of interest (e.g. A).'
        )

        form.addParam(
            'viewAtomStructSequenceAttribute',
            params.LabelParam,
            label='Display attribute over sequence: ',
            help='Display the residue attribute over the sequence '
                 'for the selected structure.'
        )

        self._defineColorScale(form)

    def _getSelectedAtomStruct(self):
        selectedName = self._structureNames[self.atomStruct.get()]

        for atomStruct in self.protocol.outputAtomStructs:
            name = os.path.splitext(
                os.path.basename(atomStruct.getFileName())
            )[0]

            if name == selectedName:
                return atomStruct

        return None

    def _showAtomStruct(self, paramName=None):
        self._atomStruct = self._getSelectedAtomStruct()

        intermediateFile = self.protocol._getExtraPath(
            'chimeraAttribute_{}.cif'.format(
                self.protocol._ATTRNAME
            )
        )

        if os.path.exists(intermediateFile):
            os.remove(intermediateFile)

        return self._showChimera(paramName)

    def _showAtomStructAttribute(self, paramName=None):
        self._atomStruct = self._getSelectedAtomStruct()

        return self._showHistogram(paramName)

    def _showAtomStructSequenceAttribute(self, paramName=None):
        self._atomStruct = self._getSelectedAtomStruct()

        return self._showSequence(paramName)

    def getAtomStructObject(self):
        if hasattr(self.protocol, 'outputAtomStructs'):
            return self._getSelectedAtomStruct()

        return super().getAtomStructObject()

    def getEnumText(self, paramName):
        if (
            paramName == 'attrName'
            and hasattr(self.protocol, 'outputAtomStructs')
        ):
            return self.protocol._ATTRNAME

        return super().getEnumText(paramName)

    def getOutSequences(self):
        if hasattr(self.protocol, 'outputSequence'):
            return self.protocol.outputSequence

        if hasattr(self.protocol, 'outputSequences'):
            return self.protocol.outputSequences

        return None