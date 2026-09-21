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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.widgets import Button

from pyworkflow.protocol import params
from pwem.wizards.wizard import ColorScaleWizardBase
from pwchem.viewers.viewer_structure_attributes import (
    SASAStructureViewer,
    plotSequenceAttribute
)

from ..protocols.protocol_deeploc import ProtDeepLoc


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


def plotLocalizationHistogramFromDataFrame(df):
    localizationCols = [
        col for col in LOCALIZATION_COLUMNS
        if col in df.columns
    ]

    if 'Protein_ID' not in df.columns:
        raise ValueError(
            "The localization file must contain a 'Protein_ID' column."
        )

    if not localizationCols:
        raise ValueError(
            'No localization probability columns were found.'
        )

    proteinIds = df['Protein_ID'].astype(str).tolist()

    nProteins = len(proteinIds)
    nLocalizations = len(localizationCols)

    if nProteins == 0:
        return

    x = np.arange(nLocalizations)
    width = 0.8 / nProteins

    fig, ax = plt.subplots(figsize=(14, 7))

    for i, proteinId in enumerate(proteinIds):
        values = df.iloc[i][localizationCols].astype(float).values

        offset = (i - (nProteins - 1) / 2) * width

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
        localizationCols,
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
                'viewDeepLoc',
                params.LabelParam,
                label='Display DeepLoc residue importance: ',
                help='Display the DeepLoc residue attribute over the sequence.'
            )

        elif hasattr(self.protocol, 'outputAtomStruct'):
            form.addParam(
                'viewLocalization',
                params.LabelParam,
                label='Display localization probabilities: ',
                help='Display the DeepLoc predicted localization probabilities.'
            )

            super()._defineParams(form)

        elif hasattr(self.protocol, 'outputAtomStructs'):
            form.addParam(
                'viewLocalization',
                params.LabelParam,
                label='Display localization probabilities: ',
                help='Display the DeepLoc predicted localization probabilities.'
            )

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
        visDic = {
            'viewLocalization': self._showLocalization
        }

        if hasattr(self.protocol, 'outputSequence'):
            visDic.update({
                'viewSequence': self._showSequenceAttrs,
                'viewDeepLoc': self._showDeepLoc
            })

        elif hasattr(self.protocol, 'outputAtomStruct'):
            visDic.update(
                super()._getVisualizeDict()
            )

        return visDic

    def _getLocalizationObject(self):
        if hasattr(self.protocol, 'outputSequence'):
            return self.protocol.outputSequence

        if hasattr(self.protocol, 'outputAtomStruct'):
            return self.protocol.outputAtomStruct

        if hasattr(self.protocol, 'outputAtomStructs'):
            return self.protocol.outputAtomStructs

        return self.protocol

    def _getLocalizationFile(self, obj):
        localizationPerc = getattr(
            obj,
            '_localizationPerc',
            None
        )

        if localizationPerc is None:
            raise FileNotFoundError(
                'No DeepLoc localization probability file was found.'
            )

        if hasattr(localizationPerc, 'get'):
            localizationPerc = localizationPerc.get()

        if not localizationPerc:
            raise FileNotFoundError(
                'No DeepLoc localization probability file was found.'
            )

        if not os.path.exists(localizationPerc):
            raise FileNotFoundError(
                f'Localization probability file not found: '
                f'{localizationPerc}'
            )

        return localizationPerc

    def _showLocalization(self, e=None):
        obj = self._getLocalizationObject()
        localizationPerc = self._getLocalizationFile(obj)

        plotLocalizationHistogram(localizationPerc)

    def _showDeepLoc(self, paramName=None):
        attrName = self.protocol._ATTRNAME
        attrDic = self.protocol.outputSequence.getAttributesDic()

        plotSequenceAttribute(
            attrDic[attrName],
            attrName=attrName
        )

    def getAtomStruct(self):
        return super().getAtomStruct()


class SequenceAttributeViewer:

    def __init__(self, sequenceData):
        self.sequenceData = sequenceData
        self.index = 0

        self.fig = plt.figure(figsize=(10, 6))
        self.ax = self.fig.add_axes([0.08, 0.20, 0.88, 0.70])

        axPrev = self.fig.add_axes([0.25, 0.05, 0.18, 0.07])
        axNext = self.fig.add_axes([0.57, 0.05, 0.18, 0.07])

        self.prevButton = Button(axPrev, '? Previous')
        self.nextButton = Button(axNext, 'Next ?')

        self.prevButton.on_clicked(self.previous)
        self.nextButton.on_clicked(self.next)

        self.update()

    def update(self):
        data = self.sequenceData[self.index]

        values = data['values']
        seqName = data['name']
        attrName = data['attrName']

        self.ax.clear()

        xs = np.arange(len(values))
        self.ax.bar(xs, values)

        self.ax.set_xlabel('Sequence position')
        self.ax.set_ylabel('{} value'.format(attrName))

        self.ax.set_title(
            '{} ? Sequence {}/{} ? {}'.format(
                attrName,
                self.index + 1,
                len(self.sequenceData),
                seqName
            )
        )

        self.ax.set_xlim(-1, len(values))

        maxY = max(values)
        self.ax.set_ylim(0, maxY + maxY / 10)

        self.fig.canvas.draw_idle()

    def next(self, event):
        if self.index < len(self.sequenceData) - 1:
            self.index += 1
            self.update()

    def previous(self, event):
        if self.index > 0:
            self.index -= 1
            self.update()