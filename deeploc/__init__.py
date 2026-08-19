# **************************************************************************
# *
# * Authors:	Blanca Pueche (blanca.pueche@cnb.csic.es)
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
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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
"""
"""

# General imports
import os, subprocess

# Scipion em imports
from pwem import Config as emConfig

# Plugin imports
from pwchem import Plugin as pwchemPlugin
from pwchem.utils import insistentRun
from scipion.install.funcs import InstallHelper

from .bibtex import _bibtexStr
from .constants import *

# Pluging variables
_logo = 'dtu_logo.png'

class Plugin(pwchemPlugin):
    """
    """

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(DEEPLOC_DIC['home'], cls.getDefaultDir(DEEPLOC_DIC))
        cls._defineVar(DEEPLOC_DIC['tar'], None)

    @classmethod
    def defineBinaries(cls, env, default=True):
        if cls.checkVarPath(DEEPLOC_DIC, 'tar'):
            cls.addDeepLocPackage(
                env,
                tarPath=cls.getVar(DEEPLOC_DIC['tar']),
                default=default
            )

        elif cls.checkVarPath(DEEPLOC_DIC, 'home'):
            cls.addDeepLocPackage(
                env,
                deeplocHome=cls.getVar(DEEPLOC_DIC['home']),
                default=default
            )

    @classmethod
    def addDeepLocPackage(cls, env, deeplocHome=None, tarPath=None, default=True):

        DEEPLOC_INSTALLED = '%s_installed' % DEEPLOC_DIC['name']

        emHome = os.path.join(
            emConfig.EM_ROOT,
            cls.getEnvName(DEEPLOC_DIC)
        )

        installationCmd = ''

        if not deeplocHome and tarPath:
            deeplocHome = emHome

            installationCmd += (
                f"mkdir -p {emHome} && "
                f"tar -xf {tarPath} -C {emHome} && "
                f"mv {emHome}/deeploc2_package/* {emHome}/ && "
                f"rm -rf {emHome}/deeploc2_package && "
            )

        elif deeplocHome != emHome:
            installationCmd += (
                f"mkdir -p {emHome} && "
                f"cp -r {deeplocHome}/* {emHome}/ && "
            )

        installer = InstallHelper(
            DEEPLOC_DIC['name'],
            packageHome=emHome,
            packageVersion=DEEPLOC_DIC['version']
        )

        installer.getCondaEnvCommand(
            DEEPLOC_DIC['name'],
            binaryVersion=DEEPLOC_DIC['version'],
            pythonVersion='3.8'
        ).addCommand(
            installationCmd +
            f"{cls.getEnvActivationCommand(DEEPLOC_DIC)} && "
            f"cd {emHome} && "
            "pip install .",
            DEEPLOC_INSTALLED
        ).addPackage(
            env,
            dependencies=['conda'],
            default=default
        )


    @classmethod
    def getDefaultDir(cls, softDic, fn=""):
        emDir = emConfig.EM_ROOT
        for file in os.listdir(emDir):
            if softDic['pattern'] in file.lower():
                foundDir = os.path.join(emDir, file, fn)
                return foundDir.rstrip('/')
        # print(f'BepiPred software could not be found in SOFTWARE directory ({emDir})')
        return os.path.join(emConfig.EM_ROOT, cls.getEnvName(softDic))

    @classmethod
    def checkVarPath(cls, softDic, var='home'):
        '''Check if a plugin variable exists and so do its path'''
        exists = False
        varValue = cls.getVar(softDic[var])
        if varValue and os.path.exists(varValue):
            exists = True
        return exists

    @classmethod
    def checkCallEnv(cls, packageDic):
        actCommand = cls.getVar(packageDic['activation'])
        try:
            if 'conda' in actCommand and not 'shell.bash hook' in actCommand:
                actCommand = f'{cls.getCondaActivationCmd()}{actCommand}'
            subprocess.check_output(actCommand, shell=True)
            envFine = True
        except subprocess.CalledProcessError as e:
            envFine = False
        return envFine

    @classmethod
    def getPluginHome(cls, path=""):
        import deeploc
        fnDir = os.path.split(deeploc.__file__)[0]
        return os.path.join(fnDir, path)

    # ---------------------------------- Protocol functions-----------------------