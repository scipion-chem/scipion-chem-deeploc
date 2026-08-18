================================
DeepLoc scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Scipion framework plugin for the use of DeepLoc 2.1.

================================
Download DeepLoc files
================================

You can download DeepLoc 2.1 in https://services.healthtech.dtu.dk/services/DeepLoc-2.1/ .

|

Once you obtain the software file (a tar.gz) you have several options to help Scipion finding it:

Option 1) Edit the scipion.conf file and add the variables:
 - DEEPLOC_TAR = <PathToDeepLoc-Tar> (deeploc-2.1.All.tar.gz)

This way, Scipion will untar and move the corresponding files to the scipion/software/em folder and install mhc-i.

Option 2) If you have unzipped the deeploc2.1 tars yourself you can either:

2.1) Move the folder (of the form deeploc2_package) to the scipion/software/em folder. Scipion will find it there.

2.2) Specify the location of the DeepLoc folder in the scipion.conf file as: DEEPLOC_HOME = <PathToDeepLoc_folder>


===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion
to run these protocols. To install the plugin, you have two options:

- **Stable version**

.. code-block::

      scipion installp -p scipion-chem-deeploc

OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

- **Developer's version**

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-deeploc.git

2. **Switch to the desired branch** (main or devel):

Scipion-chem-deeploc is constantly under development.
If you want a relatively older an more stable version, use main branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-deeploc
            git checkout devel

3. **Install**:

.. code-block::

            scipion installp -p path_to_scipion-chem-deeploc --devel

- **Tests**

To check the installation, simply run the following Scipion test:

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_dev.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_prod.svg
