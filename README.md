DelphesNtuplizer
=============

This package allows you to produce to flat Ntuples from Delphes PhaseII samples.

Table of contents
=================
  * [Clone](#clone)
  * [Initialization](#initilization)
  * [Produce Delphes Flat trees](#producing-flatrees)


Clone 
=====

If you do not attempt to contribute to this repository, simply clone it:

```bash
git clone git@github.com:ktht/DelphesNtuplizer.git
```

Initialization
==============

This package requires Delphes to be installed, and CMSSW for gcc, ROOT, FWLite and other dependencies.
Please visit section "Set up your environment" and "Set up Delphes" in https://github.com/ktht/delphes

For setting up this particular repository, run the following:

```bash
cd DelphesNtuplizer # where you cloned this repository
source ./set_env.sh
```

That's it!

Produce Delphes flat trees
==========================

The following command will produce a flat Ntuple, with 10 events, while ignoring branches related to PU (which by default are not).

```bash 
DelphesNtuplizer.py -i delphes_ntuple.root -o flat_tree.root -n 10 -p 0
```
