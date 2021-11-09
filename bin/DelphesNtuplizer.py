#!/usr/bin/env python

import sys
argv = sys.argv
sys.argv = []
import ROOT
ROOT.gSystem.ResetSignals()
sys.argv = argv

from array import array
import argparse


class TreeProducer:
    def __init__(self, pileup, prune, debug):

         # flat tree branches
         self.debug = debug
         self.prune = prune

         self.t = ROOT.TTree( "Events","Events" )
         self.maxn = 9999

         if pileup and not self.prune:
             self.vtx_size         = array( 'i', [ 0 ] )
             self.vtx_ndf          = array( 'i', self.maxn*[0] )
             self.vtx_pt2          = array( 'f', self.maxn*[ 0. ] )
             self.vtx_x            = array( 'f', self.maxn*[ 0. ] )
             self.vtx_y            = array( 'f', self.maxn*[ 0. ] )
             self.vtx_z            = array( 'f', self.maxn*[ 0. ] )

         ## main MC weight
         self.evt_nr            = array( 'l', [ 0 ] )
         self.genweight         = array( 'f', [ 0. ] )
         self.scale             = array( 'f', [ 0. ] )
         if not self.prune:
             self.scalePDF          = array( 'f', [ 0. ] )
             self.pdf1              = array( 'f', [ 0. ] )
             self.pdf2              = array( 'f', [ 0. ] )
         self.x1                = array( 'f', [ 0. ] )
         self.x2                = array( 'f', [ 0. ] )
         self.id1               = array( 'i', [ 0 ] )
         self.id2               = array( 'i', [ 0 ] )
         self.alphaQED          = array( 'f', [ 0. ] )
         self.alphaQCD          = array( 'f', [ 0. ] )
         if not self.prune:
            self.pid               = array( 'i', [ 0 ] )
         
         ## variation event weights (from "LHEEventProduct")
         if not self.prune:
             self.lheweight_size   = array( 'i', [ 0 ] )
             self.lheweight_val    = array( 'f', self.maxn*[ 0 ] )

         self.genpart_size     = array( 'i', [ 0 ] )
         self.genpart_pid      = array( 'i', self.maxn*[ 0 ] )
         self.genpart_status   = array( 'i', self.maxn*[ 0 ] )
         self.genpart_pt       = array( 'f', self.maxn*[ 0. ] )
         self.genpart_eta      = array( 'f', self.maxn*[ 0. ] )
         self.genpart_phi      = array( 'f', self.maxn*[ 0. ] )
         self.genpart_mass     = array( 'f', self.maxn*[ 0. ] )
         self.genpart_m1       = array( 'i', self.maxn*[ 0 ] )
         #self.genpart_m2       = array( 'i', self.maxn*[ 0 ] )
         self.genpart_d1       = array( 'i', self.maxn*[ 0 ] )
         self.genpart_d2       = array( 'i', self.maxn*[ 0 ] )

         self.genjet_size      = array( 'i', [ 0 ] )
         self.genjet_pt        = array( 'f', self.maxn*[ 0. ] )
         self.genjet_eta       = array( 'f', self.maxn*[ 0. ] )
         self.genjet_phi       = array( 'f', self.maxn*[ 0. ] )
         self.genjet_mass      = array( 'f', self.maxn*[ 0. ] )
         if not self.prune:
             self.genjet_deltaEta  = array( 'f', self.maxn*[ 0. ] )
             self.genjet_deltaPhi  = array( 'f', self.maxn*[ 0. ] )
             self.genjet_charge    = array( 'i', self.maxn*[ 0 ] )
             self.genjet_nCharged  = array( 'i', self.maxn*[ 0 ] )
             self.genjet_nNeutral  = array( 'i', self.maxn*[ 0 ] )
             self.genjet_nef       = array( 'f', self.maxn*[ 0. ] )
             self.genjet_cef       = array( 'f', self.maxn*[ 0. ] )

         self.genmet_pt        = array( 'f', [ 0. ] )
         self.genmet_phi       = array( 'f', [ 0. ] )

         if not self.prune:
             self.gamma_size       = array( 'i', [ 0 ] )
             self.gamma_pt         = array( 'f', self.maxn*[ 0. ] )
             self.gamma_eta        = array( 'f', self.maxn*[ 0. ] )
             self.gamma_phi        = array( 'f', self.maxn*[ 0. ] )
             self.gamma_mass       = array( 'f', self.maxn*[ 0. ] )
             self.gamma_reliso     = array( 'f', self.maxn*[ 0. ] )
             self.gamma_relisoRC   = array( 'f', self.maxn*[ 0. ] )
             self.gamma_sumPtCh    = array( 'f', self.maxn*[ 0. ] )
             self.gamma_sumPtNeu   = array( 'f', self.maxn*[ 0. ] )
             self.gamma_sumPtCPU   = array( 'f', self.maxn*[ 0. ] )
             self.gamma_sumPt      = array( 'f', self.maxn*[ 0. ] )

         self.elec_size        = array( 'i', [ 0 ] )
         self.elec_pt          = array( 'f', self.maxn*[ 0. ] )
         self.elec_eta         = array( 'f', self.maxn*[ 0. ] )
         self.elec_phi         = array( 'f', self.maxn*[ 0. ] )
         self.elec_mass        = array( 'f', self.maxn*[ 0. ] )
         self.elec_charge      = array( 'i', self.maxn*[ 0 ] )
         self.elec_dz          = array( 'f', self.maxn*[ 0. ] )
         self.elec_reliso      = array( 'f', self.maxn*[ 0. ] )
         self.elec_relisoRC    = array( 'f', self.maxn*[ 0. ] )
         self.elec_sumPtCh     = array( 'f', self.maxn*[ 0. ] )
         self.elec_sumPtNeu    = array( 'f', self.maxn*[ 0. ] )
         self.elec_sumPtCPU    = array( 'f', self.maxn*[ 0. ] )
         self.elec_sumPt       = array( 'f', self.maxn*[ 0. ] )

         self.muon_size        = array( 'i', [ 0 ] )
         self.muon_pt          = array( 'f', self.maxn*[ 0. ] )
         self.muon_eta         = array( 'f', self.maxn*[ 0. ] )
         self.muon_phi         = array( 'f', self.maxn*[ 0. ] )
         self.muon_mass        = array( 'f', self.maxn*[ 0. ] )
         self.muon_charge      = array( 'i', self.maxn*[ 0 ] )
         self.muon_dz          = array( 'f', self.maxn*[ 0. ] )
         self.muon_reliso      = array( 'f', self.maxn*[ 0. ] )
         self.muon_relisoRC    = array( 'f', self.maxn*[ 0. ] )
         self.muon_sumPtCh     = array( 'f', self.maxn*[ 0. ] )
         self.muon_sumPtNeu    = array( 'f', self.maxn*[ 0. ] )
         self.muon_sumPtCPU    = array( 'f', self.maxn*[ 0. ] )
         self.muon_sumPt       = array( 'f', self.maxn*[ 0. ] )

         self.jet_size         = array( 'i', [ 0 ] )
         self.jet_pt           = array( 'f', self.maxn*[ 0. ] )
         self.jet_eta          = array( 'f', self.maxn*[ 0. ] )
         self.jet_phi          = array( 'f', self.maxn*[ 0. ] )
         self.jet_mass         = array( 'f', self.maxn*[ 0. ] )
         if not self.prune:
             self.jet_deltaEta     = array( 'f', self.maxn*[ 0. ] )
             self.jet_deltaPhi     = array( 'f', self.maxn*[ 0. ] )
             self.jet_charge       = array( 'i', self.maxn*[ 0 ] )
             self.jet_hoe          = array( 'f', self.maxn*[ 0. ] )
             self.jet_nCharged     = array( 'i', self.maxn*[ 0 ] )
             self.jet_nNeutral     = array( 'i', self.maxn*[ 0 ] )
             self.jet_nef          = array( 'f', self.maxn*[ 0. ] )
             self.jet_cef          = array( 'f', self.maxn*[ 0. ] )
         self.jet_btag         = array( 'i', self.maxn*[ 0 ] )
         self.jet_btagAlgo     = array( 'i', self.maxn*[ 0 ] )
         self.jet_btagPhys     = array( 'i', self.maxn*[ 0 ] )
         self.jet_flavor       = array( 'i', self.maxn*[ 0 ] )
         self.jet_flavorAlgo   = array( 'i', self.maxn*[ 0 ] )
         self.jet_flavorPhys   = array( 'i', self.maxn*[ 0 ] )

         if pileup and not self.prune:
             self.rho_size         = array( 'i', [ 0 ] )
             self.rho              = array( 'f', self.maxn*[ 0. ] )

         self.met_pt           = array( 'f', [ 0. ] )
         self.met_phi          = array( 'f', [ 0. ] )

         # declare tree branches
         self.t.Branch( "event", self.evt_nr, "event/L")
         self.t.Branch( "scale", self.scale, "scale/F")
         if not self.prune:
             self.t.Branch( "scalePDF", self.scalePDF, "scalePDF/F")
             self.t.Branch( "pdf1", self.pdf1, "pdf1/F")
             self.t.Branch( "pdf2", self.pdf2, "pdf2/F")
         self.t.Branch( "x1", self.x1, "x1/F")
         self.t.Branch( "x2", self.x2, "x2/F")
         self.t.Branch( "id1", self.id1, "id1/I")
         self.t.Branch( "id2", self.id2, "id2/I")
         self.t.Branch( "alphaQED", self.alphaQED, "alphaQED/F")
         self.t.Branch( "alphaQCD", self.alphaQCD, "alphaQCD/F")
         if not self.prune:
            self.t.Branch( "PID", self.pid, "PID/I")
         
         self.t.Branch( "genweight", self.genweight, "genweight/F")
         if not self.prune:
             self.t.Branch( "nLHEWeight", self.lheweight_size, "nLHEWeight/I")
             self.t.Branch( "LHEWeight", self. lheweight_val, "LHEWeight[nLHEWeight]/I")

         if pileup and not self.prune:
             self.t.Branch( "nVertex", self.vtx_size, "nVertex/I")
             self.t.Branch( "Vertex_ndf", self.vtx_ndf, "Vertex_ndf[nVertex]/I")
             self.t.Branch( "Vertex_pt2", self.vtx_pt2, "Vertex_pt2[nVertex]/F")
             self.t.Branch( "Vertex_x", self.vtx_x, "Vertex_x[nVertex]/F")
             self.t.Branch( "Vertex_y", self.vtx_y, "Vertex_y[nVertex]/F")
             self.t.Branch( "Vertex_z", self.vtx_z, "Vertex_z[nVertex]/F")

         self.t.Branch( "nGenPart", self.genpart_size, "nGenPart/I")
         self.t.Branch( "GenPart_pdgId", self.genpart_pid, "GenPart_pdgId[nGenPart]/I")
         self.t.Branch( "GenPart_status", self.genpart_status, "GenPart_status[nGenPart]/I")
         self.t.Branch( "GenPart_pt", self.genpart_pt, "GenPart_pt[nGenPart]/F")
         self.t.Branch( "GenPart_eta", self.genpart_eta, "GenPart_eta[nGenPart]/F")
         self.t.Branch( "GenPart_phi", self.genpart_phi, "GenPart_phi[nGenPart]/F")
         self.t.Branch( "GenPart_mass", self.genpart_mass, "GenPart_mass[nGenPart]/F")
         self.t.Branch( "GenPart_genPartIdxMother", self.genpart_m1, "genpart_genPartIdxMother[nGenPart]/I")
         #self.t.Branch( "GenPart_m2", self.genpart_m2, "GenPart_m2[nGenPart]/I")
         self.t.Branch( "GenPart_d1", self.genpart_d1, "GenPart_d1[nGenPart]/I")
         self.t.Branch( "GenPart_d2", self.genpart_d2, "GenPart_d2[nGenPart]/I")

         self.t.Branch( "nGenJet", self.genjet_size, "nGenJet/I")
         self.t.Branch( "GenJet_pt", self.genjet_pt, "GenJet_pt[nGenJet]/F")
         self.t.Branch( "GenJet_eta", self.genjet_eta, "GenJet_eta[nGenJet]/F")
         self.t.Branch( "GenJet_phi", self.genjet_phi, "GenJet_phi[nGenJet]/F")
         self.t.Branch( "GenJet_mass", self.genjet_mass, "GenJet_mass[nGenJet]/F")
         if not self.prune:
             self.t.Branch( "GenJet_deltaEta", self.genjet_deltaEta, "GenJet_deltaEta[nGenJet]/F")
             self.t.Branch( "GenJet_deltaPhi", self.genjet_deltaPhi, "GenJet_deltaPhi[nGenJet]/F")
             self.t.Branch( "GenJet_charge", self.genjet_charge, "GenJet_charge[nGenJet]/I")
             self.t.Branch( "GenJet_nCharged", self.genjet_nCharged, "GenJet_nCharged[nGenJet]/I")
             self.t.Branch( "GenJet_nNeutral", self.genjet_nNeutral, "GenJet_nNeutral[nGenJet]/I")
             self.t.Branch( "GenJet_nef", self.genjet_nef, "GenJet_nef[nGenJet]/F")
             self.t.Branch( "GenJet_cef", self.genjet_cef, "GenJet_cef[nGenJet]/F")

         self.t.Branch( "GenMET_pt", self.genmet_pt, "GenMET_pt/F")
         self.t.Branch( "GenMET_phi", self.genmet_phi, "GenMET_phi/F")

         if not self.prune:
             self.t.Branch( "nPhoton", self.gamma_size, "nPhoton/I")
             self.t.Branch( "Photon_pt", self.gamma_pt, "Photon_pt[nPhoton]/F")
             self.t.Branch( "Photon_eta", self.gamma_eta, "Photon_eta[nPhoton]/F")
             self.t.Branch( "Photon_phi", self.gamma_phi, "Photon_phi[nPhoton]/F")
             self.t.Branch( "Photon_mass", self.gamma_mass, "Photon_mass[nPhoton]/F")
             self.t.Branch( "Photon_relIso", self.gamma_reliso, "Photon_relIso[nPhoton]/F")
             self.t.Branch( "Photon_relIsoRhoCorr", self.gamma_relisoRC, "Photon_relIsoRhoCorr[nPhoton]/F")
             self.t.Branch( "Photon_sumPt", self.gamma_sumPt, "Photon_sumPt[nPhoton]/F")
             self.t.Branch( "Photon_sumPtCh", self.gamma_sumPtCh, "Photon_sumPtCh[nPhoton]/F")
             self.t.Branch( "Photon_sumPtNeu", self.gamma_sumPtNeu, "Photon_sumPtNeu[nPhoton]/F")
             self.t.Branch( "Photon_sumPtCPU", self.gamma_sumPtCPU, "Photon_sumPtCPU[nPhoton]/F")

         self.t.Branch( "nElectron", self.elec_size, "nElectron/I")
         self.t.Branch( "Electron_pt", self.elec_pt, "Electron_pt[nElectron]/F")
         self.t.Branch( "Electron_eta", self.elec_eta, "Electron_eta[nElectron]/F")
         self.t.Branch( "Electron_phi", self.elec_phi, "Electron_phi[nElectron]/F")
         self.t.Branch( "Electron_mass", self.elec_mass, "Electron_mass[nElectron]/F")
         self.t.Branch( "Electron_charge", self.elec_charge, "Electron_charge[nElectron]/I")
         self.t.Branch( "Electron_dz", self.elec_dz, "Electron_dz[nElectron]/F")
         self.t.Branch( "Electron_relIso", self.elec_reliso, "Electron_relIso[nElectron]/F")
         self.t.Branch( "Electron_relIsoRhoCorr", self.elec_relisoRC, "Electron_relIsoRhoCorr[nElectron]/F")
         self.t.Branch( "Electron_sumPt", self.elec_sumPt, "Electron_sumPt[nElectron]/F")
         self.t.Branch( "Electron_sumPtCh", self.elec_sumPtCh, "Electron_sumPtCh[nElectron]/F")
         self.t.Branch( "Electron_sumPtNeu", self.elec_sumPtNeu, "Electron_sumPtNeu[nElectron]/F")
         self.t.Branch( "Electron_sumPtCPU", self.elec_sumPtCPU, "Electron_sumPtCPU[nElectron]/F")

         self.t.Branch( "nMuon", self.muon_size, "nMuon/I")
         self.t.Branch( "Muon_pt", self.muon_pt, "Muon_pt[nMuon]/F")
         self.t.Branch( "Muon_eta", self.muon_eta, "Muon_eta[nMuon]/F")
         self.t.Branch( "Muon_phi", self.muon_phi, "Muon_phi[nMuon]/F")
         self.t.Branch( "Muon_mass", self.muon_mass, "Muon_mass[nMuon]/F")
         self.t.Branch( "Muon_charge", self.muon_charge, "Muon_charge[nMuon]/I")
         self.t.Branch( "Muon_dz", self.muon_dz, "Muon_dz[nMuon]/F")
         self.t.Branch( "Muon_relIso", self.muon_reliso, "Muon_relIso[nMuon]/F")
         self.t.Branch( "Muon_relIsoRhoCorr", self.muon_relisoRC, "Muon_relIsoRhoCorr[nMuon]/F")
         self.t.Branch( "Muon_sumPt", self.muon_sumPt, "Muon_sumPt[nMuon]/F")
         self.t.Branch( "Muon_sumPtCh", self.muon_sumPtCh, "Muon_sumPtCh[nMuon]/F")
         self.t.Branch( "Muon_sumPtNeu", self.muon_sumPtNeu, "Muon_sumPtNeu[nMuon]/F")
         self.t.Branch( "Muon_sumPtCPU", self.muon_sumPtCPU, "Muon_sumPtCPU[nMuon]/F")

         self.t.Branch( "nJet", self.jet_size, "nJet/I")
         self.t.Branch( "Jet_pt", self.jet_pt, "Jet_pt[nJet]/F")
         self.t.Branch( "Jet_eta", self.jet_eta, "Jet_eta[nJet]/F")
         self.t.Branch( "Jet_phi", self.jet_phi, "Jet_phi[nJet]/F")
         self.t.Branch( "Jet_mass", self.jet_mass, "Jet_mass[nJet]/F")
         if not self.prune:
             self.t.Branch( "Jet_deltaEta", self.jet_deltaEta, "Jet_deltaEta[nJet]/F")
             self.t.Branch( "Jet_deltaPhi", self.jet_deltaPhi, "Jet_deltaPhi[nJet]/F")
             self.t.Branch( "Jet_charge", self.jet_charge, "Jet_charge[nJet]/I")
             self.t.Branch( "Jet_hoe", self.jet_hoe, "Jet_hoe[nJet]/F")
             self.t.Branch( "Jet_nCharged", self.jet_nCharged, "Jet_nCharged[nJet]/I")
             self.t.Branch( "Jet_nNeutral", self.jet_nNeutral, "Jet_nNeutral[nJet]/I")
             self.t.Branch( "Jet_nef", self.jet_nef, "Jet_nef[nJet]/F")
             self.t.Branch( "Jet_cef", self.jet_cef, "Jet_cef[nJet]/F")
         self.t.Branch( "Jet_btag", self.jet_btag,"Jet_btag[nJet]/I")
         self.t.Branch( "Jet_btagAlgo", self.jet_btagAlgo, "Jet_btagAlgo[nJet]/I")
         self.t.Branch( "Jet_btagPhys", self.jet_btagPhys, "Jet_btagPhys[nJet]/I")
         self.t.Branch( "Jet_flavor", self.jet_flavor,"Jet_flavor[nJet]/I")
         self.t.Branch( "Jet_flavorAlgo", self.jet_flavorAlgo, "Jet_flavorAlgo[nJet]/I")
         self.t.Branch( "Jet_flavorPhys", self.jet_flavorPhys, "Jet_flavorPhys[nJet]/I")

         if pileup and not self.prune:
             self.t.Branch( "nRho", self.rho_size, "nRho/I")
             self.t.Branch( "Rho", self.rho, "Rho[nRho]/F")

         self.t.Branch( "MET_pt", self. met_pt, "MET_pt/F")
         self.t.Branch( "MET_phi", self.met_phi, "MET_phi/F")

    #___________________________________________
    def processVertices(self, vertices):
        i = 0
        for item in vertices:
            self.vtx_pt2[i] = item.SumPT2
            self.vtx_ndf[i] = item.NDF
            self.vtx_x[i] = item.X   #Gamze
            self.vtx_y[i] = item.Y   #Gamze
            self.vtx_z[i] = item.Z   #Gamze
            i += 1
        self.vtx_size[0] = i


    #___________________________________________
    def processEvent(self, event, weights):
        i = 0

        self.genweight[0] = event[0].Weight
        self.evt_nr[0]    = event[0].Number
        self.scale[0]     = event[0].Scale
        if not self.prune:
            self.scalePDF[0]  = event[0].ScalePDF
            self.pdf1[0]      = event[0].PDF1
            self.pdf2[0]      = event[0].PDF2
        self.x1[0]        = event[0].X1
        self.x2[0]        = event[0].X2
        self.id1[0]       = event[0].ID1
        self.id2[0]       = event[0].ID2
        self.alphaQED[0]  = event[0].AlphaQED
        self.alphaQCD[0]  = event[0].AlphaQCD
        if not self.prune:
            self.pid[0]       = event[0].ProcessID

        if not self.prune:
            for item in weights:
                self.lheweight_val[i] = item.Weight
                i += 1
            self.lheweight_size[0] = i


    #___________________________________________
    def shifIdx(self, idx, eliminate):
        if idx < 0 or not eliminate:
            return idx
        shift = len([ i for i in eliminate if i < idx ])
        assert(shift <= idx)
        return idx - shift

    def processGenParticles(self, particles):
        eliminate = []
        if self.prune:
            for idx, item in enumerate(particles):
                if item.D1 == -1 and item.D2 == -1 and ( # final state
                    (item.M1 >= 0 and particles[item.M1].M1 < 0) or # straight from proton
                    (item.Status != 1) or # unstable
                    (abs(item.PID) not in [ 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16 ]) # keep only leptons and quarks
                ):
                    eliminate.append(idx)

        i = 0
        for j, item in enumerate(particles):
            if j not in eliminate:
                self.genpart_pid    [i] = item.PID
                self.genpart_status [i] = item.Status
                self.genpart_pt     [i] = item.PT
                self.genpart_eta    [i] = item.Eta
                self.genpart_phi    [i] = item.Phi
                self.genpart_mass   [i] = item.Mass
                self.genpart_m1     [i] = self.shifIdx(item.M1, eliminate)
                #self.genpart_m2     [i] = self.shifIdx(item.M2, eliminate)
                self.genpart_d1     [i] = self.shifIdx(item.D1, eliminate)
                self.genpart_d2     [i] = self.shifIdx(item.D2, eliminate)
                i += 1
        self.genpart_size[0] = i

    #___________________________________________
    def processGenJets(self, genjets):
        i = 0
        for item in genjets:
            self.genjet_pt       [i] = item.PT
            self.genjet_eta      [i] = item.Eta
            self.genjet_phi      [i] = item.Phi
            self.genjet_mass     [i] = item.Mass
            if not self.prune:
                self.genjet_deltaEta [i] = item.DeltaEta
                self.genjet_deltaPhi [i] = item.DeltaPhi
                self.genjet_charge   [i] = item.Charge
                self.genjet_nCharged [i] = item.NCharged
                self.genjet_nNeutral [i] = item.NNeutrals
                self.genjet_nef      [i] = item.NeutralEnergyFraction
                self.genjet_cef      [i] = item.ChargedEnergyFraction
            i += 1
        self.genjet_size[0] = i

    #___________________________________________
    def processGenMissingET(self, met):
        if met:
            self.genmet_pt    [0] = met[0].MET
            self.genmet_phi   [0] = met[0].Phi

    #___________________________________________
    def processPhotons(self, photons):
        i = 0
        for item in photons:
            self.gamma_pt      [i] = item.PT
            self.gamma_eta     [i] = item.Eta
            self.gamma_phi     [i] = item.Phi
            self.gamma_mass    [i] = item.P4().M()
            self.gamma_reliso  [i] = item.IsolationVar
            self.gamma_relisoRC[i] = item.IsolationVarRhoCorr
            self.gamma_sumPt   [i] = item.SumPt
            self.gamma_sumPtCh [i] = item.SumPtCharged
            self.gamma_sumPtNeu[i] = item.SumPtNeutral
            self.gamma_sumPtCPU[i] = item.SumPtChargedPU

            i += 1

        self.gamma_size[0] = i


    #___________________________________________
    def processElectrons(self, electrons):
        i = 0
        for item in electrons:
            self.elec_pt      [i] = item.PT
            self.elec_eta     [i] = item.Eta
            self.elec_phi     [i] = item.Phi
            self.elec_mass    [i] = item.P4().M()
            self.elec_charge  [i] = item.Charge
            self.elec_dz      [i] = item.DZ
            self.elec_reliso  [i] = item.IsolationVar
            self.elec_relisoRC[i] = item.IsolationVarRhoCorr
            self.elec_sumPt   [i] = item.SumPt
            self.elec_sumPtCh [i] = item.SumPtCharged
            self.elec_sumPtNeu[i] = item.SumPtNeutral
            self.elec_sumPtCPU[i] = item.SumPtChargedPU

            i += 1

        self.elec_size[0] = i


    #___________________________________________
    def processMuons(self, muons):
        i = 0
        for item in muons:
            self.muon_pt      [i] = item.PT
            self.muon_eta     [i] = item.Eta
            self.muon_phi     [i] = item.Phi
            self.muon_mass    [i] = item.P4().M()
            self.muon_charge  [i] = item.Charge
            self.muon_dz      [i] = item.DZ
            self.muon_reliso  [i] = item.IsolationVar
            self.muon_relisoRC[i] = item.IsolationVarRhoCorr
            self.muon_sumPt   [i] = item.SumPt
            self.muon_sumPtCh [i] = item.SumPtCharged
            self.muon_sumPtNeu[i] = item.SumPtNeutral
            self.muon_sumPtCPU[i] = item.SumPtChargedPU

            i += 1

        self.muon_size[0] = i

    #___________________________________________
    def processCHSJets(self, jets):

        i = 0
        for item in jets:
            jetp4 = item.P4()
            self.jet_pt        [i] = jetp4.Pt()
            self.jet_eta       [i] = jetp4.Eta()
            self.jet_phi       [i] = jetp4.Phi()
            self.jet_mass      [i] = jetp4.M()
            if not self.prune:
                self.jet_deltaEta  [i] = item.DeltaEta
                self.jet_deltaPhi  [i] = item.DeltaPhi
                self.jet_charge    [i] = item.Charge
                self.jet_hoe       [i] = item.EhadOverEem
                self.jet_nCharged  [i] = item.NCharged
                self.jet_nNeutral  [i] = item.NNeutrals
                self.jet_nef       [i] = item.NeutralEnergyFraction
                self.jet_cef       [i] = item.ChargedEnergyFraction
            # https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc
            self.jet_btag      [i] = item.BTag
            self.jet_btagAlgo  [i] = item.BTagAlgo
            self.jet_btagPhys  [i] = item.BTagPhys
            self.jet_flavor    [i] = item.Flavor
            self.jet_flavorAlgo[i] = item.FlavorAlgo
            self.jet_flavorPhys[i] = item.FlavorPhys

            i += 1
        self.jet_size[0] = i

    # ___________________________________________
    def processRho(self, rho):

        i = 0
        for item in rho:
            self.rho[i] = item.Rho

            i += 1

        self.rho_size[0] = i

    #___________________________________________
    def processMissingET(self, met):
        if met:
            self.met_pt    [0] = met[0].MET
            self.met_phi   [0] = met[0].Phi

    def fill(self):
        self.t.Fill()

    def write(self):
        self.t.Write()

#_______________________________________________________
def dr_match(p1, p2, drmin):
    dr = p1.P4().DeltaR(p2.P4())
    return dr < drmin


#_____________________________________________________________________________________________________________
def main():

    ROOT.gSystem.Load("libDelphes")
    try:
      ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
      ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
    except:
      pass

    bool_type = lambda s: s.lower() in [ 'True', '1' ]
    parser = argparse.ArgumentParser()
    parser.add_argument ('-i', '--input', help='input Delphes file',  default='delphes.root')
    parser.add_argument ('-o', '--output', help='output flat tree',  default='tree.root')
    parser.add_argument ('-n', '--nev', help='number of events', type=int, default=-1)
    parser.add_argument ('-d', '--debug', help='debug flag',  action='store_true',  default=False)
    parser.add_argument ('-p', '--pileup', help='input built with PU', type=bool_type, default=True)
    parser.add_argument ('-P', '--prune', help='prune', type=bool_type, default=True)

    args = parser.parse_args()

    inputFile = args.input
    outputFile = args.output
    nevents = args.nev
    debug = args.debug
    pileup = args.pileup
    prune = args.prune

    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)
    

    # Create object of class ExRootTreeReader
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    ## for now only M for electrons, LT for muons and LT for photons are defined !!
    ## should dervie new parameterisations for other working points

    branchEvent           = treeReader.UseBranch('Event')   
    branchWeight          = treeReader.UseBranch('Weight')   
    if pileup and not prune:
        branchVertex          = treeReader.UseBranch('Vertex')
    branchParticle        = treeReader.UseBranch('Particle') 
    branchGenJet          = treeReader.UseBranch('GenJet')   
    branchGenMissingET    = treeReader.UseBranch('GenMissingET')
    if not prune:
        branchPhoton          = treeReader.UseBranch('Photon')
    branchElectron        = treeReader.UseBranch('Electron')
    branchMuon            = treeReader.UseBranch('Muon')
    branchCHSJet          = treeReader.UseBranch('Jet')
    branchMissingET       = treeReader.UseBranch('MissingET')
    if pileup and not prune:
        branchRho             = treeReader.UseBranch('Rho')

    treeProducer = TreeProducer(pileup, prune, debug)

    if nevents > 0:
        numberOfEntries = nevents

    ################ Start event loop #######################
    for entry in range(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.ReadEntry(entry)

        if (entry+1)%1000 == 0:
            print(' ... processed {} events ...'.format(entry+1))

        treeProducer.processEvent(branchEvent, branchWeight)
        if pileup and not prune:
            treeProducer.processVertices(branchVertex)
        treeProducer.processGenParticles(branchParticle)
        treeProducer.processGenJets(branchGenJet)
        treeProducer.processGenMissingET(branchGenMissingET)
        treeProducer.processElectrons(branchElectron)
        treeProducer.processMuons(branchMuon)
        if not prune:
            treeProducer.processPhotons(branchPhoton)
        treeProducer.processMissingET(branchMissingET)
        treeProducer.processCHSJets(branchCHSJet)
        if pileup and not prune:
            treeProducer.processRho(branchRho)

        ## fill tree 
        treeProducer.fill()

    out_root = ROOT.TFile(outputFile,"RECREATE")
    out_root.cd()
    treeProducer.write()
 

#_______________________________________________________________________________________
if __name__ == "__main__":
    main()

