configuration=[]




# -- 0
varsets0 = {'ZTT_e3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
# -- 1 Remove var_mu1_eta
varsets1 = {'ZTT_e3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
              
# -- 2 Remove var_mu2_eta
varsets2 = {'ZTT_e3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        

# -- 3 Remove var_mu3_pT
varsets3 = {'ZTT_e3mu':['var_mu1_pT','var_mu2_pT','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
                        
# -- 4 Remove var_mu3_eta
varsets4 = {'ZTT_e3mu':['var_mu1_pT','var_mu2_pT','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
                        
                        
# -- 5 Remove var_TripletEta
varsets5 = {'ZTT_e3mu':['var_mu1_pT','var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
                        
                        
# -- 6 Remove var_mu1_pT
varsets6 = {'ZTT_e3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
                        
# -- 7 Remove var_MET_Et
varsets7 = {'ZTT_e3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}


# -- 8 Remove var_Electron_eta
varsets8 = {'ZTT_e3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}


# -- 9 Remove var_mu2_pT
varsets9 = {'ZTT_e3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}


# -- 10 Remove var_FLSignificance
varsets10 = {'ZTT_e3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT',
                         'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
# -- 11 Remove var_Phi_To_Opposite_Side
varsets11 = {'ZTT_e3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT',
                         'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
                        
# -- 12 Remove var_SVPVTauDirAngle
varsets12 = {'ZTT_e3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_Electron_pT',
                         'var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
                        
# -- 13 Remove var_Electron_pT
varsets13 = {'ZTT_e3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
                        
                        
# -- 14 Remove var_TripletPT
varsets14 = {'ZTT_e3mu':['var_Tau3MuIsolation',
                         'var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
# -- 15 Remove var_MinDistToIsoTrack
varsets15 = {'ZTT_e3mu':['var_Tau3MuIsolation',
                         'var_ThreeMuVertexChi2KF','var_DeltaPhi',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}
                        
# -- 16 Remove var_VisMass
varsets16 = {'ZTT_e3mu':['var_Tau3MuIsolation',
                         'var_ThreeMuVertexChi2KF','var_DeltaPhi',
                         'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}




configuration.append(varsets0)
configuration.append(varsets1)
configuration.append(varsets2)
configuration.append(varsets3)
configuration.append(varsets4)
configuration.append(varsets5)
configuration.append(varsets6)
configuration.append(varsets7)
configuration.append(varsets8)
configuration.append(varsets9)
configuration.append(varsets10)
configuration.append(varsets11)
configuration.append(varsets12)
configuration.append(varsets13)
configuration.append(varsets14)
configuration.append(varsets15)
configuration.append(varsets16)



selection = {
             #'var_mu1_pT': [0,60],
             #'var_mu2_pT': [0,30],
             #'var_mu3_pT': [0,20],
             #'var_TripletPT': [0,90],
             #'var_VisMass': [70,180],
             #'var_DiTauMass_Collinear': [70,180],
             #'var_ThreeMuVertexChi2KF': [0,45],
             #'var_Tau_pT': [0,100],
             'var_Muon_pT': [0,65],
             'var_Electron_pT': [0,70],
             #'var_ThreeMuVertexChi2KF': [0,80],
             'var_SVPVTauDirAngle': [0,0.4],
             'var_Tau3MuIsolation': [0,10],
             #'var_Tau_pT': [0,80],
             'var_FLSignificance': [0,20],
             'var_MinDistToIsoTrack': [0,1.45],
             'var_HPS_FL_Sig': [0,50],
             #'var_HPS_Inv_Mass_Z_Tau3mu_SpecTau': [0,180],
             'var_HPS_GJ_Angle_Ratio': [0,20],
             'var_MuonIsolation': [0,12],
             'var_4Mu_Chi2': [0,20],
             #'var_4Mu_Vertex_Disp': [0,5.0],
             'var_4Mu_MinDistToIsoTrack_mm': [0,10.0],
             'var_ElectronSumIsolation': [0,4],
            }
