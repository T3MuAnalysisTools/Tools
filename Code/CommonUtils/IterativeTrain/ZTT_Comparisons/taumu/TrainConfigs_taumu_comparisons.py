configuration=[]




# -- 0
varsets0 = {'ZTT_mu3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_Vertex_Disp', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
# -- 1 Remove var_mu3_eta
varsets1 = {'ZTT_mu3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_Vertex_Disp', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
              
# -- 2 Remove var_SVPVTauDirAngle
varsets2 = {'ZTT_mu3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_Vertex_Disp', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        

# -- 3 Remove var_DiTauMass_Collinear
varsets3 = {'ZTT_mu3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_Vertex_Disp', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
# -- 4 Remove var_mu1_pT
varsets4 = {'ZTT_mu3mu':['var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_Vertex_Disp', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
                        
# -- 5 Remove var_TripletEta
varsets5 = {'ZTT_mu3mu':['var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_Vertex_Disp', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
                        
# -- 6 Remove var_4Mu_Vertex_Disp
varsets6 = {'ZTT_mu3mu':['var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
# -- 7 Remove var_4Mu_Chi2
varsets7 = {'ZTT_mu3mu':['var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}


# -- 8 Remove var_mu3_pT
varsets8 = {'ZTT_mu3mu':['var_mu2_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}


# -- 9 Remove var_mu1_eta
varsets9 = {'ZTT_mu3mu':['var_mu2_pT','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}


# -- 10 Remove var_Phi_To_Opposite_Side
varsets10 = {'ZTT_mu3mu':['var_mu2_pT','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
# -- 11 Remove var_MET_Et
varsets11 = {'ZTT_mu3mu':['var_mu2_pT','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
# -- 12 Remove var_DeltaPhi
varsets12 = {'ZTT_mu3mu':['var_mu2_pT','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
# -- 13 Remove var_FLSignificance
varsets13 = {'ZTT_mu3mu':['var_mu2_pT','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
                        
# -- 14 Remove var_mu2_eta
varsets14 = {'ZTT_mu3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
# -- 15 Remove var_ThreeMuVertexChi2KF
varsets15 = {'ZTT_mu3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
# -- 16 Remove var_mu2_pT
varsets16 = {'ZTT_mu3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
# -- 17 Remove var_TripletPT
varsets17 = {'ZTT_mu3mu':['var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
# -- 18 Remove var_MinDistToIsoTrack
varsets18 = {'ZTT_mu3mu':['var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_VisMass',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
                        ]}
                        
                        
                        
# -- 19 Remove var_VisMass
varsets19 = {'ZTT_mu3mu':['var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_MuonIsolation', 'var_4Mu_MinDistToIsoTrack_mm'
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
configuration.append(varsets17)
configuration.append(varsets18)
configuration.append(varsets19)



selection = {'var_mu1_pT': [0,60],
             'var_mu2_pT': [0,30],
             'var_mu3_pT': [0,20],
             'var_TripletPT': [0,90],
             'var_VisMass': [70,180],
             'var_DiTauMass_Collinear': [70,180],
             'var_ThreeMuVertexChi2KF': [0,45],
             'var_Tau_pT': [0,100],
             'var_Muon_pT': [0,100],
             'var_Electron_pT': [0,80],
             'var_ThreeMuVertexChi2KF': [0,80],
             'var_SVPVTauDirAngle': [0,0.4],
             'var_Tau3MuIsolation': [0,4],
             'var_Tau_pT': [0,80],
             'var_FLSignificance': [0,60],
             'var_MinDistToIsoTrack': [0,1.5],
             'var_HPS_FL_Sig': [0,30],
             'var_HPS_Inv_Mass_Z_Tau3mu_SpecTau': [0,180],
             'var_HPS_GJ_Angle_Ratio': [0,100],
             'var_MuonIsolation': [0,5],
             'var_4Mu_Chi2': [0,70],
             'var_4Mu_Vertex_Disp': [0,1.2],
             'var_4Mu_MinDistToIsoTrack_mm': [0,1.0],
             'var_ElectronSumIsolation': [0,1.6],
            }
