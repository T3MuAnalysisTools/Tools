configuration=[]




# -- 0
varsets0 = {'ZTT_tau_CV_3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau','var_HPS_GJ_Angle_Ratio'
                        ]}
                        
# -- 1 Remove var_HPS_GJ_Angle_Ratio
varsets1 = {'ZTT_tau_CV_3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
              
# -- 2 Remove var_mu3_eta
varsets2 = {'ZTT_tau_CV_3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        

# -- 3 Remove var_DiTauMass_Collinear
varsets3 = {'ZTT_tau_CV_3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
# -- 4 Remove var_mu3_pT
varsets4 = {'ZTT_tau_CV_3mu':['var_mu1_pT','var_mu2_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
                        
# -- 5 Remove var_TripletEta
varsets5 = {'ZTT_tau_CV_3mu':['var_mu1_pT','var_mu2_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
                        
# -- 6 Remove var_mu1_pT
varsets6 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_mu1_eta','var_mu2_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
# -- 7 Remove var_mu2_eta
varsets7 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_mu1_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}


# -- 8 Remove var_Phi_To_Opposite_Side
varsets8 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_mu1_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}


# -- 9 Remove var_Tau_eta
varsets9 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_mu1_eta','var_TripletPT','var_Tau3MuIsolation',
                         'var_Tau_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}


# -- 10 Remove var_mu1_eta
varsets10 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Tau_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_MET_Et',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
# -- 11 Remove var_MET_Et
varsets11 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_Tau_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
# -- 12 Remove var_Tau_pT
varsets12 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
# -- 13 Remove var_DeltaPhi
varsets13 = {'ZTT_tau_CV_3mu':['var_mu2_pT','var_TripletPT','var_Tau3MuIsolation',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
                        
# -- 14 Remove var_mu2_pT
varsets14 = {'ZTT_tau_CV_3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau'
                        ]}
                        
                        
# -- 15 Remove var_HPS_Inv_Mass_Z_Tau3mu_SpecTau
varsets15 = {'ZTT_tau_CV_3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass',
                         'var_HPS_FL_Sig'
                        ]}
                        
# -- 16 Remove var_HPS_FL_Sig
varsets16 = {'ZTT_tau_CV_3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass'
                        ]}

# -- 17 Remove var_FLSignificance
varsets17 = {'ZTT_tau_CV_3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack',
                         'var_VisMass'
                        ]}
                        
                        
# -- 18 Remove var_MinDistToIsoTrack
varsets18 = {'ZTT_tau_CV_3mu':['var_TripletPT','var_Tau3MuIsolation',
                         'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF',
                         'var_VisMass'
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


"""
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
             'var_FLSignificance': [0,150],
             'var_MinDistToIsoTrack': [0,10.0],
             'var_HPS_FL_Sig': [0,100],
             #'var_HPS_Inv_Mass_Z_Tau3mu_SpecTau': [0,180],
             'var_HPS_GJ_Angle_Ratio': [0,90],
             'var_MuonIsolation': [0,12],
             'var_4Mu_Chi2': [0,20],
             #'var_4Mu_Vertex_Disp': [0,5.0],
             'var_4Mu_MinDistToIsoTrack_mm': [0,10.0],
             'var_ElectronSumIsolation': [0,4],
            }
"""
            
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
             'var_FLSignificance': [0,150],
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