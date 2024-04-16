configuration=[]


#This configuration includes all the vars I may be comparing

# -- 0
varsets0 = {'ZTT_tau_NoCV_3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ],
                        
            'ZTT_tau_CV_3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_HPS_FL_Sig','var_HPS_Inv_Mass_Z_Tau3mu_SpecTau','var_HPS_GJ_Angle_Ratio'
                        ],
                        
            'ZTT_mu3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Muon_pT', 'var_Muon_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_MuonIsolation','var_4Mu_Chi2', 'var_4Mu_Vertex_Disp', 'var_4Mu_MinDistToIsoTrack_mm'
                        ],
                
            'ZTT_e3mu':['var_mu1_pT','var_mu2_pT','var_mu3_pT','var_mu1_eta','var_mu2_eta','var_mu3_eta','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Electron_pT','var_Electron_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi','var_MinDistToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et', 
                         'var_VisMass', 'var_DiTauMass_Collinear',
                         'var_ElectronSumIsolation'
                        ]}




configuration.append(varsets0)

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

'''
selection = {'var_Tau3MuIsolation': [0,80],
             'var_mu1_pT': [0,80],
             'var_mu2_pT': [0,60],
             'var_mu3_pT': [0,60],
             'var_TripletPT': [0,30],
             'var_TripletEta': [0,0.4],
             'var_mu1_eta': [0,80],
             'var_mu2_eta': [0,80],
             'var_mu3_eta': [0,60],
             'var_FLSignificance': [0,60],
             'var_SVPVTauDirAngle': [0,0.4],
             'var_ThreeMuVertexChi2KF': [0,80],
             'var_DeltaPhi': [0,80],
             'var_Phi_To_Opposite_Side': [0,80],
             'var_VisMass': [0,60],
             'var_DiTauMass_Collinear': [0,60],
             'var_ThreeMuVertexChi2KF': [0,80],
             'var_Tau_pT': [0,80],
             'var_Muon_pT': [0,60],
             'var_Electron_pT': [0,60],
             'var_ElectronSumIsolation': [0,30],
             'var_MuonIsolation': [0,5],
             'var_Muon_eta': [0,60],
             'var_Tau_eta': [0,30],
             'var_Electron_eta': [0,0.4],
             'var_4Mu_Chi2': [0,70],
             'var_4Mu_MinDistToIsoTrack_mm': [0,1.0],
             'var_4Mu_Vertex_Disp': [0,6],
             'var_HPS_FL_Sig': [0,23],
             'var_HPS_Inv_Mass_Z_Tau3mu_SpecTau': [0,180],
             'var_HPS_GJ_Angle_Ratio': [0,30]
            }
'''