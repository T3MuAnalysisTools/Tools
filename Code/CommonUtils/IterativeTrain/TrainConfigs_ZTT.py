configuration=[]




# -- 1
varsets0 = {'ZTT_mu3mu':['var_Tau3MuIsolation', 'var_mu1_pT','var_mu2_pT','var_mu3_pT',
                         'var_MuonIsolation','var_Muon_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi',
                         'var_MET_Et', 'var_MET_Phi',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ],
                

            'ZTT_tau3mu':['var_Tau3MuIsolation', 'var_mu1_pT','var_mu2_pT','var_mu3_pT',
                         'var_Tau_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF',
                          'var_DeltaPhi',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ],
                
            'ZTT_e3mu':['var_Tau3MuIsolation', 'var_mu1_pT','var_mu2_pT','var_mu3_pT',
                         'var_ElectronSumIsolation','var_Electron_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_DeltaPhi',
                         'var_MET_Et', 'var_MET_Phi',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ]}
                        
# -- 1
varsets1 = {'ZTT_mu3mu':['var_TripletPT', 'var_Tau3MuIsolation',
                         'var_MuonIsolation','var_Muon_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack','var_DeltaPhi',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ],
                

            'ZTT_tau3mu':['var_TripletPT', 'var_Tau3MuIsolation',
                         'var_Tau_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack','var_DeltaPhi',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ],
                
            'ZTT_e3mu':['var_TripletPT', 'var_Tau3MuIsolation',
                         'var_ElectronSumIsolation','var_Electron_pT',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF','var_MinDistToIsoTrack','var_DeltaPhi',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ]}




configuration.append(varsets0)
configuration.append(varsets1)



selection = {'var_ThreeMuVertexChi2KF': [0,80]
            }
