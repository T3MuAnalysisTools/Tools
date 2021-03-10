configuration=[]





varsets3 = {'A':['var_flightLenSig', 'var_svpvTauAngle','var_MindcaTrackSV',
                 'var_NtracksClose','var_nsv','var_MaxD0SigBS','var_MinD0SigBS','var_Iso08MuMin'
                 ,'var_dcaTrackPV','var_MinMuonImpactAngle','var_flightLenDist'
             ]}
                




configuration.append(varsets3)




selection = {'var_max_tKink': [0,80],
              'var_MindcaTrackSV': [0,0.5],
              'var_vertexKFChi2': [0,50.0],
              'var_svpvTauAngle': [0,0.2],
              'var_flightLenSig': [0,100],
              'var_IsoKStarMass_Mu1': [0,4],
              'var_IsoKStarMass_Mu2': [0,4],
              'var_IsoKStarMass_Mu3': [0,4],
              'var_IsoPhiKKMass_Mu1': [0,4],
              'var_IsoPhiKKMass_Mu2': [0,4],
              'var_IsoPhiKKMass_Mu3': [0,4],
              'var_IsoMuMuMass_Mu1': [0,4],
              'var_IsoMuMuMass_Mu2': [0,4],
              'var_IsoMuMuMass_Mu3': [0,4],
             }
