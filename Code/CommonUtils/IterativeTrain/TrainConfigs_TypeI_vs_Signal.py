configuration=[]



# -- 0
varsets0 = {'A':['var_svpvTauAngle',
                 'var_pTMu1OverMass_TRF','var_OSSS1Angle_TRF','var_OSSS2Angle_TRF','var_cTheta_TRF_SSSS','var_cTheta_MuonOS_TauPol_TRF'
             ],
                

            'B':['var_pTMu1OverMass_TRF','var_OSSS1Angle_TRF','var_OSSS2Angle_TRF','var_cTheta_TRF_SSSS'
             ],
                
            'C':['var_pTMu1OverMass_TRF','var_OSSS1Angle_TRF','var_OSSS2Angle_TRF'
                 
             ]}



configuration.append(varsets0)
#configuration.append(varsets1)
#configuration.append(varsets2)
#configuration.append(varsets3)
#configuration.append(varsets4)
#configuration.append(varsets5)
#configuration.append(varsets6)
#configuration.append(varsets7)
#configuration.append(varsets8)
#configuration.append(varsets9)
#configuration.append(varsets10)
#configuration.append(varsets11)
#configuration.append(varsets12)
#configuration.append(varsets13)
#configuration.append(varsets14)
#configuration.append(varsets15)
#configuration.append(varsets16)
#configuration.append(varsets17)
#configuration.append(varsets18)




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
             'var_Vertex2muTrkKF' : [-1.5,60],
             }
