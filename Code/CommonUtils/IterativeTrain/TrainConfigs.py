configuration=[]



varsets1 = {'A':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_maxMuonsDca','var_MindcaTrackSV'
             ],
                
                
            'B':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_Muon1DetID','var_Muon2DetID','var_Muon3DetID' ,'var_MaxVertexPairQuality'
             ],

            
            'C':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_maxMuonsDca'
             ]}







varsets2 = {'A':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_maxMuonsDca',
             ],
                
                
            'B':['var_vertexKFChi2', 'var_svpvTauAngle', 
                 'var_Muon1DetID','var_Muon2DetID' ,
             ],

                
            'C':['var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_maxMuonsDca'
             ]}






varsets3 = {'A':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_maxMuonsDca','var_MindcaTrackSV'
             ],
                
                
            'B':['var_vertexKFChi2', 'var_svpvTauAngle', 
                 'var_Muon1DetID','var_Muon2DetID' ,'var_MindcaTrackSV'
             ],

                
            'C':['var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_maxMuonsDca','var_MindcaTrackSV'
             ]}






varsets4 = {'A':['var_vertexKFChi2',  'var_flightLenSig',
                 'var_MaxD0SigSV','var_maxMuonsDca','var_MindcaTrackSV'
             ],
                
                
            'B':['var_vertexKFChi2', 
                 'var_Muon1DetID','var_Muon2DetID' ,'var_MindcaTrackSV'
             ],

                
            'C':['var_svpvTauAngle', 
                 'var_MaxD0SigSV','var_maxMuonsDca','var_MindcaTrackSV'
             ]}




configuration.append(varsets1)
configuration.append(varsets2)
configuration.append(varsets3)
configuration.append(varsets4)



selection = {'var_max_tKink': [0,80],
              'var_MindcaTrackSV': [0,0.5],
              'var_vertexKFChi2': [0,50.0],
              'var_svpvTauAngle': [0,0.2],
              'var_flightLenSig': [0,100],
             }
