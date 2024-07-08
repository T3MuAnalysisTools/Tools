configuration=[]




# -- 0
varsets0 = {'ZTT_tau_NoCV_3mu':['var_mu3_pT','var_TripletPT','var_TripletEta','var_Tau3MuIsolation',
                         'var_Tau_pT','var_Tau_eta',
                         'var_FLSignificance', 'var_SVPVTauDirAngle','var_ThreeMuVertexChi2KF', 'var_DeltaPhi','var_MinDrToIsoTrack','var_Phi_To_Opposite_Side',
                         'var_MET_Et',
                         'var_VisMass', 'var_DiTauMass_Collinear'
                        ]}
                        
#Removed var_MET_Et
varsets1 = {'ZTT_tau_NoCV_3mu':['var_mu3_pT', 'var_TripletPT', 'var_TripletEta', 'var_Tau3MuIsolation', 'var_Tau_pT', 
            'var_Tau_eta', 'var_FLSignificance', 'var_SVPVTauDirAngle', 'var_ThreeMuVertexChi2KF', 'var_DeltaPhi', 
            'var_MinDrToIsoTrack', 'var_Phi_To_Opposite_Side', 'var_VisMass', 'var_DiTauMass_Collinear']}

#Removed var_Tau_eta
varsets2 = {'ZTT_tau_NoCV_3mu':['var_mu3_pT', 'var_TripletPT', 'var_TripletEta', 'var_Tau3MuIsolation', 'var_Tau_pT', 
            'var_FLSignificance', 'var_SVPVTauDirAngle', 'var_ThreeMuVertexChi2KF', 'var_DeltaPhi', 'var_MinDrToIsoTrack', 
            'var_Phi_To_Opposite_Side', 'var_VisMass', 'var_DiTauMass_Collinear']}

#Removed var_mu3_pT
varsets3 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_TripletEta', 'var_Tau3MuIsolation', 'var_Tau_pT', 'var_FLSignificance', 
            'var_SVPVTauDirAngle', 'var_ThreeMuVertexChi2KF', 'var_DeltaPhi', 'var_MinDrToIsoTrack', 'var_Phi_To_Opposite_Side', 
            'var_VisMass', 'var_DiTauMass_Collinear']}

#Removed var_DiTauMass_Collinear
varsets4 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_TripletEta', 'var_Tau3MuIsolation', 'var_Tau_pT', 'var_FLSignificance', 
            'var_SVPVTauDirAngle', 'var_ThreeMuVertexChi2KF', 'var_DeltaPhi', 'var_MinDrToIsoTrack', 'var_Phi_To_Opposite_Side', 
            'var_VisMass']}

#Removed var_TripletEta
varsets5 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_Tau3MuIsolation', 'var_Tau_pT', 'var_FLSignificance', 'var_SVPVTauDirAngle', 
            'var_ThreeMuVertexChi2KF', 'var_DeltaPhi', 'var_MinDrToIsoTrack', 'var_Phi_To_Opposite_Side', 'var_VisMass']}

#Removed var_Tau_pT
varsets6 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_Tau3MuIsolation', 'var_FLSignificance', 'var_SVPVTauDirAngle', 'var_ThreeMuVertexChi2KF', 
            'var_DeltaPhi', 'var_MinDrToIsoTrack', 'var_Phi_To_Opposite_Side', 'var_VisMass']}

#Removed var_Phi_To_Opposite_Side
varsets7 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_Tau3MuIsolation', 'var_FLSignificance', 'var_SVPVTauDirAngle', 'var_ThreeMuVertexChi2KF', 
            'var_DeltaPhi', 'var_MinDrToIsoTrack', 'var_VisMass']}

#Removed var_DeltaPhi
varsets8 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_Tau3MuIsolation', 'var_FLSignificance', 'var_SVPVTauDirAngle', 'var_ThreeMuVertexChi2KF', 
            'var_MinDrToIsoTrack', 'var_VisMass']}

#Removed var_SVPVTauDirAngle
varsets9 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_Tau3MuIsolation', 'var_FLSignificance', 'var_ThreeMuVertexChi2KF', 'var_MinDrToIsoTrack', 
            'var_VisMass']}

#Removed var_FLSignificance
varsets10 = {'ZTT_tau_NoCV_3mu':['var_TripletPT', 'var_Tau3MuIsolation', 'var_ThreeMuVertexChi2KF', 'var_MinDrToIsoTrack', 'var_VisMass']}


                        

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



# doesn't really matter. Selection cuts are applied as overflow bins in the analysis instance
selection = {
             #'var_mu1_pT': [0,60],
             #'var_mu2_pT': [0,30],
             #'var_mu3_pT': [0,20],
             #'var_TripletPT': [0,90],
             #'var_VisMass': [70,180],
             #'var_DiTauMass_Collinear': [70,180],
             #'var_ThreeMuVertexChi2KF': [0,45],
             #'var_Tau_pT': [0,100],
             'var_Muon_pT': [-0.01,65],
             'var_Electron_pT': [-0.01,70],
             #'var_ThreeMuVertexChi2KF': [0,80],
             'var_SVPVTauDirAngle': [-0.01,0.4],
             'var_Tau3MuIsolation': [-0.01,8],
             #'var_Tau_pT': [0,80],
             'var_FLSignificance': [-0.01,100],
             'var_MinDrToIsoTrack': [-0.01,1.45],
             'var_HPS_FL_Sig': [-0.01,50],
             #'var_HPS_Inv_Mass_Z_Tau3mu_SpecTau': [0,180],
             'var_HPS_GJ_Angle_Ratio': [-0.01,25],
             'var_MuonIsolation': [-0.01,12],
             'var_4Mu_Chi2': [-0.1,3000],
             'var_4Mu_Vertex_Disp': [-0.01,2.5],
             'var_AvgDeltaZ_3Mu_Mu_mm': [-0.1,3.0],
             'var_ElectronSumIsolation': [-0.01,4],
            }

