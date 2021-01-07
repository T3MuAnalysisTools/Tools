var_limits = {  'var_min_p':[0,10],
        'var_max_tKink': [0,80],
        'var_mindca_iso': [0,0.5],
        'var_trk_relPt': [0,10],
        'var_vertexKFChi2': [0,100.0],
        'var_svpvTauAngle': [0,0.2],
        'var_flightLenSig': [0,100],
        'var_segCompMuMin': [0.2,1],
        'var_DNNSegCompMuMin': [0.2,1],
        'var_MinD0Significance': [0,15],
        'var_pmin': [0,50],
        'var_max_cLP': [0,30],
        'var_minMatchedStations': [0,4],
        'var_trackerMuonId': [-0.4,0.4],
        'var_nMatches_mu3':[0,10],
        'muon1_seg_comp_dnn':[0,1.0],
        'muon2_seg_comp_dnn':[0,1.0],
        'muon3_seg_comp_dnn':[0,1.0],
        'Muon1_segmentCompatibility':[0,1.0],
        'Muon2_segmentCompatibility': [0,1.0],
        'Muon3_segmentCompatibility':[0,1.0],
        'abs(var_Muon3_timeAtIpInOut)': [0,100],
        'trackerMu_min_segmComp': [-0.3,0.3],
        'trackerMuonId':[-0.5,0.5],
        'globalMuon1Id': [0.1,0.8],
        'globalMuon2Id': [0.1,0.8]
        }
'''
varsets = {
        '2016vars':{
            'training':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt', 'var_minMatchedStations', 'var_segCompMuMin'],
            'spectator': ['var_tauMass', 'threeGlobal']
            },
        'TrackerMuonId':{
            'training': ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt', 'var_minMatchedStations','trackerMuonId', 'var_segCompMuMin'],
            'spectator': ['var_tauMass', 'threeGlobal']
            },
        'BDTSegmentComp':{
            'training': ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt', 'var_minMatchedStations','trackerMuonId', 'trackerMu_min_segmComp'],
            'spectator': ['var_tauMass', 'threeGlobal']
            },
        'Muon3TimeAtIp':{
            'training': ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt', 'var_minMatchedStations', 'abs(var_Muon3_timeAtIpInOut)','trackerMuonId', 'var_segCompMuMin'],
            'spectator': ['var_tauMass', 'threeGlobal']
            },
        'Muon3TimeAtIp_MinBDTSegmComp':{
            'training': ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt', 'var_minMatchedStations', 'abs(var_Muon3_timeAtIpInOut)','trackerMuonId', 'trackerMu_min_segmComp'],
            'spectator': ['var_tauMass', 'threeGlobal']
            }
        }
'''

varsets = {
         'TrackerMuonId_globalMuonId':{
            'training': ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt', 'var_minMatchedStations','trackerMuonId', 'var_segCompMuMin', 'globalMuon1Id', 'globalMuon2Id'],
            'spectator': ['var_tauMass', 'threeGlobal']
            },
        'BDTSegmentComp_globalMuonId':{
            'training': ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt', 'var_minMatchedStations','trackerMuonId', 'trackerMu_min_segmComp', 'globalMuon1Id', 'globalMuon2Id'],
            'spectator': ['var_tauMass', 'threeGlobal']
            }
        }

var_list = {
        'var_min_p':{ 'limits': [0,10], 'type': 'F' }, 
        'var_max_tKink':{ 'limits': [0,80], 'type': 'F'},
        'var_mindca_iso':{ 'limits': [0,0.5], 'type': 'F'},
        'var_trk_relPt':{ 'limits': [0,10], 'type': 'F'},
        'var_vertexKFChi2':{ 'limits': [0,100.0], 'type': 'F'},
        'var_svpvTauAngle':{ 'limits': [0,0.2], 'type': 'F'},
        'var_flightLenSig':{ 'limits': [0,100], 'type': 'F'},
        'var_segCompMuMin':{ 'limits': [0.2,1], 'type': 'F'},
        'var_DNNSegCompMuMin':{ 'limits': [0.2,1], 'type': 'F'},
        'var_MinD0Significance':{ 'limits': [0,15], 'type': 'F'},
        'var_pmin':{ 'limits': [0,50], 'type': 'F'},
        'var_max_cLP':{ 'limits': [0,30], 'type': 'F'},
        'var_minMatchedStations':{ 'limits': [0,4], 'type': 'I'},
        'var_trackerMuonId':{ 'limits': [-0.4,0.4], 'type': 'F'},
        'var_nMatches_mu3':{ 'limits':[0,10], 'type': 'I'},
        'muon1_seg_comp_dnn':{ 'limits':[0,1.0], 'type': 'F'},
        'muon2_seg_comp_dnn':{ 'limits':[0,1.0], 'type': 'F'},
        'muon3_seg_comp_dnn':{ 'limits':[0,1.0], 'type': 'F'},
        'Muon1_segmentCompatibility':{ 'limits':[0,1.0], 'type': 'F'},
        'Muon2_segmentCompatibility':{ 'limits': [0,1.0], 'type': 'F'},
        'Muon3_segmentCompatibility':{ 'limits':[0,1.0], 'type': 'F'},
        'abs(var_Muon3_timeAtIpInOut)':{ 'limits': [0,100], 'type': 'F'},
        'trackerMu_min_segmComp':{ 'limits': [-0.3,0.3], 'type': 'F'},
        'trackerMuonId':{ 'limits':[-0.5,0.5], 'type': 'F'},
        'globalMuon1Id':{ 'limits':[0.1,0.8], 'type': 'F'},
        'globalMuon2Id':{ 'limits':[0.1,0.8], 'type': 'F'}
        }
