#!/usr/bin/env python


print(__doc__)

import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.model_selection import train_test_split
import pandas, root_numpy
from sklearn.ensemble import GradientBoostingClassifier
import sklearn
from pdb import set_trace
from sklearn.externals import joblib
import argparse
import pickle



features = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
            'var_MaxD0SigBS','var_maxMuonsDca','var_minMuonsDca',
            'var_segCompMuMin','var_Iso08MuMin',
            'var_Iso08MuMax','var_MaxdeltaMuZ', 'var_MaxVertexPairQuality']




sig = pandas.DataFrame(
    root_numpy.root2array(
    '../tmva/TMVATress.root', 'TreeS_Ds',
    # 'signal.root', 'tree',
    branches = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                'var_MaxD0SigBS','var_maxMuonsDca','var_minMuonsDca',
                'var_segCompMuMin','var_Iso08MuMin',
                'var_Iso08MuMax','var_MaxdeltaMuZ', 'var_MaxVertexPairQuality'],

    selection=""
    )
    )
sig['target'] = np.ones(sig.shape[0])

print 'loaded signal'


bkg = pandas.DataFrame(
    root_numpy.root2array(
    ['../tmva/TMVATress.root'], 'TreeB',
    # 'DATI.root', 'tree',
    branches = ['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                'var_MaxD0SigBS','var_maxMuonsDca','var_minMuonsDca',
                'var_segCompMuMin','var_Iso08MuMin',
                'var_Iso08MuMax','var_MaxdeltaMuZ', 'var_MaxVertexPairQuality'],

    selection=""
    )
    )
bkg['target'] = np.zeros(bkg.shape[0])

print 'loaded bkg'


print 'Summary:'
print '# sig:', sig.shape[0]
print '# bkg:', bkg.shape[0]


X=bkg
y=sig

#iris = datasets.load_iris()

#X = iris.data
#y = iris.target
target_names = ['bkg', 'sig']

print '1'

pca = PCA(n_components=2)
print '2'
X_r = pca.fit(X).transform(X)
print '3'
#lda = LinearDiscriminantAnalysis(n_components=2)
print '4'
#X_r2 = lda.fit(X, y).transform(X)
print '5'

# Percentage of variance explained for each components
print('explained variance ratio (first two components): %s'
      % str(pca.explained_variance_ratio_))

plt.figure()
colors = ['navy', 'turquoise']
lw = 2

for color, i, target_name in zip(colors, [0, 1, 2], target_names):
    plt.scatter(X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=.8, lw=lw,
                label=target_name)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('PCA of IRIS dataset')

plt.figure()
#for color, i, target_name in zip(colors, [0, 1, 2], target_names):
#    plt.scatter(X_r2[y == i, 0], X_r2[y == i, 1], alpha=.8, color=color,
#                label=target_name)
#plt.legend(loc='best', shadow=False, scatterpoints=1)
#plt.title('LDA of IRIS dataset')

plt.show()
