#!/usr/bin/env python



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

parser = argparse.ArgumentParser()
parser.add_argument('--load', help='load pkl instead of training')

args = parser.parse_args()

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

tpr = ((sig.var_flightLenSig < 100) & (sig.var_maxMuonsDca < 0.3)  & (sig.var_MaxVertexPairQuality < 30.) & (sig.var_svpvTauAngle < 0.2) ).sum()
tpr /= float(sig.shape[0])



fpr = ((bkg.var_flightLenSig < 100) & (bkg.var_maxMuonsDca < 0.3)  & (bkg.var_MaxVertexPairQuality < 30.) & (bkg.var_svpvTauAngle < 0.2) ).sum()
fpr /= float(bkg.shape[0])


print 'Default selection -- Signal eff.: %.2f%%, Bkg. eff.: %.2f%%' % (tpr*100, fpr*100)

cutbased_file = open('cut_based.pck', 'w+')
pickle.dump((tpr, fpr), cutbased_file)
cutbased_file.close()



data = pandas.concat([sig, bkg])
train, test = train_test_split(data, test_size=0.33, random_state=42)


clf = GradientBoostingClassifier(
learning_rate=0.01, n_estimators=1000, subsample=0.8, random_state=13,
max_features=len(features), verbose=1,
min_samples_leaf=int(0.01*len(train)),
max_depth=5
)





clf.fit(train[features], train.target)
joblib.dump(clf, 'classifier.pkl', compress=True)

pred = clf.predict_proba(test[features])[:, 1]

bdt = pred.copy()



import itertools
xy = [i*j for i,j in itertools.product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
plt.plot(xy, xy, color='grey', linestyle='--')
plt.xlim([10**-5, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')

#draw baseline point
plt.plot([fpr], [tpr], label='cut based', markerfacecolor='red', marker='o', markersize=10)

plt.xscale('log')

fpr, tpr, _ = roc_curve(test.target, pred)
plt.plot(fpr, tpr, label='BDT', color='b')

roc_file = open('roc.pck', 'w+')
pickle.dump((tpr, fpr), roc_file)
roc_file.close()

plt.legend(loc='best')
plt.grid()
plt.title('ROC')
plt.savefig('roc.png')
plt.savefig('roc.pdf')
plt.clf()

decisions = []
for X,y in ((train[features], train.target), (test[features], test.target)):
    d1 = clf.decision_function(X[y>0.5]).ravel()
    d2 = clf.decision_function(X[y<0.5]).ravel()
    decisions += [d1, d2]



low = min(np.min(d) for d in decisions)
high = max(np.max(d) for d in decisions)
low_high = (low,high)
bins = 50

plt.hist(decisions[0],
         color='r', alpha=0.5, range=low_high, bins=bins,
         histtype='stepfilled', normed=True,
         label='S (train)')
plt.hist(decisions[1],
         color='b', alpha=0.5, range=low_high, bins=bins,
         histtype='stepfilled', normed=True,
         label='B (train)')



hist, bins = np.histogram(decisions[2],
                          bins=bins, range=low_high, normed=True)
scale = len(decisions[2]) / sum(hist)
err = np.sqrt(hist * scale) / scale

width = (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='S (test)')

hist, bins = np.histogram(decisions[3],
                          bins=bins, range=low_high, normed=True)
scale = len(decisions[2]) / sum(hist)
err = np.sqrt(hist * scale) / scale

plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='B (test)')

plt.xlabel("BDT output")
plt.ylabel("Arbitrary units")
plt.legend(loc='best')
plt.savefig('overtrain.pdf')
plt.clf()
