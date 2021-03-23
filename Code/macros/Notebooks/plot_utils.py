PLOT_DIR = './plots/'

# compare ROC curves
def compare_rocs(decisions_arr, y_test_arr, event_weights_test_arr, titles, plot_name):
    plt.clf()
    
    for i, decision in enumerate(decisions_arr):
        fpr, tpr, thresholds = roc_curve(y_test_arr[i], decision, sample_weight=event_weights_test_arr[i])
        roc_auc = auc(fpr, tpr)
        plt.plot(tpr, 1-fpr, lw=1, label=titles[i]+' ROC (area = %0.2f)'%(roc_auc))


    plt.plot([1, 0], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Background rejection')
    plt.ylabel('Signal acceptance')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower left")
    plt.grid()
    plt.figure(figsize=(6.2,6))
    plt.tight_layout()
    cms = plt.text(
        0.05, 0.9, u"CMS $\it{Simulations}$",
        fontsize=18, fontweight='bold',
        horizontalalignment='left', 
        verticalalignment='bottom', 
        transform=ax.transAxes
    )
    plt.savefig(PLOT_DIR+plot_name, format='svg', bbox_inches='tight')

# function to compare train and test predictions
def compare_train_test(model, model_type, X_train, y_train, X_test, y_test, event_weights_train, event_weights_test, plot_name, bins=30):
    plt.clf()
    decisions = [[], [], [], []]
    weights = [[],[],[],[]]
    if (model_type=='dnn'):
        predictions_train = model.predict(X_train)
        predictions_test = model.predict(X_test)
    elif (model_type=='bdt'):
        predictions_train = model.decision_function(X_train)
        predictions_test = model.decision_function(X_test)

    for j in range(len(X_train)):
        if(y_train[j]>0.5):
            decisions[0].append(predictions_train[j].item())
            weights[0].append(event_weights_train[j].item())
        else:
            decisions[1].append(predictions_train[j].item())
            weights[1].append(event_weights_train[j].item())

    for j in range(len(X_test)):
        if(y_test[j]>0.5):
            weights[2].append(event_weights_test[j].item())
            decisions[2].append(predictions_test[j].item())
        else:
            weights[3].append(event_weights_test[j].item())
            decisions[3].append(predictions_test[j].item())

    low = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    low_high = (low,high)

    plt.hist(decisions[0],
             color='b', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density=True,
             label='S (train)', weights=weights[0])
    plt.hist(decisions[1],
             color='r', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', density = True,
             label='B (train)', weights=weights[1])

    hist, bins = np.histogram(decisions[2],
                              bins=bins, range=low_high, density=True, weights=weights[2])
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale

    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='S (test)')

    hist, bins = np.histogram(decisions[3],
                              bins=bins, range=low_high, density=True, weights=weights[3])
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale

    plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='B (test)')

    if (model_type=='dnn'): plt.xlabel("NN output")
    elif (model_type=='bdt'): plt.xlabel("BDT output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc='upper center')
    #plt.show()
    plt.savefig(PLOT_DIR+model_type+'_'+plot_name, format='svg')
    
