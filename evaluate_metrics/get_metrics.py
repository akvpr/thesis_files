import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, confusion_matrix
from sklearn.metrics import roc_curve, auc, precision_recall_curve





def compare_classifiers_same(fig_save_path, include_all, save_fig):

    print('Comparing other vs own. Include all: ' + str(include_all))

    #first go through the real labels
    variant_file = "../data/testVariants_p2"

    var_label_dict = {}

    with open(variant_file, 'r') as file:

        for line in file:

            line = line.strip()
            fields = line.split('\t')
            #print(fields)
            variant_line = fields[0] + ' ' + fields[1]
            var_label_dict[variant_line] = fields[2]


    #load the MVP predictions
    prediction_file = "mvp/tested_variants.txt"
    mvp_var_score_dict = {}
    total_variants = 0
    total_matches = 0
    with open(prediction_file, 'r') as pred_file:

        for line in pred_file:

            line = line.strip()
            fields = line.split(':')
            if len(fields) < 10:
                print('MVP: Missing prediction; skip')
                continue

            variant_line = fields[0] + ' ' + fields[1]
            #print(fields)
            label = fields[-1]
            score = float(fields[-2])
            #print(fields)
            if variant_line not in var_label_dict:
                print(variant_line)
                print('MVP: not in dict')
            else:

                total_variants += 1
                if var_label_dict[variant_line] == 'neutral':
                    dict_label = 0
                else:
                    dict_label = 1

                mvp_var_score_dict[variant_line] = [score,  dict_label]

                if var_label_dict[variant_line] == label:
                    total_matches += 1

    print('MVP - Total matches is ' + str(total_matches) + ' out of ' + str(total_variants))


    #load the P2's predictions
    prediction_file = "PON-P2_results/p2_test_predictions.txt"
    p2_var_score_dict = {}
    total_variants = 0
    total_matches = 0
    with open(prediction_file, 'r') as p2_file:

        p2_file.readline()

        for line in p2_file:

            line = line.strip()
            fields = line.split('\t')
            variant_line = fields[0]
            #print(fields)
            label = fields[5].lower()
            score = float(fields[3])
            #print(fields)
            if variant_line not in var_label_dict:
                print(variant_line)
                print('PONP2: not in dict')
                continue
            else:
                #if the label is unknown, we can either skip this variant or include the score for it anyway
                if label == 'unknown':
                    if include_all:
                        pass
                    else:
                        continue

                total_variants += 1
                if var_label_dict[variant_line] == 'neutral':
                    dict_label = 0
                else:
                    dict_label = 1

                p2_var_score_dict[variant_line] = [score,  dict_label]

                if var_label_dict[variant_line] == label:
                    total_matches += 1

    print('PONP2 - Total matches is ' + str(total_matches) + ' out of ' + str(total_variants))

    #also load our own predictions
    own_pred_file = "own_predictions.txt"
    own_var_score_dict = {}
    total_variants = 0
    total_matches = 0
    with open(own_pred_file, 'r') as own_file:

        for line in own_file:

            line = line.strip()
            fields = line.split('\t')
            seq_name = fields[0]
            mut_str = fields[1]

            variant_line = seq_name + ' ' + mut_str
            label = fields[3].lower()
            score = float(fields[2])

            if variant_line not in var_label_dict:
                print(variant_line)
                print('Own: not in dict')

            else:

                total_variants += 1
                if var_label_dict[variant_line] == 'neutral':
                    dict_label = 0
                else:
                    dict_label = 1

                own_var_score_dict[variant_line] = [score,  dict_label]

                if var_label_dict[variant_line] == label:
                    total_matches += 1

    print('Model - Total matches is ' + str(total_matches) + ' out of ' + str(total_variants))

    #FATHMM
    true_labels = []
    fathmm_predictions = []
    fathmm_scores = []
    fathmm_var_score_dict = {}
    with open("fathmm/fathmm_predictions_neutral.tab", 'r') as file:
        file.readline()
        for line in file:
            #we want to get the score. field 4. prediction is field 3. some might be missing; this is indicated in field 3
            fields = line.split('\t')
            #remove empty tabs
            fields = [field for field in fields if field]

            if fields[3] == 'TOLERATED':
                fathmm_predictions.append(0)
            elif fields[3] == 'DAMAGING':
                fathmm_predictions.append(1)
            else:
                continue

            if var_label_dict[variant_line] == 'neutral':
                dict_label = 0
            else:
                dict_label = 1

            variant_line = fields[1] + ' ' + fields[2]
            fathmm_var_score_dict[variant_line] = [score,  dict_label]
            score = float(fields[4])
            fathmm_scores.append(score)

            true_labels.append(0)

    with open("fathmm/fathmm_predictions_pathogenic.tab", 'r') as file:
        file.readline()
        for line in file:
            #we want to get the score. field 4. prediction is field 3. some might be missing; this is indicated in field 3
            fields = line.split('\t')
            #remove empty tabs
            fields = [field for field in fields if field]

            if fields[3] == 'TOLERATED':
                fathmm_predictions.append(0)
            elif fields[3] == 'DAMAGING':
                fathmm_predictions.append(1)
            else:
                continue

            if var_label_dict[variant_line] == 'neutral':
                dict_label = 0
            else:
                dict_label = 1

            variant_line = fields[1] + ' ' + fields[2]
            fathmm_var_score_dict[variant_line] = [score,  dict_label]

            score = float(fields[4])
            fathmm_scores.append(score)
            
            true_labels.append(1)

    true_labels = np.array(true_labels)
    fathmm_predictions = np.array(fathmm_predictions)
    fathmm_scores = np.array(fathmm_scores)


    with open("ALL_test_variants.txt", "w") as file:
        for variant in var_label_dict:
            
            out_line = variant

            if variant in own_var_score_dict:
                own_pred = own_var_score_dict[variant][0]
            else:
                own_pred = 'missing'
            if variant in mvp_var_score_dict:
                mvp_pred = mvp_var_score_dict[variant][0]
            else:
                mvp_pred = 'missing'
            if variant in p2_var_score_dict:
                p2_pred = p2_var_score_dict[variant][0]
            else:
                p2_pred = 'missing'
            if variant in fathmm_var_score_dict:
                fathmm_pred = fathmm_var_score_dict[variant][0]
            else:
                fathmm_pred = 'missing'

            out_line = "%s\t%s\t%s\t%s\t%s\t%s" % (variant, var_label_dict[variant], own_pred, mvp_pred, p2_pred, fathmm_pred)
            file.write(out_line + '\n')
            

    #then go through that file..

    true_labels = []
    predictions_own_test = []
    predictions_mvp_test = []
    predictions_p2_test = []
    predictions_fathmm_test = []

    with open("ALL_test_variants.txt", "r") as file:
        for line in file:
            line = line.strip()
            fields = line.split('\t')
            label = fields[1]
            own_pred = fields[2]
            mvp_pred = fields[3]
            p2_pred = fields[4]
            fathmm_pred = fields[5]

            #only get the ones where we have all
            if not 'missing' in fields:
                if label == 'neutral':
                    true_labels.append(0)
                else:
                    true_labels.append(1)
                own_pred = float(fields[2])
                mvp_pred = float(fields[3])
                p2_pred = float(fields[4])
                fathmm_pred = float(fields[5])
                predictions_own_test.append(own_pred)
                predictions_mvp_test.append(mvp_pred)
                predictions_p2_test.append(p2_pred)
                predictions_fathmm_test.append(fathmm_pred)


    print(len(true_labels))

    if include_all:
        fig_save_path = fig_save_path + '_include_all'

    #true_labels = np.array(true_labels)
    #predictions_own_test = np.array(predictions_own_test)
    #print(true_labels)
    #print(predictions_own_test)

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Compute micro-average ROC curve and ROC area
    fpr["test_mvp"], tpr["test_mvp"], _ = roc_curve(true_labels, predictions_mvp_test, pos_label=1)
    roc_auc["test_mvp"] = auc(fpr["test_mvp"], tpr["test_mvp"])
    fpr["test_PONP2"], tpr["test_PONP2"], _ = roc_curve(true_labels, predictions_p2_test, pos_label=1)
    roc_auc["test_PONP2"] = auc(fpr["test_PONP2"], tpr["test_PONP2"])
    fpr["test_own"], tpr["test_own"], _ = roc_curve(true_labels, predictions_own_test, pos_label=1)
    roc_auc["test_own"] = auc(fpr["test_own"], tpr["test_own"])
    fpr["test_FATHMM"], tpr["test_FATHMM"], _ = roc_curve(true_labels, predictions_fathmm_test, pos_label=0)
    roc_auc["test_FATHMM"] = auc(fpr["test_FATHMM"], tpr["test_FATHMM"])

    plt.figure()
    lw = 2
    plt.plot(fpr["test_mvp"], tpr["test_mvp"],
             lw=lw, label='MVP (area = %0.2f)' % roc_auc["test_mvp"])
    plt.plot(fpr["test_PONP2"], tpr["test_PONP2"],
             lw=lw, label='PONP2 (area = %0.2f)' % roc_auc["test_PONP2"])
    plt.plot(fpr["test_own"], tpr["test_own"],
             lw=lw, label='Pathopred (area = %0.2f)' % roc_auc["test_own"])
    plt.plot(fpr["test_FATHMM"], tpr["test_FATHMM"],
             lw=lw, label='FATHMM (area = %0.2f)' % roc_auc["test_FATHMM"])

    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.grid()
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    if save_fig:
        plt.savefig(fig_save_path + '_ROC.png')

    #also do one with log axis
    plt.figure()
    lw = 2
    plt.semilogx(fpr["test_mvp"], tpr["test_mvp"],
             lw=lw, label='MVP (area = %0.2f)' % roc_auc["test_mvp"])
    plt.semilogx(fpr["test_PONP2"], tpr["test_PONP2"],
             lw=lw, label='PONP2 (area = %0.2f)' % roc_auc["test_PONP2"])
    plt.semilogx(fpr["test_own"], tpr["test_own"],
             lw=lw, label='Pathopred (area = %0.2f)' % roc_auc["test_own"])
    plt.semilogx(fpr["test_FATHMM"], tpr["test_FATHMM"],
             lw=lw, label='FATHMM (area = %0.2f)' % roc_auc["test_FATHMM"])
    
    plt.grid()
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    if save_fig:
        plt.savefig(fig_save_path + '_ROC_log')

    #also get PR curve
    prec = dict()
    rec = dict()
    prec["test_mvp"], rec["test_mvp"], _ = precision_recall_curve(true_labels, predictions_mvp_test, pos_label=1)
    prec["test_own"], rec["test_own"], _ = precision_recall_curve(true_labels, predictions_own_test, pos_label=1)
    prec["test_PONP2"], rec["test_PONP2"], _ = precision_recall_curve(true_labels, predictions_p2_test, pos_label=1)
    prec["test_FATHMM"], rec["test_FATHMM"], _ = precision_recall_curve(true_labels, predictions_fathmm_test, pos_label=0)
    plt.figure()
    lw = 2
    plt.plot(rec["test_mvp"], prec["test_mvp"], 
             lw=lw, label='MVP')
    plt.plot(rec["test_PONP2"], prec["test_PONP2"],
             lw=lw, label='PONP2')
    plt.plot(rec["test_own"], prec["test_own"], 
             lw=lw, label='Pathopred')
    plt.plot(rec["test_FATHMM"], prec["test_FATHMM"],
             lw=lw, label='FATHMM')
    plt.grid()

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PRC')
    plt.legend(loc="lower right")
    if save_fig:
        plt.savefig(fig_save_path + '_PRC')



def compare_classifiers(fig_save_path, include_all, save_fig):

    print('Comparing P2 vs own. Include all: ' + str(include_all))

    #first go through the real labels
    variant_file = "../data/testVariants_p2"

    var_label_dict = {}

    with open(variant_file, 'r') as file:

        for line in file:

            line = line.strip()
            fields = line.split('\t')
            #print(fields)
            variant_line = fields[0] + ' ' + fields[1]
            var_label_dict[variant_line] = fields[2]

    #load the predictor's predictions
    prediction_file = "PON-P2_results/p2_test_predictions.txt"
    p2_var_score_dict = {}
    total_variants = 0
    total_matches = 0
    with open(prediction_file, 'r') as p2_file:

        p2_file.readline()

        for line in p2_file:

            line = line.strip()
            fields = line.split('\t')
            variant_line = fields[0]
            #print(fields)
            label = fields[5].lower()
            score = float(fields[3])
            #print(fields)
            if variant_line not in var_label_dict:
                print(variant_line)
                print('not in dict')
            else:
                #if the label is unknown, we can either skip this variant or include the score for it anyway
                if label == 'unknown':
                    if include_all:
                        pass
                    else:
                        continue

                total_variants += 1
                if var_label_dict[variant_line] == 'neutral':
                    dict_label = 0
                else:
                    dict_label = 1

                p2_var_score_dict[variant_line] = [score,  dict_label]

                if var_label_dict[variant_line] == label:
                    total_matches += 1

    print('Total matches is ' + str(total_matches) + ' out of ' + str(total_variants))

    #also load our own predictions; only those that were predicted on above
    own_pred_file = "own_predictions.txt"
    own_var_score_dict = {}
    total_variants = 0
    total_matches = 0
    with open(own_pred_file, 'r') as own_file:

        for line in own_file:

            line = line.strip()
            fields = line.split('\t')
            seq_name = fields[0]
            mut_str = fields[1]

            variant_line = seq_name + ' ' + mut_str
            label = fields[3].lower()
            score = float(fields[2])

            if variant_line not in var_label_dict:
                print(variant_line)
                print('not in dict')
            elif variant_line not in p2_var_score_dict:
                #print('not predicted by other predictor')
                #print(variant_line)
                continue
            else:

                total_variants += 1
                if var_label_dict[variant_line] == 'neutral':
                    dict_label = 0
                else:
                    dict_label = 1

                own_var_score_dict[variant_line] = [score,  dict_label]

                if var_label_dict[variant_line] == label:
                    total_matches += 1

    print('Total matches is ' + str(total_matches) + ' out of ' + str(total_variants))

    #we do have some variants missing from our own predictions, likely because of error status flag on the sample.
    #let's exclude these ones so we only compare the exact same variants
    to_del = []
    for variant_line in p2_var_score_dict:
        if variant_line not in own_var_score_dict:
            #print(variant_line)
            to_del.append(variant_line)
    for var in to_del:
        del p2_var_score_dict[var]

    print(len(p2_var_score_dict.keys()))
    print(len(own_var_score_dict.keys()))
            

    #get a vector for predicted score and true labels so we can calculate the ROC
    prediction_matrix = np.array(p2_var_score_dict.values())
    predictions_p2_test = prediction_matrix[:,0]
    labels_p2_test = prediction_matrix[:,1]

    prediction_matrix = np.array(own_var_score_dict.values())
    predictions_own_test = prediction_matrix[:,0]
    labels_own_test = prediction_matrix[:,1]

    #print(predictions_p2_test)
    #print(labels_p2_test)

    if include_all:
        fig_save_path = fig_save_path + '_include_all'

    #Get test metrics
    c0_indices = labels_p2_test == 0
    c1_indices = labels_p2_test == 1
    num_neg = sum(c0_indices.astype(int))
    num_pos = sum(c1_indices.astype(int))
    #threshold is 0.5 to assign class
    predictions_copy = (predictions_p2_test > 0.5).astype(int)
    confmat = confusion_matrix(labels_p2_test, predictions_copy)
    TN = confmat[0][0]
    FN = confmat[1][0]
    TP = confmat[1][1]
    FP = confmat[0][1]
    PPV = float(TP)/(TP+FP)
    NPV = float(TN)/(TN+FN)
    sens = float(TP)/(TP+FN)
    spec = float(TN)/(FP+TN)
    acc = float((TP+TN))/(TP+TN+FP+FN)
    mcc = matthews_corrcoef(labels_p2_test, predictions_copy)
    nmcc = float((1+mcc))/2
    OPM = float(( (PPV+NPV)*(sens+spec)*(acc+nmcc) )) / 8
    bacc = ( float(TP) / num_pos + float(TN) / num_neg ) / 2
    accmetrics_test = [PPV,NPV,sens,spec,acc,mcc,OPM, bacc]
    colNames = ['PPV', 'NPV', 'Sens', 'Spec', 'Acc', 'MCC', 'OPM', 'BACC']
    df = pd.DataFrame([accmetrics_test], columns=colNames)
    print("##### On test set : PON-P2 ######")
    print(df)

    #Get test metrics
    c0_indices = labels_own_test == 0
    c1_indices = labels_own_test == 1
    num_neg = sum(c0_indices.astype(int))
    num_pos = sum(c1_indices.astype(int))
    #threshold is 0.5 to assign class
    predictions_copy = (predictions_own_test > 0.5).astype(int)
    confmat = confusion_matrix(labels_own_test, predictions_copy)
    TN = confmat[0][0]
    FN = confmat[1][0]
    TP = confmat[1][1]
    FP = confmat[0][1]
    PPV = float(TP)/(TP+FP)
    NPV = float(TN)/(TN+FN)
    sens = float(TP)/(TP+FN)
    spec = float(TN)/(FP+TN)
    acc = float((TP+TN))/(TP+TN+FP+FN)
    mcc = matthews_corrcoef(labels_own_test, predictions_copy)
    nmcc = float((1+mcc))/2
    OPM = float(( (PPV+NPV)*(sens+spec)*(acc+nmcc) )) / 8
    bacc = ( float(TP) / num_pos + float(TN) / num_neg ) / 2
    accmetrics_test = [PPV,NPV,sens,spec,acc,mcc,OPM, bacc]
    colNames = ['PPV', 'NPV', 'Sens', 'Spec', 'Acc', 'MCC', 'OPM', 'BACC']
    df = pd.DataFrame([accmetrics_test], columns=colNames)
    print("##### On test set : own predictor ######")
    print(df)


    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Compute micro-average ROC curve and ROC area
    fpr["test_ponp2"], tpr["test_ponp2"], _ = roc_curve(labels_p2_test, predictions_p2_test, pos_label=1)
    roc_auc["test_ponp2"] = auc(fpr["test_ponp2"], tpr["test_ponp2"])
    fpr["test_own"], tpr["test_own"], _ = roc_curve(labels_own_test, predictions_own_test, pos_label=1)
    roc_auc["test_own"] = auc(fpr["test_own"], tpr["test_own"])

    plt.figure()
    lw = 2
    plt.plot(fpr["test_ponp2"], tpr["test_ponp2"], color='darkorange',
             lw=lw, label='ROC test_ponp2 (area = %0.2f)' % roc_auc["test_ponp2"])
    plt.plot(fpr["test_own"], tpr["test_own"], color='lightblue',
             lw=lw, label='ROC test_own (area = %0.2f)' % roc_auc["test_own"])

    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.grid()
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    if save_fig:
        plt.savefig(fig_save_path + '_ROC.png')

    #also do one with log axis
    plt.figure()
    lw = 2
    plt.semilogx(fpr["test_ponp2"], tpr["test_ponp2"], color='darkorange',
             lw=lw, label='ROC test_ponp2 (area = %0.2f)' % roc_auc["test_ponp2"])
    plt.semilogx(fpr["test_own"], tpr["test_own"], color='lightblue',
             lw=lw, label='ROC test_own (area = %0.2f)' % roc_auc["test_own"])
    plt.grid()
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    if save_fig:
        plt.savefig(fig_save_path + '_ROC_log')

    #also get PR curve
    prec = dict()
    rec = dict()
    prec["test_ponp2"], rec["test_ponp2"], _ = precision_recall_curve(labels_p2_test, predictions_p2_test, pos_label=1)
    prec["test_own"], rec["test_own"], _ = precision_recall_curve(labels_own_test, predictions_own_test, pos_label=1)
    plt.figure()
    lw = 2
    plt.plot(rec["test_ponp2"], prec["test_ponp2"], color='darkorange',
             lw=lw, label='PR test_ponp2')
    plt.plot(rec["test_own"], prec["test_own"], color='lightblue',
             lw=lw, label='PR test_own')
    plt.grid()

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PRC')
    plt.legend(loc="lower right")
    if save_fig:
        plt.savefig(fig_save_path + '_PRC')


if __name__ == "__main__":

    dataset_name = 'ALL_same'
    include_all = True
    save_fig = True

    #compare_classifiers(dataset_name, include_all, save_fig)
    compare_classifiers_same(dataset_name, include_all, save_fig)
