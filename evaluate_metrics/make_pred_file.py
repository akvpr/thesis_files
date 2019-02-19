import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#don't need the gpu for predictions, save it for other stuff
import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = ""
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.metrics import matthews_corrcoef, confusion_matrix
from tables import * 
from keras.models import load_model

import load_h5


def ref_aa_match(identifier, mut_pos, ref_aa):

    prot_folder = "../proteins_p2/"

    #open the refseq fasta sequence
    with open(prot_folder + identifier + '.fasta', 'r') as prot_file:
        prot_file.readline()
        full_sequence = prot_file.read()
        full_sequence = full_sequence.strip()
        full_sequence = full_sequence.replace('\n', '')

    if len(full_sequence) < mut_pos:
        return False

    return full_sequence[mut_pos-1] == ref_aa


def run_predictions():

    #go through the h5 file and set up own output predictions file
    h5path = "../../p2_presplit.h5"
    h5file = open_file(h5path, mode="r", title="All data")
    check_flag = True
    use_extract = False

    #load our prediction model
    model = load_model('../models/p2_la4_dr7.hdf5')

    prediction_output_file = 'own_predictions.txt'
    output_handle = open(prediction_output_file, 'w')

    category = 'test'

    allnodes = h5file.list_nodes('/' + category)
    print('Reading data for category ' + category)
    i = 0
    for node in allnodes:
        
        stable = node.properties
        
        if check_flag:
            #first check if this sample is ok or error
            status = stable.col('statusFlag')
            if status[0] == 'error':
                continue

        mut_str = stable.col('mutStr')[0]
        seq_name = stable.col('seqName')[0]
        if not ref_aa_match( seq_name, int(mut_str[1:-1]), mut_str[0] ):
            print(seq_name + ' and variant position does not match')
            continue

        refalt = np.array(node.refalt.read())
        freq = np.array(node.freqs.read())
        sinfo = np.array(node.sinfo.read())
        pentr = np.array(node.pentr.read())
        seq = np.array(node.seq.read())

        if use_extract:
            ss = np.array(node.ss_extract.read())
            if ss.shape[0] < 21:
                to_pad = 21 - ss.shape[0]
                ss = np.pad(ss, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        else:
            ss = np.array(node.ssmatrix.read())

        #expand to prepare for prediction
        refalt = np.expand_dims(refalt,axis=0)
        freq = np.expand_dims(freq, axis=0)
        sinfo = np.expand_dims(sinfo, axis=0)
        pentr = np.expand_dims(pentr, axis=0)
        seq = np.expand_dims(seq, axis=0)
        ss = np.expand_dims(ss, axis=0)

        labelText = stable.col('label')
        labelText = labelText[0]
        seqname = stable.col('seqName')
        if labelText == "neutral":
                label = 0
        elif labelText == "pathogenic":
                label = 1

        seq_name = stable.col('seqName')[0]
        mut_str = stable.col('mutStr')[0]
        #print(seq_name + ' - ' + mut_str)

        indata = [refalt, freq, sinfo, pentr, seq, ss]

        prediction = model.predict(indata, batch_size=1)[0][0]
        if prediction < 0.5:
            pred_label = 'neutral'
        else:
            pred_label = 'pathogenic'
            
        output_line = seq_name + '\t' + mut_str + '\t' + str(prediction) + '\t' + pred_label
        print(output_line)
        output_handle.write(output_line + '\n')
        




def run_predictions_sequence():

    #go through the h5 file and set up own output predictions file
    h5path = "/scratch_ssd/Data/p2_sequence.h5"
    h5file = open_file(h5path, mode="r", title="All data")
    check_flag = True
    use_extract = False

    #load our prediction model
    model = load_model('../models/p2_sequence.h5')

    prediction_output_file = 'own_predictions_sequence.txt'
    output_handle = open(prediction_output_file, 'w')

    category = 'test'

    allnodes = h5file.list_nodes('/' + category)
    print('Reading data for category ' + category)
    i = 0
    for node in allnodes:
        
        stable = node.properties

        #first check so we have a match
        mut_str = stable.col('mutStr')[0]
        seq_name = stable.col('seqName')[0]
        if not ref_aa_match( seq_name, int(mut_str[1:-1]), mut_str[0] ):
            print(seq_name + ' and variant position does not match')
            continue
        
        if check_flag:
            #first check if this sample is ok or error
            status = stable.col('statusFlag')
            if status[0] == 'error':
                continue


        refalt = np.array(node.refalt.read())
        blosum = np.array(node.blosum.read())
        seq = np.array(node.seq.read())
        ss = np.array(node.ssmatrix.read())

        mut_pos = int(mut_str[1:-1])
        seq = seq[mut_pos-11:mut_pos+10, :]

        if seq.shape[0] < 21:
            to_pad = 21 - seq.shape[0]
            seq = np.pad(seq, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))


        #expand to prepare for prediction
        blosum = np.expand_dims(blosum, axis=0)
        refalt = np.expand_dims(refalt,axis=0)
        seq = np.expand_dims(seq, axis=0)
        ss = np.expand_dims(ss, axis=0)

        labelText = stable.col('label')
        labelText = labelText[0]
        seqname = stable.col('seqName')
        if labelText == "neutral":
                label = 0
        elif labelText == "pathogenic":
                label = 1

        seq_name = stable.col('seqName')[0]
        mut_str = stable.col('mutStr')[0]
        #print(seq_name + ' - ' + mut_str)

        indata = [refalt, blosum, seq, ss]

        prediction = model.predict(indata, batch_size=1)[0][0]
        if prediction < 0.5:
            pred_label = 'neutral'
        else:
            pred_label = 'pathogenic'
            
        output_line = seq_name + '\t' + mut_str + '\t' + str(prediction) + '\t' + pred_label
        print(output_line)
        output_handle.write(output_line + '\n')
        



def evaluate_predictions():

    h5path = "../../p2_presplit.h5"
    h5file = open_file(h5path, mode="r", title="All data")
    check_flag = True
    use_extract = True
    fig_save_path = "newtest"

    #load our prediction model
    model = load_model('../models/p2_la4_dr7.hdf5')

    indata_train, indata_val, indata_test, train_labels, val_labels, test_labels = load_h5.load_data_presplit_components(h5path, True)
    #indata_train, indata_val, indata_test, train_labels, val_labels, test_labels = load_h5.load_data_presplit_seqver(h5path, True)
    
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    c0_indices = train_labels == 0
    c1_indices = train_labels == 1
    num_neg = sum(c0_indices.astype(int))
    num_pos = sum(c1_indices.astype(int))

    predictions_train = model.predict(indata_train, batch_size=128)
    #threshold is 0.5 to assign class
    predictions_copy = (predictions_train > 0.5).astype(int)
    confmat = confusion_matrix(train_labels, predictions_copy)

    TN = confmat[0][0]
    FN = confmat[1][0]
    TP = confmat[1][1]
    FP = confmat[0][1]
    #or..
    #TN, FP, FN, TP = confmat.ravel()
    PPV = float(TP)/(TP+FP)
    NPV = float(TN)/(TN+FN)
    sens = float(TP)/(TP+FN)
    spec = float(TN)/(FP+TN)
    acc = float((TP+TN))/(TP+TN+FP+FN)
    mcc = matthews_corrcoef(train_labels, predictions_copy)
    nmcc = float((1+mcc))/2
    OPM = float(( (PPV+NPV)*(sens+spec)*(acc+nmcc) )) / 8
    bacc = ( float(TP) / num_pos + float(TN) / num_neg ) / 2
    fpr["train"], tpr["train"], _ = roc_curve(train_labels, predictions_train, pos_label=1)
    roc_auc["train"] = auc(fpr["train"], tpr["train"])

    accmetrics_train = [PPV,NPV,sens,spec,acc,mcc,OPM, bacc, roc_auc['train']]
    colNames = ['PPV', 'NPV', 'Sens', 'Spec', 'Acc', 'MCC', 'OPM', 'BACC', 'AUC']
    df = pd.DataFrame([accmetrics_train], columns=colNames)
    print("##### On train set ######")
    print(df)

    c0_indices = val_labels == 0
    c1_indices = val_labels == 1
    num_neg = sum(c0_indices.astype(int))
    num_pos = sum(c1_indices.astype(int))

    predictions_val = model.predict(indata_val, batch_size=128)
    #threshold is 0.5 to assign class
    predictions_copy = (predictions_val > 0.5).astype(int)
    confmat = confusion_matrix(val_labels, predictions_copy)
    TN = confmat[0][0]
    FN = confmat[1][0]
    TP = confmat[1][1]
    FP = confmat[0][1]
    PPV = float(TP)/(TP+FP)
    NPV = float(TN)/(TN+FN)
    sens = float(TP)/(TP+FN)
    spec = float(TN)/(FP+TN)
    acc = float((TP+TN))/(TP+TN+FP+FN)
    mcc = matthews_corrcoef(val_labels, predictions_copy)
    nmcc = float((1+mcc))/2
    OPM = float(( (PPV+NPV)*(sens+spec)*(acc+nmcc) )) / 8
    bacc = ( float(TP) / num_pos + float(TN) / num_neg ) / 2

    fpr["val"], tpr["val"], _ = roc_curve(val_labels, predictions_val, pos_label=1)
    roc_auc["val"] = auc(fpr["val"], tpr["val"])
    accmetrics_train = [PPV,NPV,sens,spec,acc,mcc,OPM, bacc, roc_auc['val']]
    colNames = ['PPV', 'NPV', 'Sens', 'Spec', 'Acc', 'MCC', 'OPM', 'BACC', 'AUC']
    df = pd.DataFrame([accmetrics_train], columns=colNames)
    print("##### On val set ######")
    print(df)

    c0_indices = test_labels == 0
    c1_indices = test_labels == 1
    num_neg = sum(c0_indices.astype(int))
    num_pos = sum(c1_indices.astype(int))
    
    predictions_test = model.predict(indata_test, batch_size=128)
    #threshold is 0.5 to assign class
    predictions_copy = (predictions_test > 0.5).astype(int)
    confmat = confusion_matrix(test_labels, predictions_copy)
    TN = confmat[0][0]
    FN = confmat[1][0]
    TP = confmat[1][1]
    FP = confmat[0][1]
    PPV = float(TP)/(TP+FP)
    NPV = float(TN)/(TN+FN)
    sens = float(TP)/(TP+FN)
    spec = float(TN)/(FP+TN)
    acc = float((TP+TN))/(TP+TN+FP+FN)
    mcc = matthews_corrcoef(test_labels, predictions_copy)
    nmcc = float((1+mcc))/2
    OPM = float(( (PPV+NPV)*(sens+spec)*(acc+nmcc) )) / 8
    bacc = ( float(TP) / num_pos + float(TN) / num_neg ) / 2

    fpr["test"], tpr["test"], _ = roc_curve(test_labels, predictions_test, pos_label=1)
    roc_auc["test"] = auc(fpr["test"], tpr["test"])
    accmetrics_train = [PPV,NPV,sens,spec,acc,mcc,OPM, bacc, roc_auc['test']]
    colNames = ['PPV', 'NPV', 'Sens', 'Spec', 'Acc', 'MCC', 'OPM', 'BACC', 'AUC']
    df = pd.DataFrame([accmetrics_train], columns=colNames)
    print("##### On test set ######")
    print(df)

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Compute micro-average ROC curve and ROC area
    fpr["val"], tpr["val"], _ = roc_curve(val_labels, predictions_val, pos_label=1)
    roc_auc["val"] = auc(fpr["val"], tpr["val"])
    fpr["test"], tpr["test"], _ = roc_curve(test_labels, predictions_test, pos_label=1)
    roc_auc["test"] = auc(fpr["test"], tpr["test"])

    plt.figure()
    lw = 2
    plt.plot(fpr["test"], tpr["test"], color='darkorange',
             lw=lw, label='ROC test (area = %0.2f)' % roc_auc["test"])
    plt.plot(fpr["val"], tpr["val"], color='gold',
            lw=lw, label='ROC val (area = %0.2f)' % roc_auc["val"])

    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.grid()
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    plt.savefig(fig_save_path + '_ROC.png')

    #also do one with log axis
    plt.figure()
    lw = 2
    plt.semilogx(fpr["test"], tpr["test"], color='darkorange',
             lw=lw, label='ROC test (area = %0.2f)' % roc_auc["test"])
    plt.semilogx(fpr["val"], tpr["val"], color='gold',
            lw=lw, label='ROC val (area = %0.2f)' % roc_auc["val"])
    #plt.semilogx([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.grid()
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    plt.savefig(fig_save_path + '_ROC_log')

    #also get PR curve
    prec = dict()
    rec = dict()
    prec["val"], rec["val"], _ = precision_recall_curve(val_labels, predictions_val, pos_label=1)
    prec["test"], rec["test"], _ = precision_recall_curve(test_labels, predictions_test, pos_label=1)
    plt.figure()
    lw = 2
    plt.plot(rec["test"], prec["test"], color='darkorange',
             lw=lw, label='PR test')
    plt.plot(rec["val"], prec["val"], color='gold',
            lw=lw, label='PR val')
    plt.grid()

    #plt.plot([1, 0], [1, 0], color='navy', lw=lw, linestyle='--')
    #plt.xlim([0.0, 1.0])
    #plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PRC')
    plt.legend(loc="lower right")
    plt.savefig(fig_save_path + '_PRC')



if __name__ == "__main__":

    run_predictions()
    #run_predictions_sequence()
    #evaluate_predictions()
