import numpy as np
from tables import * 
from sklearn.model_selection import train_test_split

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



def load_data_category_component(h5file, category, check_flag, use_extract):


    #load train data
    cat_refalt = []
    cat_freq = []
    cat_sinfo = []
    cat_pentr = []
    cat_seq = []
    cat_ss = []
    cat_labels = []

    allnodes = h5file.list_nodes('/' + category)
    print('Reading data for category ' + category)
    i = 0
    for node in allnodes:

        stable = node.properties

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


        labelText = stable.col('label')
        labelText = labelText[0]
        if labelText == "neutral":
                cat_labels.append(0)
        elif labelText == "pathogenic": 
                cat_labels.append(1)

        #freq = np.expand_dims(freq, axis=-1)
        #sinfo = np.expand_dims(sinfo, axis=-1)
        #pentr = np.expand_dims(pentr, axis=-1)
        #seq = np.expand_dims(seq, axis=-1)
        #ss = np.expand_dims(ss, axis=-1)

        cat_refalt.append(refalt)
        cat_freq.append(freq)
        cat_sinfo.append(sinfo)
        cat_pentr.append(pentr)
        cat_seq.append(seq)
        cat_ss.append(ss)

    cat_refalt = np.array(cat_refalt)
    cat_freq = np.array(cat_freq)
    cat_sinfo = np.array(cat_sinfo)
    cat_pentr = np.array(cat_pentr)
    cat_seq = np.array(cat_seq)
    cat_ss = np.array(cat_ss)
    cat_labels = np.array(cat_labels)

    print(cat_refalt.shape)
    print(cat_freq.shape)
    print(cat_sinfo.shape)
    print(cat_pentr.shape)
    print(cat_seq.shape)
    print(cat_ss.shape)
    print(cat_labels.shape)

    print('Done reading data for category ' + category)

    numPath = (cat_labels == 1).sum()
    numNeut = (cat_labels == 0).sum()
    print('Pathogenic : ' + str(numPath))
    print('Neutral : ' + str(numNeut))

    indata_cat = [cat_refalt, cat_freq, cat_sinfo, cat_pentr, cat_seq, cat_ss]

    return indata_cat, cat_labels


def load_data_category_feat(h5file, category, check_flag):

    #load train data (this one has only train data; split this up further below..)
    cat_ss3 = []
    cat_ss6 = []
    cat_rsa = []
    cat_dihed = []
    cat_feat = []
    cat_labels = []
    print('Reading all data')
    allnodes = h5file.list_nodes('/training')
    i = 0
    for node in allnodes:

        stable = node.properties

        if check_flag:
            #first check if this sample is ok or error
            status = stable.col('statusFlag')
            #print(status)
            if status[1] == 'error':
                continue
        
        ss3 = np.array(node.ss3.read())
        ss6 = np.array(node.ss6.read())
        rsa = np.array(node.rsa.read())
        dihed = np.array(node.dihed.read())
        feat = np.array(node.feat.read())
            
        cat_ss3.append(ss3)
        cat_ss6.append(ss6)
        cat_rsa.append(rsa)
        cat_dihed.append(dihed)
        cat_feat.append(feat)
        
        labelText = stable.col('label')
        #one sample has a double insertion here
        labelText = labelText[0]
        seqname = stable.col('seqName')
        if labelText == "neutral":
                cat_labels.append(0)
        elif labelText == "pathogenic":
                cat_labels.append(1)
        
        i += 1

    cat_ss3 = np.array(cat_ss3)
    cat_ss6 = np.array(cat_ss6)
    cat_rsa = np.array(cat_rsa)
    cat_dihed = np.array(cat_dihed)
    cat_feat = np.array(cat_feat)
    cat_labels = np.array(cat_labels)

    print('Done reading all data')

    h5file.close()

    numPath = (cat_labels == 1).sum()
    numNeut = (cat_labels == 0).sum()
    print('Pathogenic : ' + str(numPath))
    print('Neutral : ' + str(numNeut))

    #indata_cat = [cat_ss3, cat_ss6, cat_rsa, cat_dihed, cat_feat]
    return cat_ss3, cat_ss6, cat_rsa, cat_dihed, cat_feat, cat_labels



def load_data_presplit_components(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    indata_train, train_labels = load_data_category_component(h5file, 'training', check_flag, False)
    indata_val, val_labels = load_data_category_component(h5file, 'validation', check_flag, False)
    indata_test, test_labels = load_data_category_component(h5file, 'test', check_flag, False)

    return indata_train, indata_val, indata_test, train_labels, val_labels, test_labels



def load_data_presplit_shuffle(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    indata_train, train_labels = load_data_category_component(h5file, 'training', check_flag, False)
    indata_val, val_labels = load_data_category_component(h5file, 'validation', check_flag, False)
    indata_test, test_labels = load_data_category_component(h5file, 'test', check_flag, False)

    #concatenate val and test, then split again....
    pool_refalt = np.concatenate( (indata_val[0], indata_test[0]), axis=0 )
    pool_freq = np.concatenate( (indata_val[1], indata_test[1]), axis=0 )
    pool_sinfo = np.concatenate( (indata_val[2], indata_test[2]), axis=0 )
    pool_pentr = np.concatenate( (indata_val[3], indata_test[3]), axis=0 )
    pool_seq = np.concatenate( (indata_val[4], indata_test[4]), axis=0 )
    pool_ss = np.concatenate( (indata_val[5], indata_test[5]), axis=0 )
    pool_labels = np.concatenate( (val_labels, test_labels), axis=0)

    print('All')
    print(pool_refalt.shape)

    refalt_s_val, refalt_s_test, freq_s_val, freq_s_test, sinfo_s_val, sinfo_s_test, pentr_s_val, pentr_s_test, seq_s_val, seq_s_test, ss_s_val, ss_s_test, labels_s_val, labels_s_test = train_test_split(pool_refalt, pool_freq, pool_sinfo, pool_pentr, pool_seq, pool_ss, pool_labels, test_size=0.5, random_state=0)

    print('\nSplit')
    print(refalt_s_val.shape)
    print(refalt_s_test.shape)

    indata_val = [refalt_s_val, freq_s_val, sinfo_s_val, pentr_s_val, seq_s_val, ss_s_val]
    indata_test = [refalt_s_test, freq_s_test, sinfo_s_test, pentr_s_test, seq_s_test, ss_s_test]


    return indata_train, indata_val, indata_test, train_labels, labels_s_val, labels_s_test




def load_data_all(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    all_ss3, all_ss6, all_rsa, all_dihed, all_feat, all_labels = load_data_category_feat(h5file, 'training', check_flag)

    #then we split it up
    print('Splitting data')

    print('All')
    print(all_feat.shape)

    ss3_s_trainval, ss3_s_test, ss6_s_trainval, ss6_s_test, rsa_s_trainval, rsa_s_test, dihed_s_trainval, dihed_s_test, feat_s_trainval, feat_s_test, labels_s_trainval, labels_s_test = train_test_split(all_ss3, all_ss6, all_rsa, all_dihed, all_feat, all_labels, test_size=0.2, random_state=0, stratify=all_labels)

    print('\nSplit')
    print(feat_s_trainval.shape)
    print(feat_s_test.shape)

    #further split training into 10% validation set as well
    ss3_s_train, ss3_s_val, ss6_s_train, ss6_s_val, rsa_s_train, rsa_s_val, dihed_s_train, dihed_s_val, feat_s_train, feat_s_val, labels_s_train, labels_s_val = train_test_split(ss3_s_trainval, ss6_s_trainval, rsa_s_trainval, dihed_s_trainval, feat_s_trainval, labels_s_trainval, test_size=0.2, random_state=0, stratify=labels_s_trainval)
    print('\nSplit again')
    print(feat_s_train.shape)
    print(feat_s_val.shape)

    numPath_train = (labels_s_train == 1).sum()
    numNeut_train = (labels_s_train == 0).sum()
    numPath_val = (labels_s_val == 1).sum()
    numNeut_val = (labels_s_val == 0).sum()
    numPath_test = (labels_s_test == 1).sum()
    numNeut_test = (labels_s_test == 0).sum()
    print('Pathogenic train: ' + str(numPath_train))
    print('Neutral train: ' + str(numNeut_train))
    print('Pathogenic val: ' + str(numPath_val))
    print('Neutral val: ' + str(numNeut_val))
    print('Pathogenic test: ' + str(numPath_test))
    print('Neutral test: ' + str(numNeut_test))

    indata_train = [ss3_s_train, ss6_s_train, rsa_s_train, dihed_s_train, feat_s_train]
    indata_val = [ss3_s_val, ss6_s_val, rsa_s_val, dihed_s_val, feat_s_val]
    indata_test = [ss3_s_test, ss6_s_test, rsa_s_test, dihed_s_test, feat_s_test]

    #indata_train = np.concatenate((ss3_s_train, ss6_s_train, rsa_s_train, dihed_s_train, feat_s_train), axis = 2)
    #indata_val = np.concatenate((ss3_s_val, ss6_s_val, rsa_s_val, dihed_s_val, feat_s_val), axis = 2)
    #indata_test = np.concatenate((ss3_s_test, ss6_s_test, rsa_s_test, dihed_s_test, feat_s_test), axis = 2)

    return indata_train, indata_val, indata_test, labels_s_train, labels_s_val, labels_s_test

def load_data_train_test(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    train_ss3, train_ss6, train_rsa, train_dihed, train_feat, train_labels = load_data_category_feat(h5file, 'training', check_flag)
    test_ss3, test_ss6, test_rsa, test_dihed, test_feat, test_labels = load_data_category_feat(h5file, 'test', check_flag)

    h5file.close()

    numPath = (train_labels == 1).sum()
    numNeut = (train_labels == 0).sum()
    print('Pathogenic train: ' + str(numPath))
    print('Neutral train: ' + str(numNeut))
    numPath = (test_labels == 1).sum()
    numNeut = (test_labels == 0).sum()
    print('Pathogenic test: ' + str(numPath))
    print('Neutral test: ' + str(numNeut))

    print('All')
    print(train_feat.shape)

    ss3_s_train, ss3_s_val, ss6_s_train, ss6_s_val, rsa_s_train, rsa_s_val, dihed_s_train, dihed_s_val, feat_s_train, feat_s_val, labels_s_train, labels_s_val = train_test_split(train_ss3, train_ss6, train_rsa, train_dihed, train_feat, train_labels, test_size=0.2, random_state=0, stratify=train_labels)

    print('\nSplit')
    print(feat_s_train.shape)
    print(feat_s_val.shape)

    indata_train = [ss3_s_train, ss6_s_train, rsa_s_train, dihed_s_train, feat_s_train]
    indata_val = [ss3_s_val, ss6_s_val, rsa_s_val, dihed_s_val, feat_s_val]
    indata_test = [test_ss3, test_ss6, test_rsa, test_dihed, test_feat]

    #indata_train = np.concatenate((ss3_s_train, ss6_s_train, rsa_s_train, dihed_s_train, feat_s_train), axis = 2)
    #indata_val = np.concatenate((ss3_s_val, ss6_s_val, rsa_s_val, dihed_s_val, feat_s_val), axis = 2)
    #indata_test = np.concatenate((test_ss3, test_ss6, test_rsa, test_dihed, test_feat), axis = 2)
    #print(indata_train.shape)
    #print(indata_val.shape)
    #print(indata_test.shape)

    return indata_train, indata_val, indata_test, labels_s_train, labels_s_val, test_labels





def load_data_train_test_components(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    train_data_vec, train_labels = load_data_category_component(h5file, 'training', check_flag, False)
    test_data_vec, test_labels = load_data_category_component(h5file, 'test', check_flag, False)

    h5file.close()

    train_refalt = train_data_vec[0]
    train_freq = train_data_vec[1]
    train_sinfo = train_data_vec[2]
    train_pentr = train_data_vec[3] 
    train_seq = train_data_vec[4]
    train_ss = train_data_vec[5]

    print('All')
    print(train_refalt.shape)

    refalt_s_train, refalt_s_val, freq_s_train, freq_s_val, sinfo_s_train, sinfo_s_val, pentr_s_train, pentr_s_val, seq_s_train, seq_s_val, ss_s_train, ss_s_val, labels_s_train, labels_s_val = train_test_split(train_refalt, train_freq, train_sinfo, train_pentr, train_seq, train_ss, train_labels, test_size=0.2, random_state=0, stratify=train_labels)

    print('\nSplit')
    print(refalt_s_train.shape)
    print(refalt_s_val.shape)

    indata_train = [refalt_s_train, freq_s_train, sinfo_s_train, pentr_s_train, seq_s_train, ss_s_train]
    indata_val = [refalt_s_val, freq_s_val, sinfo_s_val, pentr_s_val, seq_s_val, ss_s_val]


    return indata_train, indata_val, test_data_vec, labels_s_train, labels_s_val, test_labels



def load_data_train_test_mod7(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    train_data_vec, train_labels = load_data_category_component(h5file, 'training', check_flag, True)
    test_data_vec, test_labels = load_data_category_component(h5file, 'test', check_flag, True)

    h5file.close()

    train_refalt = train_data_vec[0]
    train_freq = train_data_vec[1]
    train_sinfo = train_data_vec[2]
    train_pentr = train_data_vec[3] 
    train_seq = train_data_vec[4]
    train_ss = train_data_vec[5]

    print('All')
    print(train_refalt.shape)

    refalt_s_train, refalt_s_val, freq_s_train, freq_s_val, sinfo_s_train, sinfo_s_val, pentr_s_train, pentr_s_val, seq_s_train, seq_s_val, ss_s_train, ss_s_val, labels_s_train, labels_s_val = train_test_split(train_refalt, train_freq, train_sinfo, train_pentr, train_seq, train_ss, train_labels, test_size=0.2, random_state=0, stratify=train_labels)

    print('\nSplit')
    print(refalt_s_train.shape)
    print(refalt_s_val.shape)

    indata_train = [refalt_s_train, freq_s_train, sinfo_s_train, pentr_s_train, seq_s_train, ss_s_train]
    indata_val = [refalt_s_val, freq_s_val, sinfo_s_val, pentr_s_val, seq_s_val, ss_s_val]
    indata_test = [test_refalt, test_freq, test_sinfo, test_pentr, test_seq, test_ss]


    return indata_train, indata_val, indata_test, labels_s_train, labels_s_val, test_labels



def load_data_all_components(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    train_data_vec, train_labels = load_data_category_component(h5file, 'training', check_flag, False)

    train_refalt = train_data_vec[0]
    train_freq = train_data_vec[1]
    train_sinfo = train_data_vec[2]
    train_pentr = train_data_vec[3] 
    train_seq = train_data_vec[4]
    train_ss = train_data_vec[5]

    h5file.close()

    numPath = (train_labels == 1).sum()
    numNeut = (train_labels == 0).sum()
    print('Pathogenic train: ' + str(numPath))
    print('Neutral train: ' + str(numNeut))

    train_refalt = train_data_vec[0]
    train_freq = train_data_vec[1]
    train_sinfo = train_data_vec[2]
    train_pentr = train_data_vec[3] 
    train_seq = train_data_vec[4]
    train_ss = train_data_vec[5]

    print('All')
    print(train_refalt.shape)

    refalt_s_train, refalt_s_val, freq_s_train, freq_s_val, sinfo_s_train, sinfo_s_val, pentr_s_train, pentr_s_val, seq_s_train, seq_s_val, ss_s_train, ss_s_val, labels_s_train, labels_s_val = train_test_split(train_refalt, train_freq, train_sinfo, train_pentr, train_seq, train_ss, train_labels, test_size=0.2, random_state=0, stratify=train_labels)

    print('\nSplit')
    print(refalt_s_train.shape)
    print(refalt_s_val.shape)

    indata_train = [train_refalt, train_freq, train_sinfo, train_pentr, train_seq, train_ss]
    #labels_s_val = []
    #test_labels = []
    indata_val = [refalt_s_val, freq_s_val, sinfo_s_val, pentr_s_val, seq_s_val, ss_s_val]
    indata_test = [test_refalt, test_freq, test_sinfo, test_pentr, test_seq, test_ss]
    #indata_val = []
    #indata_test = []


    return indata_train, indata_val, indata_test, train_labels, labels_s_val, test_labels


######################



def load_data_train_test_mod6(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    #load train data

    train_refalt = []
    train_freq = []
    train_sinfo = []
    train_pentr = []
    train_seq = []

    train_ss_seq = []
    train_ss_sinfo = []
    train_ss_pentr = []

    train_labels = []

    allnodes = h5file.list_nodes('/training')
    print('Reading train data')
    i = 0
    for node in allnodes:

        if i > 1000:
            continue

        i += 1
        
        stable = node.properties
        
        if check_flag:
            #first check if this sample is ok or error
            status = stable.col('statusFlag')
            if status[1] == 'error':
                continue

        refalt = np.array(node.refalt.read())
        refalt = np.expand_dims(refalt,axis=0)
        freq = np.array(node.freqs.read())
        sinfo = np.array(node.sinfo.read())
        pentr = np.array(node.pentr.read())
        seq = np.array(node.seq.read())

        #get the mutpos to get around position to test
        mut_str = stable.col('mutStr')[0]
        mut_pos = int(mut_str[1:-1])

        ss_seq = np.array(node.ss_seq.read())
        ss_sinfo = np.array(node.ss_sinfo.read())
        ss_pentr = np.array(node.ss_pentr.read())

        ss_seq = ss_seq[mut_pos-11:mut_pos+10,:]
        ss_sinfo = ss_sinfo[mut_pos-11:mut_pos+10,:]
        ss_pentr = ss_pentr[mut_pos-11:mut_pos+10,:]

        if ss_seq.shape[0] < 21:
            to_pad = 21 - ss_seq.shape[0]
            ss_seq = np.pad(ss_seq, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
            ss_sinfo = np.pad(ss_sinfo, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
            ss_pentr = np.pad(ss_pentr, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))


        labelText = stable.col('label')
        #one sample has a double insertion here
        labelText = labelText[0]
        seqname = stable.col('seqName')
        if labelText == "neutral":
                train_labels.append(0)
        elif labelText == "pathogenic":
                train_labels.append(1)

        train_refalt.append(refalt)
        train_freq.append(freq)
        train_sinfo.append(sinfo)
        train_pentr.append(pentr)
        train_seq.append(seq)

        train_ss_seq.append(ss_seq)
        train_ss_sinfo.append(ss_sinfo)
        train_ss_pentr.append(ss_pentr)


    train_refalt = np.array(train_refalt)
    train_freq = np.array(train_freq)
    train_sinfo = np.array(train_sinfo)
    train_pentr = np.array(train_pentr)
    train_seq = np.array(train_seq)

    

    train_ss_seq = np.array(train_ss_seq)
    train_ss_sinfo = np.array(train_ss_sinfo)
    train_ss_pentr = np.array(train_ss_pentr)

    train_labels = np.array(train_labels)

    print(train_refalt.shape)
    print(train_freq.shape)
    print(train_sinfo.shape)
    print(train_pentr.shape)
    print(train_seq.shape)

    print(train_ss_seq.shape)
    print(train_ss_sinfo.shape)
    print(train_ss_pentr.shape)

    print(train_labels.shape)

    print('Done reading training data')

    #load test data
    test_refalt = []
    test_freq = []
    test_sinfo = []
    test_pentr = []
    test_seq = []

    test_ss_seq = []
    test_ss_sinfo = []
    test_ss_pentr = []

    test_labels = []

    allnodes = h5file.list_nodes('/test')
    print('Reading test data')
    for node in allnodes:

        stable = node.properties
        
        if check_flag:
            #first check if this sample is ok or error
            status = stable.col('statusFlag')
            if status[1] == 'error':
                continue
        
        labelText = stable.col('label')
        #one sample has a double insertion here
        labelText = labelText[0]
        if labelText == "neutral":
                test_labels.append(0)
        elif labelText == "pathogenic":
                test_labels.append(1)
        
        refalt = np.array(node.refalt.read())
        refalt = np.expand_dims(refalt,axis=0)
        freq = np.array(node.freqs.read())
        sinfo = np.array(node.sinfo.read())
        pentr = np.array(node.pentr.read())
        seq = np.array(node.seq.read())

        #get the mutpos to get around position to test
        mut_str = stable.col('mutStr')[0]
        mut_pos = int(mut_str[1:-1])

        ss_seq = np.array(node.ss_seq.read())
        ss_sinfo = np.array(node.ss_sinfo.read())
        ss_pentr = np.array(node.ss_pentr.read())

        ss_seq = ss_seq[mut_pos-11:mut_pos+10]
        ss_sinfo = ss_sinfo[mut_pos-11:mut_pos+10]
        ss_pentr = ss_pentr[mut_pos-11:mut_pos+10]

        if ss_seq.shape[0] < 21:
            to_pad = 21 - ss_seq.shape[0]
            ss_seq = np.pad(ss_seq, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
            ss_sinfo = np.pad(ss_sinfo, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
            ss_pentr = np.pad(ss_pentr, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))



        test_refalt.append(refalt)
        test_freq.append(freq)
        test_sinfo.append(sinfo)
        test_pentr.append(pentr)
        test_seq.append(seq)

        test_ss_seq.append(ss_seq)
        test_ss_sinfo.append(ss_sinfo)
        test_ss_pentr.append(ss_pentr)

        
    test_refalt = np.array(test_refalt)
    test_freq = np.array(test_freq)
    test_sinfo = np.array(test_sinfo)
    test_pentr = np.array(test_pentr)
    test_seq = np.array(test_seq)

    test_ss_seq = np.array(test_ss_seq)
    test_ss_sinfo = np.array(test_ss_sinfo)
    test_ss_pentr = np.array(test_ss_pentr)

    test_labels = np.array(test_labels)

    print(test_refalt.shape)
    print(test_freq.shape)
    print(test_sinfo.shape)
    print(test_pentr.shape)
    print(test_seq.shape)

    print(test_ss_seq.shape)
    print(test_ss_sinfo.shape)
    print(test_ss_pentr.shape)

    print(test_labels.shape)

    print('Done reading test data')

    h5file.close()


    numPath = (train_labels == 1).sum()
    numNeut = (train_labels == 0).sum()
    print('Pathogenic train: ' + str(numPath))
    print('Neutral train: ' + str(numNeut))
    numPath = (test_labels == 1).sum()
    numNeut = (test_labels == 0).sum()
    print('Pathogenic test: ' + str(numPath))
    print('Neutral test: ' + str(numNeut))

    print('All')
    print(train_refalt.shape)

    refalt_s_train, refalt_s_val, freq_s_train, freq_s_val, sinfo_s_train, sinfo_s_val, pentr_s_train, pentr_s_val, seq_s_train, seq_s_val, \
        ss_seq_s_train, ss_seq_s_val, ss_sinfo_s_train, ss_sinfo_s_val, ss_pentr_s_train, ss_pentr_s_val, labels_s_train, labels_s_val = train_test_split(train_refalt, train_freq, train_sinfo, train_pentr, train_seq, 
            train_ss_seq, train_ss_sinfo, train_ss_pentr, train_labels, test_size=0.2, random_state=0, stratify=train_labels)

    print('\nSplit')
    print(refalt_s_train.shape)
    print(refalt_s_val.shape)

    #expand labels??
    labels_s_train = np.expand_dims(labels_s_train, axis=1)
    labels_s_train = np.expand_dims(labels_s_train, axis=1)
    labels_s_val = np.expand_dims(labels_s_val, axis=1)
    labels_s_val = np.expand_dims(labels_s_val, axis=1)
    test_labels = np.expand_dims(test_labels, axis=1)
    test_labels = np.expand_dims(test_labels, axis=1)

    indata_train = [refalt_s_train, freq_s_train, sinfo_s_train, pentr_s_train, seq_s_train, ss_seq_s_train, ss_sinfo_s_train, ss_pentr_s_train]
    indata_val = [refalt_s_val, freq_s_val, sinfo_s_val, pentr_s_val, seq_s_val, ss_seq_s_val, ss_sinfo_s_val, ss_pentr_s_val]
    indata_test = [test_refalt, test_freq, test_sinfo, test_pentr, test_seq, test_ss_seq, test_ss_sinfo, test_ss_pentr]


    return indata_train, indata_val, indata_test, labels_s_train, labels_s_val, test_labels












def load_data_category_seqver(h5file, category, check_flag, use_extract):


    #load train data
    cat_refalt = []
    cat_blosum = []
    cat_seq = []
    cat_ss = []
    cat_labels = []

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
        #refalt[0, -2] = 0
        #refalt = np.expand_dims(refalt,axis=0)
        blosum = np.array(node.blosum.read())
        seq = np.array(node.seq.read())
        ss = np.array(node.ssmatrix.read())

        mut_pos = int(mut_str[1:-1])
        seq = seq[mut_pos-11:mut_pos+10, :]

        if seq.shape[0] < 21:
            to_pad = 21 - seq.shape[0]
            seq = np.pad(seq, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))

        #print(seq.shape)

        labelText = stable.col('label')
        labelText = labelText[0]
        if labelText == "neutral":
                cat_labels.append(0)
        elif labelText == "pathogenic": 
                cat_labels.append(1)


        cat_refalt.append(refalt)
        cat_blosum.append(blosum)
        cat_seq.append(seq)
        cat_ss.append(ss)

    cat_refalt = np.array(cat_refalt)
    cat_blosum = np.array(cat_blosum)
    cat_seq = np.array(cat_seq)
    cat_ss = np.array(cat_ss)
    cat_labels = np.array(cat_labels)

    print(cat_refalt.shape)
    print(cat_blosum.shape)
    print(cat_seq.shape)
    print(cat_ss.shape)
    print(cat_labels.shape)

    print('Done reading data for category ' + category)

    numPath = (cat_labels == 1).sum()
    numNeut = (cat_labels == 0).sum()
    print('Pathogenic : ' + str(numPath))
    print('Neutral : ' + str(numNeut))

    indata_cat = [cat_refalt, cat_blosum, cat_seq, cat_ss]

    return indata_cat, cat_labels


def load_data_presplit_seqver(h5path, check_flag):

    h5file = open_file(h5path, mode="r", title="All data")

    indata_train, train_labels = load_data_category_seqver(h5file, 'training', check_flag, False)
    indata_val, val_labels = load_data_category_seqver(h5file, 'validation', check_flag, False)
    indata_test, test_labels = load_data_category_seqver(h5file, 'test', check_flag, False)

    return indata_train, indata_val, indata_test, train_labels, val_labels, test_labels
