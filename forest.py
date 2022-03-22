import csv
import pandas as pd
import math

df = pd.read_csv("mutations.csv")
tmp = df.sum(axis=0)
i=0
for columnName in df.columns:
    if i!=0:
        if tmp[i] < 5:                                  #drop mutations with 2 or less samples 
            df.drop(columns=[columnName], inplace=True)
    i=i+1
tmp2 = df.sum(axis=1)
for j in range(0, df.shape[0]):
    if j!=0:
        if tmp2[j] < 1:                                 #drop rows not affected by selected columns
            df.drop([j], inplace=True)
df.reset_index(inplace=True, drop=True)

class Node:                                             #node
    def __init__(self, data):
        self.left = None
        self.right = None
        self.data = data

def printTree(node, level=0):
        if node != None:
            printTree(node.left, level + 1)
            print(' ' * 4 * level + '->', node.data)
            printTree(node.right, level + 1)

def topMutation(df):                  #returns top mutation 
    rows = df.shape[0]
    tmpList=[]
    for(columnName, columnData) in df.iteritems():
        tp = 0                          #for each column initialize its variables to 0
        tn = 0 
        fp = 0
        fn = 0 
        for i in range(rows):
            if df['patients'][i][0] == 'C' and df[columnName][i] == 1:
                tp+=1
            elif df['patients'][i][0] != 'C' and df[columnName][i] == 0:
                tn+=1
            elif df['patients'][i][0] == 'C' and df[columnName][i] == 0:
                fn+=1
            elif df['patients'][i][0] != 'C' and df[columnName][i] == 1:
                fp+=1
        tmp = [columnName, tp-fp]                          #get the tp-fp values for each mutation
        tmpList.append(tmp)     #put all mutations and their tp-fp values into a list that we can sort

    #output list holds all mutations and sorts them by tp-fp
    from operator import itemgetter
    outputlist = sorted(tmpList, key=itemgetter(1), reverse=True)
    return outputlist[0][0]

def createTree(df):
    mutation = topMutation(df)
    left = df.copy()
    right = df.copy()
    for i in range(0, df.shape[0]):
        if df[mutation][i] == 1:
            right.drop([i], inplace=True)
        else:
            left.drop([i], inplace=True)
    left.reset_index(inplace=True, drop=True)
    right.reset_index(inplace=True, drop=True)

    root = Node(mutation)
    root.left = splitNode(left)
    root.right = splitNode(right)
    return root

def splitNode(df):
    cancer = 0
    noncancer = 0
    for i in range(0, df.shape[0]):
        if df['patients'][i][0] == 'C':
            cancer=cancer+1
        else:
            noncancer=noncancer+1

    if cancer/(noncancer+0.0000000000001) > 2.5:           #no further splits when samples have a clear ratio
        return Node('C')
    elif noncancer/(cancer+0.000000000001) > 2.5:
        return Node('NC')
    elif df.shape[0] < 11:
        if cancer > noncancer:
            return Node('C')
        else:
            return Node('NC')
    else:
        mutation = topMutation(df)
        left = df.copy()
        right = df.copy()
        for i in range(0, df.shape[0]):         #split node 
            if df[mutation][i] == 1:
                right.drop([i], inplace=True)
            else:
                left.drop([i], inplace=True)
        left.reset_index(inplace=True, drop=True)
        right.reset_index(inplace=True, drop=True)
        root = Node(mutation)                   #recursively call splitNode 
        root.left = splitNode(left)
        root.right = splitNode(right)
        return root

def predictValue(root, row):       #sample is ran through the tree and determined to be C or NC
    if root.data == 'C':
        return 'C'
    if root.data == 'NC':
        return 'NC'
    if row[root.data] == 1:
        return predictValue(root.left, row)
    else:
        return predictValue(root.right, row)

def getFeatureStats(feature, forests):
    root = 0 
    secondPick = 0
    thirdPick = 0
    for i in range(0, 3):
        for j in range(0, 100):
            if forests[i][j].data == feature:
                root=root+1
            if forests[i][j].right != None:
                if forests[i][j].right.data == feature:
                    secondPick = secondPick+1 
                if forests[i][j].right.right != None:
                    if forests[i][j].right.right.data == feature:
                        thirdPick = thirdPick+1
    return [feature, root, secondPick, thirdPick]
            

def testTree(root, test_data):          #print tree and test it
    printTree(root)
    test_rows = test_data.shape[0]
    accur = []
    sens = []
    spec = []
    prec = []
    miss = []
    fdr = []
    fomr = []
    tn=0
    fp=0
    fn=0
    tp=0
    #find test_data row prediction
    #give tn fp fn tp value
    for i in range(0, test_rows):
        prediction = predictValue(root, test_data.loc[i])
        if prediction == 'C':
            if test_data['patients'][i][0] == 'C':
                tp=tp+1
            else:
                fp=fp+1
        else:
            if test_data['patients'][i][0] != 'C':
                tn=tn+1
            else:
                fn=fn+1

    #prints the tree classifiers
    print('\n')
    print('out of bag # ', test_data.shape[0])
    for i in range(0, test_data.shape[0]):
        print(test_data['patients'][i], ' ', end='')
    print('\n')
    print('tp: ', tp, ' tn: ', tn, ' fp: ', fp, ' fn: ', fn)                #output tests
    #o	Accuracy
    accur.append((tp+tn)/(tp+tn+fp+fn+0.0000000000000001)*100)
    print('Accuracy: ', (tp+tn)/(tp+tn+fp+fn+0.000000000000000000001)*100)
    #o	Sensitivity
    sens.append((tp)/(tp+fn+0.000000000000000000001)*100)
    print('Sensitivity: ', (tp)/(tp+fn+0.000000000000001)*100)
    #o	Specificity
    spec.append((tn)/(tn+fp+0.00000000000000001)*100)
    print('Specificity: ', (tn)/(tn+fp+0.0000000000000000001)*100)
    #o	Precision
    prec.append((tp)/(tp+fp+0.000000000000000001)*100)
    print('Precision: ', (tp)/(tp+fp+0.000000000000000000001)*100)
    #o	Miss rate
    miss.append((fn)/(fn+tp+0.0000000000000001)*100)
    print('Miss rate: ', (fn)/(fn+tp+0.00000000000000001)*100)
    #o	False discovery rate
    fdr.append((fp)/(fp+tp+0.0000000000000000000001)*100)
    print('False Discovery Rate: ', (fp)/(fp+tp+0.000000000000000000000001)*100)
    #o	False omission rate
    fomr.append((fn)/(fn+tn+0.0000000000001)*100)
    print('False Omission Rate: ', (fn)/(fn+tn+0.00000000000000000000001)*100)
    print('\n')
    return [tp, tn, fp, fn, accur, sens, spec, prec, miss, fdr, fomr]

def randomForest(df):               #generate a random forest
    n = math.floor((df.shape[0]*2)/3) 
    forest = []
    for i in range(0, 100):              
        bootstrap_dataset = df.sample(n, frac=None, replace=True, axis=0, ignore_index=True)
        tmp = pd.concat([df, bootstrap_dataset])
        tmp.drop_duplicates(subset=['patients'], keep=False, inplace=True)
        bootstrap_dataset.reset_index(inplace=True, drop=True)
        tmp.reset_index(inplace=True, drop=True)
        tree = createTree(bootstrap_dataset)
        testTree(tree, tmp)
        forest.append(tree)                      #every iteration a tree is made and added to forest
    return forest

def testForest(forest, test_data):
    
    votes = []     #[cancer, noncancer]    
    accur = []
    sens = []
    spec = []
    prec = []
    miss = []
    fdr = []
    fomr = []
    tn=0
    fp=0
    fn=0
    tp=0                            
    for j in range(0, test_data.shape[0]):
        cancer=0
        noncancer=0
        for i in range(0, 100):
            prediction = predictValue(forest[i], test_data.loc[j])
            if prediction == 'C':
                cancer=cancer+1
                if test_data['patients'][j][0] == 'C':
                    tp=tp+1
                else:
                    fp=fp+1
            else:
                noncancer=noncancer+1
                if test_data['patients'][j][0] != 'C':
                    tn=tn+1
                else:
                    fn=fn+1
        if cancer > noncancer:
            votes.append('C')
        else:
            votes.append('NC')
    print('\n')
    print('test data# ', test_data.shape[0])
    for i in range(0, test_data.shape[0]):
        print(test_data['patients'][i], ' Prediction: ', votes[i], end='')
    print('\n')
    print('tp: ', tp, ' tn: ', tn, ' fp: ', fp, ' fn: ', fn)                #output tests
    #o	Accuracy
    accur.append((tp+tn)/(tp+tn+fp+fn+0.0000000000000001)*100)
    print('Accuracy: ', (tp+tn)/(tp+tn+fp+fn+0.000000000000000000001)*100)
    #o	Sensitivity
    sens.append((tp)/(tp+fn+0.000000000000000000001)*100)
    print('Sensitivity: ', (tp)/(tp+fn+0.000000000000001)*100)
    #o	Specificity
    spec.append((tn)/(tn+fp+0.00000000000000001)*100)
    print('Specificity: ', (tn)/(tn+fp+0.0000000000000000001)*100)
    #o	Precision
    prec.append((tp)/(tp+fp+0.000000000000000001)*100)
    print('Precision: ', (tp)/(tp+fp+0.000000000000000000001)*100)
    #o	Miss rate
    miss.append((fn)/(fn+tp+0.0000000000000001)*100)
    print('Miss rate: ', (fn)/(fn+tp+0.00000000000000001)*100)
    #o	False discovery rate
    fdr.append((fp)/(fp+tp+0.0000000000000000000001)*100)
    print('False Discovery Rate: ', (fp)/(fp+tp+0.000000000000000000000001)*100)
    #o	False omission rate
    fomr.append((fn)/(fn+tn+0.00000000000001)*100)
    print('False Omission Rate: ', (fn)/(fn+tn+0.00000000000000000000001)*100)
    print('\n')
    return [tp, tn, fp, fn, accur, sens, spec, prec, miss, fdr, fomr]

def threeFold(df):                    #3 fold cross validation
    n = math.floor(df.shape[0]/3)
    groupOne = df.sample(n, frac=None, replace=False, axis=0, ignore_index=False)
    df = pd.concat([df, groupOne])
    df.drop_duplicates(subset=['patients'], keep=False, inplace=True)
    groupTwo = df.sample(n, frac=None, replace=False, axis=0, ignore_index=False)
    df = pd.concat([df, groupTwo])
    df.drop_duplicates(subset=['patients'], keep=False, inplace=True)
    groupThree = df

    combo1 = pd.concat([groupOne, groupTwo])
    combo2 = pd.concat([groupTwo, groupThree])
    combo3 = pd.concat([groupOne, groupThree])
    groupOne.reset_index(inplace=True, drop=True)
    groupTwo.reset_index(inplace=True, drop=True)
    groupThree.reset_index(inplace=True, drop=True)                     #create three groups 
    combo1.reset_index(inplace=True, drop=True)
    combo2.reset_index(inplace=True, drop=True)
    combo3.reset_index(inplace=True, drop=True)

    stats= []
    forestOne = randomForest(combo1)                        #build a forest with each group 
    stats.append(testForest(forestOne, groupThree))
    forestTwo = randomForest(combo2)
    stats.append(testForest(forestTwo, groupOne))
    forestThree = randomForest(combo3)
    stats.append(testForest(forestThree, groupTwo))
    forests = [forestOne, forestTwo, forestThree]
    statz = []
    for column in df.columns:
        statz.append(getFeatureStats(column, forests))
    from operator import itemgetter
    feature1 = sorted(statz, key=itemgetter(1), reverse=True)
    feature2 = sorted(statz, key=itemgetter(2), reverse=True)
    feature3 = sorted(statz, key=itemgetter(3), reverse=True)

    #print in table form the top 10 mutations by accuracy
    from tabulate import tabulate
    print(tabulate(feature1[0:10]))
    print(tabulate(feature2[0:10]))
    print(tabulate(feature3[0:10]))


    print('Averages')
    tp=(stats[0][0]+stats[1][0]+stats[2][0])/3
    tn=(stats[0][1]+stats[1][1]+stats[2][1])/3
    fp=(stats[0][2]+stats[1][2]+stats[2][2])/3
    fn=(stats[0][3]+stats[1][3]+stats[2][3])/3
    print('tp: ', tp, ' tn: ', tn, ' fp: ', fp, ' fn: ', fn)
    print('Accuracy: ', (tp+tn)/(tp+tn+fp+fn+0.000000000000000000001)*100)
    print('Sensitivity: ', (tp)/(tp+fn+0.000000000000001)*100)
    print('Specificity: ', (tn)/(tn+fp+0.0000000000000000001)*100)
    print('Precision: ', (tp)/(tp+fp+0.000000000000000000001)*100)
    print('Miss rate: ', (fn)/(fn+tp+0.00000000000000001)*100)
    print('False Discovery Rate: ', (fp)/(fp+tp+0.000000000000000000000001)*100)
    print('False Omission Rate: ', (fn)/(fn+tn+0.00000000000000000000001)*100)
    
threeFold(df)


