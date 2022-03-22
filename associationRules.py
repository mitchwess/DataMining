import pandas as pd


df = pd.read_csv("mutations.csv")
F1 = pd.DataFrame(df.sum(axis=0)).T

i=0
for columnName in F1[0:]:
    if i != 0:
        if F1[columnName][0] < 4:                                  #drop mutations with 3 or less samples 
            F1.drop(columns=[columnName], inplace=True)
    i=i+1

F2 = []
i=0
for columnName1 in F1[0:]:
    j=0
    for columnName2 in F1[0:]:                  #for every possible pair
        if columnName1 != columnName2:
            if i!= 0:
                if j!=0:
                    count = 0
                    for sample in range(0, df.shape[0]):           #count how many samples they coexit in 
                        if df[columnName1][sample] == 1 and df[columnName2][sample] == 1:
                            count = count +1
                    if count >= 4:                              #add to pairs if they share 4 or more samples                         
                            F2.append([columnName1, columnName2, F1[columnName1][0], F1[columnName2][0], count])
            j=j+1
    i=i+1

association = []
for i in range(len(F2)):                #for every pair make an association rule for the pair
    totalSamples = df.shape[0]
    support = F2[i][4]/totalSamples
    confidence = F2[i][4]/F2[i][2]
    sXc = support*confidence
    association.append([F2[i][0], F2[i][1], support*100, confidence*100, sXc])  #append support confidence and support*confidence

from operator import itemgetter             #print output
supportResult = sorted(association, key=itemgetter(2), reverse=True)
confidenceResult = sorted(association, key=itemgetter(3), reverse=True)
from tabulate import tabulate 
print('SUPPORT')
print(tabulate(supportResult[0:30]))
print('\n')   
print('CONFIDENCE')
print(tabulate(confidenceResult[0:30]))     

#support threshold above 3.1%
#confidence threshold above 70%
commonRules = []
for i in range(len(association)):           #get rules that pass requirements 
        if association[i][2] > 3.1:
            if association[i][3] > 60:
                commonRules.append(association[i])
from operator import itemgetter
rankings = sorted(commonRules, key=itemgetter(4), reverse=True)
from tabulate import tabulate
print('SELECTED RULES RANKED BY SUPPORT * CONFIDENCE')
print(tabulate(rankings[0:10]))

triples = []
for i in range(len(rankings)):          #joining
    for j in range(len(rankings)):
        if rankings[i][1] != rankings[j][1]: 
            if rankings[i][0] == rankings[j][0]:
                triples.append([rankings[i][0], rankings[i][1], rankings[j][1]])
        if rankings[i][0] != rankings[j][0]:
            if rankings[i][1] == rankings[j][1]:
                triples.append([rankings[i][1], rankings[i][0], rankings[j][0]])
        if rankings[i][0] != rankings[j][1]:
            if rankings[i][1] == rankings[j][0]:
                triples.append([rankings[i][1], rankings[i][0], rankings[j][1]])
        if rankings[i][1] != rankings[j][0]:
            if rankings[i][1] == rankings[j][0]:
                triples.append([rankings[i][1], rankings[i][1], rankings[j][0]])

def testPair(item1, item2, list):
    for i in range(len(list)):
        if list[i][0] == item1:
            if list[i][1] == item2: 
                if list[i][4] < 4:
                    return False
                else:
                    return True
F3 = []
for i in range(len(triples)):
    if testPair(triples[i][0], triples[i][1], F2) == True:     #checks all subsets have 4 or more samples
        if testPair(triples[i][0], triples[i][2], F2) == True:
            if testPair(triples[i][2], triples[i][1], F2) == True:
                F3.append(triples[i])

def findResults(abc):
    onetwo = 0
    count = 0
    all = 0
    for j in range(df.shape[0]):       #finds support and confidence for a given touple 
        if df[abc[0]][j] == 1:
            if df[abc[1]][j] == 1:
                onetwo = onetwo+1
                if df[abc[2]][j] == 1:
                    count=count+1
    support = count/df.shape[0]
    confidence = count/onetwo
    sXc = support*confidence
    return [abc[0], abc[1], abc[2], support*100, confidence*100, sXc]

rules = []
for i in range(len(F3)):
    abc = [F3[i][0], F3[i][1], F3[i][2]]            #each feature gets to be a consequent
    bca = [F3[i][1], F3[i][2], F3[i][0]]
    cab = [F3[i][2], F3[i][0], F3[i][1]]
    rules.append(findResults(abc))
    rules.append(findResults(bca))
    rules.append(findResults(cab))

from operator import itemgetter             #print output
supportResult = sorted(rules, key=itemgetter(3), reverse=True)
confidenceResult = sorted(rules, key=itemgetter(4), reverse=True)
from tabulate import tabulate 
print('SUPPORT')
print(tabulate(supportResult[0:30]))
print('\n')   
print('CONFIDENCE')
print(tabulate(confidenceResult[0:30]))

finalRules = []
for i in range(len(rules)):           #get rules that pass requirements 
        if rules[i][3] > 2.4:
            if rules[i][4] > 60:
                finalRules.append(rules[i])

from operator import itemgetter
tmp = sorted(finalRules, key=itemgetter(5), reverse=True)
from tabulate import tabulate
print('FINAL RULES RANK')
print(tabulate(tmp[0:6]))

