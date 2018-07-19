# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:52:12 2018

@author: Alan
"""
import tellurium as te
import numpy as np
import math
import matplotlib.pyplot as plt
import os

#def prints(what): Used when making prog
#    print(what.keys())
#    for j in what.keys():
#        print(what[j])
        
def RandomGRN(genes,Probabilities, InitVals, Perturbation):
    
    
    Interactions = ["SingleAct", "SingleRep", "and", "or", "Nand", "Nor", "Counter"]
    GRN = {'Rate': [], 'mRNA':[], 'Prot': [],'PDeg':[],'mDeg': [],'Leak': [],'TRa':[], 'TRb':[], 'DualDissoc':[], 'Dissoc':[], 'HillCoeff': []}
    GRNInt = {'TF':[],'TFs':[]}
    StringList=[]
    def ModelInitNames():
        for k in range(genes):
            num = str(k+1)
            GRN['Leak'].append("L" + num)
            GRN['Rate'].append('R' + num +' := ')
            GRN['mDeg'].append('Dm' + num)
            GRN['PDeg'].append('Dp' + num)
            GRN['mRNA'].append('M' + num)
            GRN['Prot'].append('P' + num)
            GRN['TRa'].append('TRa' + num)
            GRN['TRb'].append('TRb' + num)
            GRN['HillCoeff'].append('H' + num)
            GRN['Dissoc'].append('K' + num)
            GRN['DualDissoc'].append('Kd' + num)
        return
    def AssignProteins():
        nums =[int(math.ceil(float(genes)/2)), int(math.floor(float(genes)/2))]
        RanAct = np.random.choice(GRN['Prot'], nums[np.random.randint(2)])
        GRNInt.update({'activators': RanAct})
        index1 = np.arange(len(GRNInt['activators']))
        index2 = np.sort(-index1)
        for m in range(len(index1)):
            ActivatorCur =GRNInt['activators'][m]
            Activators1 = GRNInt['activators'][index2[m]:]
            Activators2 =  GRNInt['activators'][:index1[m]]
            if m > len(index1)-2:
                Activators1 = []
            if ActivatorCur in Activators1 or ActivatorCur in Activators2:
                AssignProteins()
        return
    def DetermineInteractors():
        Proteins = []
        ProbabilityMatrix = []
        ProbSingle = Probabilities[0]
        ProbDouble = Probabilities[1]
        ProbCounter = Probabilities[2]
        for h in range(len(GRN['Prot'])):
            Proteins.append(GRN['Prot'][h])
        for s in GRNInt['activators']:#Issue generates too many interactions when activators has multiple regulatory interactions per gene
            if s in Proteins:
                Proteins.remove(s)
            GRNInt.update({'repressors': Proteins})
        Psing = ProbSingle/2
        Pdub = ProbDouble/4
        ProbabilityMatrix = [Psing, Psing, Pdub, Pdub, Pdub, Pdub, ProbCounter]
        GRNInt.update({'interactions': np.random.choice(Interactions,genes,True,ProbabilityMatrix)})
        for Inter in range(len(GRNInt['interactions'])):  
            reg = GRNInt['interactions'][Inter]    
            RegType = 'repressors'
            if 'Single' in reg :
                if 'Act' in reg:
                    RegType = 'activators'
                GRNInt['TF'].append([reg, [GRNInt[RegType][np.random.randint(len(GRNInt[RegType]))]]])
            else:
                Tps={'activators':[], 'repressors':[]}
                for m in Tps.keys():
                    n=1
                    Tps[m].append(np.random.choice(GRNInt['activators'],2))
                    while Tps[m][0][n-1] == Tps[m][0][n]:
                        Tps[m][0][n] =np.random.choice(GRNInt['activators'])
                TypeA = 'activators'
                TypeR = 'repressors'
                if reg == "and" or reg == "or":
                    GRNInt['TF'].append([reg, [Tps[TypeA][0][0],Tps[TypeA][0][1]]])
                elif reg == "Nand" or reg == "Nor":
                    GRNInt['TF'].append([reg, [Tps[TypeR][0][0],Tps[TypeR][0][1]]])
                else:
                    GRNInt['TF'].append([reg, [Tps[TypeA][0][0],Tps[TypeR][0][1]]] )
        return 
    def ValueGeneration():
        Species=''
        Values = {'M':[InitVals[0]], 'P':[InitVals[1]],'Dp':[InitVals[2]],'Dm':[InitVals[3]],'L':[InitVals[4]],'TRa':[InitVals[5]], 'TRb':[InitVals[6]],'Kd':[InitVals[7]], 'K':[InitVals[8]], 'H': [InitVals[9]]}
        Keys = Values.keys()
                    
            
        for i in Keys:
            for n in np.arange(1, genes+1):
                val = 0
                if i == 'M' or i=='4P':
                    Values[i].append(val)
                    Species += i + str(n) + ' = ' + str(val) + '; \n'
                else:
                    while val <= 0:
                        val = round(np.random.normal(Values[i][0],.25),3)
#                        val = np.random.normal(Values[i][0],.25)
#                        if i == 'M' or i=='P':
#                            val = round(val,0)
#                        else:
#                            val = round(val,3)
                Values[i].append(val)
                Species += i + str(n) + ' = ' + str(val) + '; \n'
        StringList.append(Species)
        return
    #%% Identify proteins regulating genes
    def ChooseProtein():
        for i in range(genes):
            TF=[]
            string = GRNInt['TF'][i][1]
            for n in string:
                TF.append(int(n.replace("P","")))
            GRNInt['TFs'].append(TF)
    #print(GRNInt['TFs'])
        return
    
    
    #%%Generate Rxn Rates
    def GetReactionRates():
        ReactionRates=[]
        WorkingString=''
        for i in range(genes):
            Reg = GRNInt['TF'][i][0]
            
            WorkingString += 'Rm_' + str(i+1) + ' := '
            TF1 = GRNInt['TFs'][i][0] - 1
            P1 = GRN['Prot'][TF1]
            K1 =  GRN['Dissoc'][TF1]
            HCoeff = GRN['HillCoeff'][i]
            num = ''
            denom = ''
            frac = ''
            if 'Single' in Reg:
                frac = '(' + P1 + '*' + K1 + ')'
                if 'Act' in Reg:
                    Hill = '(' + frac + '/(1 + ' + frac + '))'
                else:
                    Hill = '(1 /(1 + ' + frac + '))'
                WorkingString = GRN['Leak'][i] + '+' + GRN['TRa'][i] + '*' + Hill + '-' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i] + '\n'
                
            elif len(GRNInt['TFs'][i]) >1 :
                TF2 = GRNInt['TFs'][i][1] - 1
                
                P2 = GRN['Prot'][TF2]
                K2 =  GRN['Dissoc'][TF2]
                K3 = GRN['DualDissoc'][TF2]
                if 'Counter' in Reg:
                    num = K1 + '*' + P1
                    denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + '+' + K3 + '*' + P1 + '*' + P2 + ')'#K3 is a place holder.. need to generate dissoc associated with dual reg 
                    eq = '(' + num  + '/' + denom + ')'
                            
                elif  not('N' in Reg) and("and" in Reg or 'or' in Reg):
                    if 'and' in Reg:
                        num = '(' + K1 + '*' + K2 + '*' + P1 + '*' + P2  + ')'
                        denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + '+' + K1 + '*' + K2 + '*' + P1 + '*' + P2 + ')'
                        eq = '(' + num  + '/' + denom + ')'
                    else:
                        num = '(' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + ')'
                        denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2  + ')'
                        eq = '(' + num  + '/' + denom + ')'
                else:
                    denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + '+' + K3+ '*' + P1 + '*' + P2 + ')'
                    if 'and' in Reg:
                        num = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + ')'
                        eq = '(' + num  + '/' + denom + ')'
                    else:
                        eq =  '(1' + '/' +  denom + ')'
                WorkingString = GRN['Leak'][i] + '+' + GRN['TRa'][i] + '*' + eq + '-' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i]+ '\n' 
        
            ReactionRates.append(WorkingString)
        return ReactionRates
    #%%Events-- attempt to perturb a user defined genes trancription at the associated time
    def GetEvents():
        EventString = ''
#        for m in range(genes):
#            EventString += 'Em' + str(m+1) + ': at(' + GRN['mRNA'][m] + '<0): ' + GRN['mRNA'][m] + '=0;\n'
#            EventString += 'Ep' + str(m+1) + ': at(' + GRN['Prot'][m] + '<0): ' + GRN['Prot'][m] + '=0;\n'    Old Catch. Previous issue.
        counter = 0
        for k in Perturbation.keys():
            counter +=1
            EventString += 'Ed' + str(counter) + ': at(time >'  + str(Perturbation[k][0]) +'): ' + k + '= ' + str(39) +';\n' 
        StringList.append(EventString)
        return
    #%% Reaction Model Generation
    def GetReactions():
        ReactionM =''
        ReactionP =''
        ReactionSumString= 'model Random_GRN()\n' 
        
        for i in np.arange(1, genes+1):
            M = GRN['mRNA'][i-1]
            P = GRN['Prot'][i-1]
            ReactionM = 'Rm' + str(i) + ':' + '$N' + '=> ' + M + ';' + Rates[i-1] + '\n'
            ReactionP = 'Rp' + str(i) + ':' + '$A' + '=> ' + P + ';' + GRN['TRb'][i-1] + '*' + M + '-' + P + '*' + GRN['PDeg'][i-1] + '\n'
            ReactionSumString += ReactionM +ReactionP
        StringList.append(ReactionSumString)
          
        
    
    def GetModelString():
        antStr = ''
        for i in StringList:
            antStr += i
        antStr += 'end'
        return antStr 
    #%%
    ModelInitNames()
    AssignProteins()
    DetermineInteractors()
    ChooseProtein()
    Rates = GetReactionRates()
    GetReactions()
    ValueGeneration()
#    GetEvents()
    antStr = GetModelString()
    return antStr
#RandomGRN(10,[0.4,0.4,0.2], [0,0,0.5,0.75,0.01,0.5,15,0.65,0.65,1])

#%%
def GetModel(genes,Probabilities, InitVals, Name):
    tries = 0
    while tries < 10:
       try:
          antStr = RandomGRN(genes,Probabilities, InitVals,{'Tra1':[10]})
          r = te.loadAntimonyModel(antStr)
          print r.steadyState()
          break
       except:
           tries = tries + 1
           print "Failed to find steady state"
       
#    r.draw(layout='fdp')  
#    print(antStr)
    return r

# fileName == entire path
def writecsvFile (fileName, r, data):
    
    names = r.getFloatingSpeciesIds()
    header_string = names[0] 
    for name in names[1-len(names):]:
        header_string = header_string + ',' + name 
 
    fh = open(fileName, "wb")
    fh.write(header_string)
    fh.write('\n')
    for x in xrange(data.shape[0]):
        for y in xrange(data.shape[1]):
            val = data[x][y]
            s1 = "{:6.5f}".format(val)
            fh.write (s1)
            if y != data.shape[1] - 1:
                fh.write (',')
        fh.write('\n')
    fh.close()
def Output(data, Name,r,rSeed,result,NoisyResult):
    
    
    for i in range(len(data)):
        if i == 1 and data[i] == 'y':
            writecsvFile(fileName + Name + '_Results.csv',r,result)
        elif data[i] == 'y':
            writecsvFile(fileName + Name + '_Noisy_result.csv',r,NoisyResult)
    if rSeed != 0:
        fh = open(fileName + Name + '_Seed.txt', 'wb')
        fh.write('Random Seed = ' + str(rSeed))
    sbmlStr = r.getSBML()
    te.saveToFile (fileName + Name + '_Model.xml', sbmlStr)
#%%
#antStr = RandomGRN(10,[0.4,0.4,0.2], [0,0,0.5,0.75,0.01,0.5,15,0.65,0.45,1])
#r.draw(layout='fdp')
#Values = {'M':[], 'P':[],'Dp':[],'Dm':[],'L':[],'TRa':[], 'TRb':[],'Kd':[], 'K':[], 'H': []}
def RunModel(genes,tmax,InitProb,Percent,Name,DataOut,rSeed, Steps):
    if rSeed != 0:
        np.random.seed (rSeed)
    #    Param = []
#    t=0.19142857142857142
    Percent /= 100
#    for p in np.linspace(.05,.85,4):
    InitParams = [0,0,0.5,0.9,0.8,30,30,0.2,0.5,1]
    r = GetModel(genes,InitProb, InitParams, Name) 
    tries = 0
#    Param = [[r.getGlobalParameterIds()]]
    while tries < 10:
        try:
            r.reset()
            result = r.simulate(0,tmax, Steps)
            NoisyResult = np.zeros([len(result[:,0]), len(result[0,:])])
            for k in range(len(result[:,0])):
                for i in range(len(result[0])):
                    CurrVal = -1
                    while CurrVal < 0:
                        CurrVal = result[k,i] + np.random.normal(0,Percent*result[k,i])
                    NoisyResult[k,i] = CurrVal
            vars = {'P':[],'M':[]}
            for e in vars.keys():
                plt.figure()    
                for k in np.arange(1,genes+1):
                    vars[e].append(e + str(k))
                    for i in range(len(vars[e])):
                        if e=='P':
                            Max = np.max(result[:,1:genes+1])
                            Min = np.min(result[:,1:genes+1])
                            plt.grid(color='k', linestyle='-', linewidth=.4)
                            plt.ylim(Min,Max*1.1)
                            plt.xlim(0,tmax)
                            plt.yticks(np.arange(Min,Max*1.1,Max/12))
                            plt.plot (result[:,0],result[:,i+1], label = vars[e][i])
                            plt.legend(vars[e])
                        else:
                           Max = np.max(result[:,genes+1:])
                           Min = np.min(result[:,genes+1:])
                           plt.grid(color='k', linestyle='-', linewidth=.4)
                           plt.ylim(Min,Max*1.1)
                           plt.xlim(0,tmax)
                           plt.yticks(np.arange(Min,Max*1.1,Max/12))
                           plt.plot (result[:,0],result[:,i+genes+1], label = vars[e][i])
                           plt.legend(vars[e])             
                         
            vars = {'P':[],'M':[]}             
            for e in vars.keys():
                plt.figure()    
                for k in np.arange(1,genes+1):
                    vars[e].append(e + str(k))
                    for i in range(len(vars[e])):
                        if e=='P':
                            Max = np.max(NoisyResult[:,1:genes+1])
                            Min = np.min(NoisyResult[:,1:genes+1])
                            plt.grid(color='k', linestyle='-', linewidth=.4)
                            plt.ylim(Min,Max*1.1)
                            plt.xlim(0,tmax)
                            plt.yticks(np.arange(Min,Max*1.1,Max/12))
                            plt.plot (NoisyResult[:,0],NoisyResult[:,i+1], label = vars[e][i])
                            plt.legend(vars[e])
                        else:
                            Max = np.max(NoisyResult[:,genes+1:])
                            Min = np.min(NoisyResult[:,genes+1:])
                            plt.grid(color='k', linestyle='-', linewidth=.4)
                            plt.ylim(Min,Max*1.1)
                            plt.xlim(0,tmax)
                            plt.yticks(np.arange(Min,Max*1.1,Max/12))
                            plt.plot (NoisyResult[:,0],NoisyResult[:,i+genes+1], label = vars[e][i])
                            plt.legend(vars[e])
             
            #r.draw(layout='fdp')      
            
#            Param.append(r.getGlobalParameterValues())
            break
        except:
             tries = tries + 1
             print "Failed to solve Model"
    
    Output(DataOut, Name, r, rSeed, result, NoisyResult)
    
    #if tries == 10:
        #     tries = 0
        #     r = GetModel(genes,InitProb, InitParams, Name)
while True:
    try:
        # Note: Python 2.x users should use raw_input, the equivalent of 3.x's input
#        break
        Name = str(raw_input("Please enter the name of your network(str): "))
        fileName = 'c:\\' + Name  +'\\'
        os.mkdir(fileName)
        genes = int(raw_input("Please enter the number of genes in your network(int): "))
        InitProb=[]
        InitProb.append(float(raw_input("Please enter the probability of Single regulation.(float): ")))
        InitProb.append(float(raw_input("Please enter the probability of Double regulation.(float): ")))
        CountProb = round(1-np.sum(InitProb),3)
        print('The probabilty of Counter regulation will be ' + str(CountProb) + '\n' )
        InitProb.append(CountProb)
        
        while any(i<0 for i in InitProb) == True or np.sum(InitProb) !=1.0:
             print("Sorry, the sum of the probabilities must equal 1 or at least one is negative.")
             InitProb[0] = float(raw_input("Please enter the probability of Single regulation.(float): "))
             InitProb[1] = float(raw_input("Please enter the probability of Double regulation.(float): "))
             CountProb = round(1-np.sum(InitProb[0:2]),3)
             InitProb[2] = CountProb
             print('Your new probabilty of Counter regulation will be ' + str(CountProb) + '\n' )
#             break
        tMax = int(raw_input("Please enter the duration of your network simulation(int): "))
        Steps = int(raw_input("Please enter number of steps for your network simulation(int): "))
        NLevel = float(raw_input("Please enter the percentage of the noise level(float): "))
        DataOut = []
        DataOut.append(str(raw_input("Would you like simulated data exported(y/n): ")))
        DataOut.append(str(raw_input("Would you like noisy simulated data exported(y/n): ")))
        rSeed = int (raw_input ('Please enter an optional random seed (0 for no seed): '))
        #RunModel(genes,tMax,InitProb,NLevel,Name,DataOut,rSeed)
    except WindowsError:
        print("The directory exists, choose another name.")
    except ValueError:
        print("There was an issue with your last entry.")
        #better try again... Return to the start of the loop
        continue
    else:
        #age was successfully parsed!
        #we're ready to exit the loop.
        break
RunModel(genes,tMax,InitProb,NLevel,Name,DataOut,rSeed, Steps)
##%%     
#print r.getRatesOfChange()
#
#print r.model.getFloatingSpeciesInitConcentrations()
#
#print r.getGlobalParameterValues()
#
#print r.getGlobalParameterIds()
#
#print r.getFloatingSpeciesIds()
#
#print r.getFloatingSpeciesConcentrations()