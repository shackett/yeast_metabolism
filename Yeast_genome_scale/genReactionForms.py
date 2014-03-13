# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 11:09:06 2013

This will create reaction formulas identical to the one described in the
Liebermeister paper. The input is the conversion table used for Juns C script

@author: vitoz@ethz.ch
"""

import sys
import re
# Definitions

def nr2id(inp):
    inp = str(inp)
    while len(inp) < 4:
        inp = '0'+inp
        
    return(inp)

def tab2reactions(convtab):
    """ Parses the table in a transparent structure that can easely translated into reaction
forms:
rxnList[rxnID](
['rct']: contains reactands fields [subsid] with subfields ['stoi'],['actsit']
['act']: contains modulator fields [modid] with subfields
['actsit'](0 for noncomp, n for comp),
['hill'](hill coefficient, n > 0: fixed, n == 0: to be determined)
['inh']: contains modulator fields [modid] with subfields like act+
['subtype'] for inh(competitive, noncompetitive or uncompetitive) and for act(MM or CC)
MM = (1 / ((Kd/L)^n + 1), CC = (1 + (L/D))^n
['rev']: 0 or 1)
"""
    
    rxnList = dict()
    
    for line in convtab:
        
        line = str(line)
        line = line.replace('\n','')
        line = line.replace('"','')
        vec = line.split('\t') # 0 ReactionID (r_xxxx), 1 Reversibility (bool), 2 type(rct, act, inh),
        #3 substrateId (t_xxxx) , 4 bindingSite (double),5 stoi(int),6 hill (int),
        
        rxnID = vec[0].replace('r_','')
        
        if not(rxnID in rxnList.keys()):
            rxnList[rxnID] = dict()
            rxnList[rxnID]['rev'] = vec[1]
            rxnList[rxnID]['rct'] = dict()
            rxnList[rxnID]['act'] = dict()
            rxnList[rxnID]['inh'] = dict()
        
        # go through the substrate, if not yet present in the dictionary add them
        # each substrate is a dictionary with two values: active site and stoichiometry
        
        metID = nr2id(vec[3].replace('t_',''))
        if vec[2] == 'rct':
            if not(metID in rxnList[rxnID].keys()):
                rxnList[rxnID]['rct'][metID] = dict()
                rxnList[rxnID]['rct'][metID]['actsit'] = float(vec[4])
                rxnList[rxnID]['rct'][metID]['stoi'] = int(vec[5])
            else:
                raise NameError(metID + ' twice in table for reaction ' + rxnID)
        
        # activator and inhibitor dictionaries are added to their respective
        # rection-specific dictionaries
               
        # activators
        elif vec[2] == 'act':
            rxnList[rxnID]['act'][metID] = dict()
            rxnList[rxnID]['act'][metID]['hill'] = float(vec[6])
            rxnList[rxnID]['act'][metID]['actsit'] = float(vec[4])
            rxnList[rxnID]['act'][metID]['subtype'] = vec[7].lower()
        
        # inhibitors 
        elif vec[2][0:3] == 'inh':
            rxnList[rxnID]['inh'][metID] = dict()
            rxnList[rxnID]['inh'][metID]['hill'] = float(vec[6])
            rxnList[rxnID]['inh'][metID]['actsit'] = float(vec[4])
            rxnList[rxnID]['inh'][metID]['subtype'] = vec[7].lower()
            
        else:
            print(vec[2] + ' not a valid value for type.')
        
    return(rxnList)


# Write rate laws:

def rateLaws(rxnList,formModus):
    """uses as an input the rxList
formModus:
- rm: Reversible Menten with substrate inhibition (=Juns reactionforms)
- cc: Convinience Kinetics (does only support uncompetitive inhibition)
"""
      
    eqList = dict()
    texeqList = dict()
    
    for rxn in rxnList.keys():
        subst_num = str()
        prod_num = str()
        subst_km = str()
        prod_km =str() # only needed if irreversible backwards
        
        # create first the numerator, as it is invariant over all rate laws
        for rct in rxnList[rxn]['rct'].keys():
            tID = 't_' + rct
            km = 'K_r_' + rxn + '_' + tID
            #numerator
            t_num = tID+'^'+str(abs(rxnList[rxn]['rct'][rct]['stoi']))
            # process substrates
            if rxnList[rxn]['rct'][rct]['stoi'] < 0:
                if (subst_num == str()):
                    subst_num = t_num
                    subst_km = km + '^' + str(abs(rxnList[rxn]['rct'][rct]['stoi']))
                else:
                    subst_num = subst_num + '*' + t_num
                    subst_km = subst_km + '*' + km + '^' + str(abs(rxnList[rxn]['rct'][rct]['stoi']))
            else:
                if prod_num == str():
                    prod_num = t_num
                    prod_km = km + '^' + str(abs(rxnList[rxn]['rct'][rct]['stoi']))
                else:
                    prod_num = prod_num + '*' + t_num
                    prod_km = prod_km + '*' + km + '^' + str(abs(rxnList[rxn]['rct'][rct]['stoi']))
        
        # Create denominator in the different reaction form cases
        if formModus == 'rm':
            subSites = dict()
            prodSites = dict()
            denSites = dict()
            for rct in rxnList[rxn]['rct'].keys():
                tDen = '(t_' + rct + '/' +'K_r_' + rxn + '_t_' + rct + ')'
                for i in range(abs(rxnList[rxn]['rct'][rct]['stoi'])):
                    if rxnList[rxn]['rct'][rct]['stoi'] < 0:
                        subSites[str(rxnList[rxn]['rct'][rct]['actsit']+i*0.1)]=tDen
                    else:
                        prodSites[str(rxnList[rxn]['rct'][rct]['actsit']+i*0.1)]=tDen
            # sort all sites in a decreasing order
            for site in sorted(set(prodSites.keys()+subSites.keys()),reverse=1):
                if site in denSites.keys():
                    1
                elif (site in prodSites.keys()) and (site in subSites.keys()):
                    denSites[site] = '(1+' + subSites[site] + '+' + prodSites[site]
                elif (site in prodSites.keys()) and float(site).is_integer():
                    denSites[site] = '(1+' + prodSites[site]
                elif (site in subSites.keys()) and float(site).is_integer():
                    denSites[site] = '(1+' + subSites[site]
                elif site in subSites.keys():
                    intSite = str(float(site).__floordiv__(1))
                    if intSite in prodSites.keys():
                        denSites[intSite] = '(' + prodSites[intSite] + '+' + '(1+' + subSites[site] + ')*(1+' + subSites[intSite] +')'
                    else:
                        denSites[site] = '(1+' + subSites[site]
                elif site in prodSites.keys():
                    intSite = str(float(site).__floordiv__(1))
                    if intSite in subSites.keys():
                        denSites[intSite] = '(' + subSites[intSite]+ '+' + '(1+' + prodSites[site] + ')*(1+'+ prodSites[intSite] + ')'
                    else:
                        denSites[site] = '(1+' + prodSites[site]
                        
            # Add the inhibitors to the denominators
            for rct in rxnList[rxn]['inh'].keys():
                tDen = '(t_' + rct + '/' +'K_r_' + rxn + '_t_' + rct + ')'
                if rxnList[rxn]['inh'][rct]['hill'] == 0:
                    tDen = tDen + '^' + 'KH_r_' + rxn + '_t_' + rct
                elif rxnList[rxn]['inh'][rct]['hill'] != 1:
                    tDen = tDen + '^' + str(rxnList[rxn]['inh'][rct]['hill'])
                
                  
                site = str(rxnList[rxn]['inh'][rct]['actsit'])
                if not(site in denSites.keys()):
                    site = str(float(site).__floordiv__(1))
                if rxnList[rxn]['inh'][rct]['subtype'] == 'uncompetitive': # uncompetitive inh
                    denSites['i'+rct] = '(1+' + tDen
                elif rxnList[rxn]['inh'][rct]['subtype'] == 'competitive':# comp inhibition
                    denSites[site] = denSites[site] + '+' + tDen
                elif rxnList[rxn]['inh'][rct]['subtype'] == 'noncompetitive':# noncomp inhibition
                    denSites[site] = denSites[site] + '+' + tDen
                    denSites['i'+rct] = '(1+' + tDen
                else: raise NameError('invalid inhibitor type (valid: uncompetitive,competitive,noncompetitive) in ' + rxn)
            
            # stitch the denominator together
            denTerm = ''
            for site in denSites.keys():
                denTerm = denTerm + denSites[site] + ')*'
            denTerm = denTerm[:-1]
                
        # convinience kinetics denominator
        elif formModus == 'cc':
            subSites = dict()
            prodSites = dict()
            for rct in rxnList[rxn]['rct'].keys():
                tDen = '(t_' + rct + '/' +'K_r_' + rxn + '_t_' + rct + ')'
                for i in range(abs(rxnList[rxn]['rct'][rct]['stoi'])):
                    if i == 0:
                        t_den = '(1'
                    t_den = t_den + '+' + tDen + '^' + str(i+1)
                if rxnList[rxn]['rct'][rct]['stoi'] < 0:
                    subSites[rct] = t_den
                else:
                    prodSites[rct] = t_den
                    
            denTerm = ''
            for site in subSites.keys():
                denTerm = denTerm + subSites[site] + ')*'
            
            denTerm = denTerm[:-1] + '+'
            for site in prodSites.keys():
                denTerm = denTerm + prodSites[site] + ')*'
            denTerm = denTerm[:-1]
            # Add the inhibitors
            for rct in rxnList[rxn]['inh']:
                tDen = '(t_' + rct + '/' +'K_r_' + rxn + '_t_' + rct + ')'
                if rxnList[rxn]['inh'][rct]['hill'] == 0:
                    tDen = tDen + '^' + 'KH_r_' + rxn + '_t_' + rct
                elif rxnList[rxn]['inh'][rct]['hill'] != 1:
                    tDen = tDen + '^' + str(rxnList[rxn]['inh'][rct]['hill'])
                denTerm = '(' + denTerm + ')*(1+' + tDen + ')'
            
            
        # Stitch the equations together:
        # reversible:
        if rxnList[rxn]['rev'] == '0':
            eq = '1/('+subst_km +')*(' + subst_num + '-' + prod_num + '/Keqr_' + rxn + ')/(' + denTerm + ')'
            texeq = '\\frac{\\frac{1}{' + subst_km + '}*\\left('+ subst_num + '-' + '\\frac{' + prod_num + '}{Keqr_' + rxn + '}\\right)}{' + denTerm + '}'
            
        elif rxnList[rxn]['rev'] == '1':#irreversible
            eq = '(('+ subst_num +')/('+ subst_km +'))/(' + denTerm + ')'
            texeq = '\\frac{\\frac{'+ subst_num +'}{'+ subst_km +'}}{' + denTerm + '}'
        else: #irreversibly backwards
            eq = '-(('+ prod_num +')/('+ prod_km +'))/(' + denTerm + ')'
            texeq = '-\\frac{\\frac{'+ prod_num+'}{'+ prod_km +'}}{' + denTerm + '}'
        # Add the activator term:
        for rct in rxnList[rxn]['act']:
            tDen = '(t_' + rct + '/' +'K_r_' + rxn + '_t_' + rct + ')'
            if rxnList[rxn]['act'][rct]['hill'] == 0:
                tDen = tDen + '^' + 'KH_r_' + rxn + '_t_' + rct
            elif rxnList[rxn]['act'][rct]['hill'] != 1:
                tDen = tDen + '^' + str(rxnList[rxn]['act'][rct]['hill'])
            eq = '(1+' + tDen + ')*' + eq
            texeq = '(1+' + tDen + ')*' + texeq
        eq = re.sub('\^1(?![0-9\.])','',eq)
        eq = re.sub('\\)\\(',')*(',eq)
        texeq = re.sub('\^1(?![0-9\.])','',texeq)
        texeq = re.sub('\\)\\(',')*(',texeq)
        
        # do some latex editing
        texeq = re.sub('\((?=t\_.{4}\/K.{14}\))',r'\\frac{',texeq)
        texeq = re.sub('\/(?=K.{14}\))','}{',texeq)
        texeq = re.sub('(?<=}{K.{14})\)','}',texeq)
        texeq = re.sub('(?<!left)\(',r'\\left(',texeq)
        texeq = re.sub('(?<!right)\)',r'\\right)',texeq)
        texeq = texeq.replace('_','\_')
        texeq = re.sub(r'\\',r'\\\\',texeq)
        eqList[rxn] = eq
        texeqList[rxn] = texeq
    return(dict(eqList = eqList,texeqList = texeqList))

def EqList2R(eqLists):
    """ Loop through the file, where each line corresponds to one reaction.
Remove brackets and reformat it to an R file that stores all reactions
in a variable rxnf, which is a list with the reactions, named by their
reaction ID"""
    eqList = eqLists['eqList']
    texeqList = eqLists['texeqList']
    for rxn in eqList.keys():
        line = "rxnf[[\'r_"+ rxn + '-'+ str(mode) + "\']]$rxnForm <- as.formula( ~I("+eqList[rxn]+')+0)\n'
        sys.stdout.write(line)
        line = "rxnf[[\'r_"+ rxn + '-'+ str(mode) + "\']]$rxnFormTex <- '"+texeqList[rxn]+"'\n"
        sys.stdout.write(line)



# Script start

#outname = '/home/vitoz/Dropbox/Rabino/git/vzcode/liebermeisterRxn'

""" This script should get the input from the standardinput, d.h. it can be called from the command line with
it should get its input with the following columns
ReactionID (r_xxxx), reversibility (bool, 0= reversible)type(rct, act, inh), substrateId (t_xxxx) , bindingSite (double), stoi(int), hill (int),
]"""


#inpTable = sys.stdin.readlines()
#inpTable = open('./pythonTest.py','r').readlines
inpTable = open('./test3.py','r').readlines()
mode = inpTable[0].replace('"','').replace('\n','') # cc or rm
inpTable = inpTable[1:]
rxnList = tab2reactions(inpTable)
eqLists = rateLaws(rxnList,mode)
EqList2R(eqLists)

