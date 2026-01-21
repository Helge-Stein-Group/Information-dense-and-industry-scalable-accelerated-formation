
import os
import pandas as pd
import numpy as np
import h5py
import hdfdict
import matplotlib.pyplot as plt
import re
from datetime import datetime

def getCyclelimitListEOLPuls(cycleindex:list, stepindex:list,current:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        cycleindexlist = []
        stepindexlist = []
        currentlist = []

        for l in range(len(current)-1):
            cycleindexlist.append(cycleindex[l])
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l]!=stepindex[l-1] and l >0 and stepindex[l+1]==stepindex[l]and stepindex[l]!=2 and stepindex[l]!=1 ):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist

def getCyclelimitListEOLPulsLongterm(cycleindex:list, stepindex:list,current:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        cycleindexlist = []
        stepindexlist = []
        currentlist = []

        for l in range(len(current)-1):
            cycleindexlist.append(cycleindex[l])
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l]!=stepindex[l-1] and l >0 and stepindex[l+1]==stepindex[l]and stepindex[l]!=2 and stepindex[l]!=1 and stepindex[l]!=3 and stepindex[l]!=5 and stepindex[l]!=7 and stepindex[l]!=8 and stepindex[l]!=10):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist

def getCyclelimitListEOLCCCV(cycleindex:list, stepindex:list,current:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        cycleindexlist = []
        stepindexlist = []
        currentlist = []
        for l in range(len(current)-1):
            cycleindexlist.append(cycleindex[l])
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l]!=stepindex[l-1] and l >0 and stepindex[l+1]==stepindex[l] and stepindex[l]!=1 and stepindex[l]!=2 ):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist

def getCyclelimitListEOLCCCVLongterm(cycleindex:list, stepindex:list,current:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        cycleindexlist = []
        stepindexlist = []
        currentlist = []
        for l in range(len(current)-1):
            cycleindexlist.append(cycleindex[l])
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l]!=stepindex[l-1] and l >0 and stepindex[l+1]==stepindex[l]and stepindex[l]!=2 and stepindex[l]!=1 and stepindex[l]!=3 and stepindex[l]!=5 and stepindex[l]!=7 and stepindex[l]!=8 and stepindex[l]!=10):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist

def getCyclelimitListEOLCCCV4xFormation(cycleindex:list, stepindex:list,current:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        cycleindexlist = []
        stepindexlist = []
        currentlist = []

        for l in range(len(current)-1):
            cycleindexlist.append(cycleindex[l])
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l]!=stepindex[l-1] and l >0 and stepindex[l+1]==stepindex[l]):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist


def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval
        
def natural_keys(text):
     return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

def getChannelData():
        '''
        Reading in csv file of the testname, poping some parameters to avaoid discspace problems and saving remaining as a dictionary.

        Passed Parameters
        ------------------
        datapath(string): Path to data choosen before when using getpossibletestnames()
        ------------------

        Returns
        ------------------
        datadict(dictionary): Dictionary containing all read-in data from the csv outputfile from arbin.
        ------------------
        '''
        content = os.listdir()
        print(str(content))
        df = pd.read_csv(str(content[0]), delimiter=',')
        df.pop('Date_Time')
        df.pop('Step_Time(s)')
        #print(df)   
        df = df.astype(np.float64)
        datadict = df.to_dict()
        if len(content) > 1:
                print('There are more csv files in.')
                content.sort(key=natural_keys)
                for i in range(1,len(content)):
                    df2 = pd.read_csv(str(content[i]), delimiter=',')  
                    df2.pop('Date_Time')
                    df2.pop('Step_Time(s)')
                    df = pd.concat([df, df2])
                df.set_index('Data_Point', inplace = True)
                completeDict = df.to_dict()
                return completeDict
        else:
            return datadict

def serperaterawData(data:dict):
        '''
        Serperates and saves the raw data saved in a dictionary to lists of each in future used parameter.

        Passed Parameters
        ------------------
        data(Dictioanry): Dictionary containing all information from the read in csv result file from arbin.
        ------------------

        Returns
        ------------------
        testTime(list of int): List containing increasing testtime(int) read-in from csv-datafile.
        cycleIndex(list of int): List containing the corresponging cycleindex to each datapoint read-in from csv-file.
        I(list of float): List containing all Currentvalues(float) read-in from csv-datafile.
        V(list of float): List containing all Voltagevalues(float) read-in from csv-datafile.
        chargeEnergy(list of float): List containing all Chargeenergyvalues(loat) read-in from csv-datafile.
        dischargeEnergy(list of float): List containing all Dischargeenergyvalues(float) read-in from csv-datafile.
        ------------------
        '''
        testTime=list(data['Test_Time(s)'].values())
        #datapoint = list(data['Data_Point'].values())                
        cycleIndex =list(data['Cycle_Index'].values())
        stepIndex = list(data['Step_Index'].values())
        I = list(data['Current(A)'].values())
        V = list(data['Voltage(V)'].values())
        chargeEnergy = list(data['Charge_Energy(Wh)'].values())
        dischargeEnergy = list(data['Discharge_Energy(Wh)'].values())
        #chargeCapacity = list(data['Charge_Capacity(Ah)'].values())
        #dischargeCapacity = list(data['Discharge_Capacity(Ah)'].values())
        return  testTime,cycleIndex,stepIndex,I,V,chargeEnergy,dischargeEnergy
        
def getandsafeCycledata(cycleList:list,testTime:list,current:list, voltage: list, chargeEnergy:list, dischargeEnergy:list,analysisDict:dict, cyclenamelist:list):
        '''
        Slices and saves data for each cycle.
        
        
        Passed Parameters
        ------------------
        cycleList(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        testTime(list of int): List containing increasing testtime(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(float) read-in from csv-datafile.
        voltage(list of float): List containing all Voltagevalues(float) read-in from csv-datafile.
        chargeEnergy(list of float): List containing all Chargeenergyvalues(loat) read-in from csv-datafile.
        dischargeEnergy(list of float): List containing all Dischargeenergyvalues(float) read-in from csv-datafile.
        analysisDict(Dictioanry): Global Dictionary to save all raw and split up Data
        cyclenamelist(list of int): List containing the limiting points of the charge/discharge cycles also including CV step.
        ------------------

        Returns
        ------------------
       
        ------------------
        
        '''
        for i in range(len(cycleList)):
            cycletesttime = []
            if i == 0:
                minlimit = 0
            else:
                minlimit = cycleList[i-1]
            correction = testTime[minlimit]
            for l in range(minlimit,cycleList[i]):
                cycletesttime.append(testTime[l]-correction)
            currentlist = []
            voltlist = []
            chargeenergylist = []
            dischargeenergylist= []
            cycleend = cycleList[i]
            cyclestart = cycleList[i-1]
            if i == 0:
                cycleend = cycleList[0]
                cyclestart = 0
            chargeEdiff = chargeEnergy[cyclestart]
            dchargeEdiff = dischargeEnergy[cyclestart]
            for l in range(cyclestart,cycleend):
                currentlist.append(current[l])
                voltlist.append(voltage[l])
                chargeenergylist.append(abs(chargeEnergy[l] - chargeEdiff))
                dischargeenergylist.append(abs(dischargeEnergy[l]- dchargeEdiff))
            capacity = np.trapz(currentlist,cycletesttime) * (1000/3600)
            capacityrounded = np.around(capacity,4)
            try:
                if currentlist[0] != currentlist[1]:
                    del currentlist[0]
                    del voltlist[0]
                    del chargeenergylist[0]
                    del dischargeenergylist[0]
                    del cycletesttime[0]
                if currentlist[-1] != currentlist[-2]:
                    del currentlist[-1]
                    del voltlist[-1]
                    del chargeenergylist[-1]
                    del dischargeenergylist[-1]
                    del cycletesttime[-1]
            except:
                pass
            analysisDict['split'][cyclenamelist[i]]['t'] = np.array(cycletesttime)
            analysisDict['split'][cyclenamelist[i]]['I'] = np.array(currentlist) 
            analysisDict['split'][cyclenamelist[i]]['V'] = np.array(voltlist)
            analysisDict['split'][cyclenamelist[i]]['chargeE'] = np.array(chargeenergylist)
            analysisDict['split'][cyclenamelist[i]]['dischargeE'] = np.array(dischargeenergylist)
            analysisDict['split'][cyclenamelist[i]]['C'] = [capacityrounded]
            #analysisDict['split']['Capacities(mAh)'].append(capacityrounded)

def saveDatatoDic(testTime:list,stepindex:list,current:list, voltage: list, chargeEnergy:list, dischargeEnergy:list, analysisDict:dict, cycleindexlist:list):
        '''
        Saves raw and within calcualted data of one test as numpy arrays to the analysisDict.

        Passed Parameters
        ------------------
        testTime(list of int): List containing increasing testtime(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(float) read-in from csv-datafile.
        voltage(list of float): List containing all Voltagevalues(float) read-in from csv-datafile.
        chargeEnergy(list of float): List containing all Chargeenergyvalues(loat) read-in from csv-datafile.
        dischargeEnergy(list of float): List containing all Dischargeenergyvalues(float) read-in from csv-datafile.
        analysisDict(Dictioanry): Global Dictionary to save all raw and split up Data
        cycleindexlist(list of int): List contaoning all Cycleinideces read in from csv file.
        ------------------
                
        Returns
        ------------------

        ------------------
        '''
        analysisDict['raw']['Cycle_Index'] = np.array(cycleindexlist)
        analysisDict['raw']['Step_Index'] = np.array(stepindex)
        analysisDict['raw']['Testtime(s)'] = np.array(testTime)
        analysisDict['raw']['Voltage(V)'] = np.array(voltage)
        analysisDict['raw']['Current(A)'] = np.array(current)
        analysisDict['raw']['charge_Energy(Wh)'] = np.array(chargeEnergy)
        analysisDict['raw']['discharge_Energy(Wh)'] = np.array(dischargeEnergy)

def changeDatatype(datadict:dict):
        '''
        Changes data to lists and floats.

        Passed Parameters
        ------------------
        datadict(nested dict): Containing all data from a hdf5file.
        ------------------
                
        Returns
        ------------------

        ------------------
        '''
        for key in datadict:
            if key == 'raw':
                for cycle in datadict[key]:
                    o = list(datadict[key][cycle])
                    datadict[key][cycle] = o
                    for i in range(len(datadict[key][cycle])):
                        datadict[key][cycle][i] = datadict[key][cycle][i].astype(np.float64)
            if key == 'split':
                for cycle in (datadict[key]):
                    for param in (datadict[key][cycle]):
                        try:
                            l = list(datadict[key][cycle][param])
                            datadict[key][cycle][param] = l
                            for i in range(len(datadict[key][cycle][param])):
                                datadict[key][cycle][param][i] = datadict[key][cycle][param][i].astype(np.float64)
                        except:
                            o = list(datadict[key][cycle])
                            datadict[key][cycle] = o
                            for i in range(len(datadict[key][cycle])):
                                datadict[key][cycle][i] = datadict[key][cycle][i].astype(np.float64)
    
def save_dict_to_hdf5(dic, filename, hdf5path):
        '''
        Saves a dictioanry to a hdf5file.

        This function was copied from stackoverflow.
        '''
        path = os.getcwd()
        os.chdir(hdf5path)
        with h5py.File(filename, 'a') as h5file:
            recursively_save_dict_contents_to_group(h5file, '/', dic)
        os.chdir(path)

def recursively_save_dict_contents_to_group(h5file, path, dic):
            '''
            Saves dictionary content to groups.

            This function was copied from stackoverflow.
            '''
            # argument type checking
            if not isinstance(dic, dict):
                raise ValueError("must provide a dictionary")        
            if not isinstance(path, str):
                raise ValueError("path must be a string")
            if not isinstance(h5file, h5py._hl.files.File):
                raise ValueError("must be an open h5py file")
            # save items to the hdf5 file
            for key, item in dic.items():
                #print(key,item)
                key = str(key)
                if isinstance(item, list):
                    item = np.array(item)
                    #print(item)
                if not isinstance(key, str):
                    raise ValueError("dict keys must be strings to save to hdf5")
                # save strings, numpy.int64, and numpy.float64 types
                if isinstance(item, (np.int64, np.float64, str, float, float, float,int)):
                    #print( 'here' )
                    h5file[path + key] = item
                    #print(h5file[path + key])
                    #print(item)
                    if not h5file[path + key].value == item:
                        raise ValueError('The data representation in the HDF5 file does not match the original dict.')
                # save numpy arrays
                elif isinstance(item, np.ndarray):            
                    try:
                        h5file[path + key] = item
                    except:
                        item = np.array(item).astype('|S32')      # S32 defines you length of reserved diskspace and max number of letters
                        h5file[path + key] = item
                    #if not np.array_equal(h5file[path + key].value, item):
                    #   raise ValueError('The data representation in the HDF5 file does not match the original dict.')
                # save dictionaries
                elif isinstance(item, dict):
                    recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
                # other types cannot be saved and will result in an error
                else:
                    #print(item)
                    raise ValueError('Cannot save %s type.' % type(item))


def EOLC58HDF5files():

    pathdata = r"C:\Users\Leon Fischer\Desktop\PhD\Data\INFORM\EOL 20C\EOL eff C20 pulsed batch 1"

    hdf5path = pathdata

    os.chdir(pathdata)
    folder = os.listdir()
    namecriteria = 'test'

    for content in folder:
        if namecriteria in content:
            print(content)
            origpath = os.getcwd()
            os.chdir(content)
            #print(content)
            hdf5name = content.split('_')[0]+'_' +content.split('_')[1]+'_' +content.split('_')[2] + '.hdf5'
            data = getChannelData()
            os.chdir(origpath)
            testTime,cycleIndex,stepIndex,I,V,chargeEnergy,dischargeEnergy = serperaterawData(data)
            cyclelist = getCyclelimitListEOLPuls(cycleIndex,stepIndex,I)
            analysisDict={}
            cyclenamelist = []
            analysisDict['raw'] = {}
            analysisDict['split'] = {}
            for i in range(1,len(cyclelist)+1):
                analysisDict['split']['Step_'+str(i)] = {} 
                cyclenamelist.append('Step_'+str(i))
            print(len(cyclelist))
            getandsafeCycledata(cyclelist,testTime,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclenamelist )
            for i in range(1,len(analysisDict['split'].keys())):
                print('Step_'+str(i))
                t = analysisDict['split']['Step_'+str(i)]['t']
                v = analysisDict['split']['Step_'+str(i)]['V']
                plt.plot(t,v)
                plt.show()
            saveDatatoDic(testTime,stepIndex,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclelist)
            #print(analysisDict['raw'].keys())
            #t = analysisDict['raw']['Testtime(s)']
            #v = analysisDict['raw']['Voltage(V)']
            changeDatatype(analysisDict)  # Error of ValueError: [TypeError('cannot convert dictionary update sequence element #0 to a sequence') is because of datatype np.arra in beginning
            save_dict_to_hdf5(analysisDict,hdf5name, hdf5path)
            analysisDict.clear()
            cyclenamelist.clear()
            os.chdir(origpath)

def EOLC58Longtermfile():
    pathdata = r"C:\Users\Leon Fischer\Desktop\PhD\Data\INFORM\EOL 20C\C5 pulse\Longterm"

    hdf5path = r"C:\Users\Leon Fischer\Desktop\PhD\Data\INFORM\EOL 20C\C5 pulse\Longterm"

    os.chdir(pathdata)
    folder = os.listdir()
    namecriteria = 'InFormPulsedC5'

    for content in folder:
        if namecriteria in content:
            print(content)
            origpath = os.getcwd()
            os.chdir(content)
            #print(content)
            hdf5name = content.split('_')[0]+'_' +content.split('_')[1]+'_' +content.split('_')[2] + '.hdf5'
            data = getChannelData()
            os.chdir(origpath)
            testTime,cycleIndex,stepIndex,I,V,chargeEnergy,dischargeEnergy = serperaterawData(data)
            cyclelist = getCyclelimitListEOLPulsLongterm(cycleIndex,stepIndex,I)
            analysisDict={}
            cyclenamelist = []
            analysisDict['raw'] = {}
            analysisDict['split'] = {}
            for i in range(1,len(cyclelist)+1):
                analysisDict['split']['Step_'+str(i)] = {} 
                cyclenamelist.append('Step_'+str(i))
            print(len(cyclelist))
            getandsafeCycledata(cyclelist,testTime,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclenamelist )
            #for i in range(1,len(analysisDict['split'].keys())):
            #    print('Step_'+str(i))
            #    t = analysisDict['split']['Step_'+str(i)]['t']
            #    v = analysisDict['split']['Step_'+str(i)]['V']
            #    plt.plot(t,v)
            #    plt.show()
            saveDatatoDic(testTime,stepIndex,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclelist)
            #print(analysisDict['raw'].keys())
            #t = analysisDict['raw']['Testtime(s)']
            #v = analysisDict['raw']['Voltage(V)']
            changeDatatype(analysisDict)  # Error of ValueError: [TypeError('cannot convert dictionary update sequence element #0 to a sequence') is because of datatype np.arra in beginning
            save_dict_to_hdf5(analysisDict,hdf5name, hdf5path)
            analysisDict.clear()
            cyclenamelist.clear()
            os.chdir(origpath)

pathdata = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C20 CC - eff C20 pulsed ref\300Cycles"

hdf5path = pathdata
namecriteria = '_InForm'

def EOLCCCVHDF5files():
    os.chdir(pathdata)
    folder = os.listdir()
    for content in folder:
        if namecriteria in content:
            print(content)
            origpath = os.getcwd()
            os.chdir(content)
            #print(content)
            hdf5name = content.split('_')[0]+'_' +content.split('_')[1]+'_' +content.split('_')[2] + '.hdf5'
            data = getChannelData()
            os.chdir(origpath)
            testTime,cycleIndex,stepIndex,I,V,chargeEnergy,dischargeEnergy = serperaterawData(data)
            cyclelist = getCyclelimitListEOLCCCV(cycleIndex,stepIndex,I)
            analysisDict={}
            cyclenamelist = []
            analysisDict['raw'] = {}
            analysisDict['split'] = {}
            for i in range(1,len(cyclelist)+1):
                analysisDict['split']['Step_'+str(i)] = {} 
                cyclenamelist.append('Step_'+str(i))
            print(len(cyclelist))
            getandsafeCycledata(cyclelist,testTime,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclenamelist )
            #for i in range(1,len(analysisDict['split'].keys())):
            #    print('Step_'+str(i))
            #    t = analysisDict['split']['Step_'+str(i)]['t']
            #    v = analysisDict['split']['Step_'+str(i)]['V']
            #    plt.plot(t,v)
            #    plt.show()
            saveDatatoDic(testTime,stepIndex,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclelist)
            #print(analysisDict['raw'].keys())
            #t = analysisDict['raw']['Testtime(s)']
            #v = analysisDict['raw']['Voltage(V)']
            changeDatatype(analysisDict)  # Error of ValueError: [TypeError('cannot convert dictionary update sequence element #0 to a sequence') is because of datatype np.arra in beginning
            save_dict_to_hdf5(analysisDict,hdf5name, hdf5path)
            analysisDict.clear()
            cyclenamelist.clear()
            os.chdir(origpath)

def EOLCCCCVLongtermfile():

    os.chdir(pathdata)
    folder = os.listdir()
    for content in folder:
        if namecriteria in content:
            print(content)
            origpath = os.getcwd()
            os.chdir(content)
            #print(content)
            hdf5name = content.split('_')[0]+'_' + content.split('_')[1]+'_' +content.split('_')[2]+'_new' + '.hdf5'
            data = getChannelData()
            os.chdir(origpath)
            testTime,cycleIndex,stepIndex,I,V,chargeEnergy,dischargeEnergy = serperaterawData(data)
            cyclelist = getCyclelimitListEOLCCCVLongterm(cycleIndex,stepIndex,I)
            analysisDict={}
            cyclenamelist = []
            analysisDict['raw'] = {}
            analysisDict['split'] = {}
            for i in range(1,len(cyclelist)+1):
                analysisDict['split']['Step_'+str(i)] = {} 
                cyclenamelist.append('Step_'+str(i))
            print(len(cyclelist))
            getandsafeCycledata(cyclelist,testTime,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclenamelist )
            #for i in range(1,len(analysisDict['split'].keys())):
            #    current = analysisDict['split']['Step_'+str(i)]['I']
            #    print('Step_'+str(i))
            #    print(current[0])
            #    t = analysisDict['split']['Step_'+str(i)]['t']
            #    v = analysisDict['split']['Step_'+str(i)]['V']
            #    plt.plot(t,v)
            #    plt.show()
            saveDatatoDic(testTime,stepIndex,I,V,chargeEnergy,dischargeEnergy,analysisDict,cyclelist)
            #print(analysisDict['raw'].keys())
            #t = analysisDict['raw']['Testtime(s)']
            #v = analysisDict['raw']['Voltage(V)']
            changeDatatype(analysisDict)  # Error of ValueError: [TypeError('cannot convert dictionary update sequence element #0 to a sequence') is because of datatype np.arra in beginning
            save_dict_to_hdf5(analysisDict,hdf5name, hdf5path)
            analysisDict.clear()
            cyclenamelist.clear()
            os.chdir(origpath)




######################### ZSW data

def getChannelDataZSW():
        def atof(text):
            try:
                retval = float(text)
            except ValueError:
                retval = text
            return retval
        def natural_keys(text):
            return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

        '''
        Reading in csv file of the testname, poping some parameters to avaoid discspace problems and saving remaining as a dictionary.

        Passed Parameters
        ------------------
        datapath(string): Path to data choosen before when using getpossibletestnames()
        ------------------

        Returns
        ------------------
        datadict(dictionary): Dictionary containing all read-in data from the csv outputfile from arbin.
        ------------------
        '''
        content = os.listdir()
        #print(str(content))
        df = pd.read_csv(str(content[0]),index_col=False, delimiter=';')
        '''l =df.pop('Schritt')
        l =df.pop('Status')
        l =df.pop('Zeitstempel')
        l =df.pop('Schritt Zeit')
        l =df.pop('Zyklus')
        l =df.pop('Zyklusebene')
        l =df.pop('Prozedur')
        l =df.pop('AhAkku')
        l =df.pop('AhStep')
        l =df.pop('Energie')
        l =df.pop('WhStep')
        l =df.pop('PreOhm')
        l =df.pop('Ri_Pulse')
        l =df.pop('TypK')
        l =df.pop('AliveCnt')
        l =df.pop('Rin')
        l =df.pop('Temperatur')'''
        l =df.pop('time')
        l =df.pop('cycle')
        l =df.pop('channel')
        l =df.pop('cell_temp1')
        l =df.pop('cell_temp2')
        l =df.pop('charge_count_pos')
        l =df.pop('charge_count_neg')
        l =df.pop('charge_count_sum')
        l =df.pop('energy_count_pos')
        l =df.pop('energy_count_neg')
        l =df.pop('energy_count_sum')
        l =df.pop('power')
        l =df.pop('resistance')
        l =df.pop('soc')
        l =df.pop('soc_approx')
        #df.rename(columns={'Progr. Zeit': 'time', 'Spannung':"voltage", 'Strom': 'current'}, inplace=True)
        df.rename(columns={'offset': 'time', 'stepcnt':'stepindex'}, inplace=True)
        time_list= []
        for volt in df['time']:
            t = volt/1000
            time_list.append(float(t))
        series = pd.Series(time_list)
        df['time'] = series
        ''' volt_list= []
        for volt in df['voltage']:
            t = volt.replace(",", ".")
            volt_list.append(float(t))
        series = pd.Series(volt_list)
        df['voltage'] = series
        curr_list= []
        for curr in df['current']:
            t = curr.replace(",", ".")
            curr_list.append(float(t)/1000)
        series = pd.Series(curr_list)
        df['current'] = series'''
        #print(df)   
        #df['Zeitstempel'] = pd.to_datetime(df['Zeitstempel'])
        # Define the reference start time
        #reference_time = pd.to_datetime(df['Zeitstempel'][0])
        #df['Zeitstempel'] = np.round((df['Zeitstempel'] - reference_time).dt.total_seconds(),2)
        #df = df.astype(np.float64)  84744.132
        datadict = df.to_dict()
        if len(content) > 1:
                print('There are more csv files in.')
                for i in range(1,len(content)):
                    print(content[i])
                    df2 = pd.read_csv(str(content[i]),index_col=False, delimiter=';')
                    l = df2.pop('cycle')
                    l = df2.pop('offset')
                    l = df2.pop('channel')
                    l = df2.pop('charge_count_neg')
                    l = df2.pop('charge_count_pos')
                    l = df2.pop('power')
                    l = df2.pop('energy_count_sum')
                    l = df2.pop('energy_count_neg')
                    l = df2.pop('energy_count_pos')
                    l = df2.pop('charge_count_sum')
                    l = df2.pop('resistance')
                    l = df2.pop('soc')
                    l = df2.pop('soc_approx')
                    l = df2.pop('cell_temp1')
                    l = df2.pop('cell_temp2')
                    l = df2.pop('time')
                    l = df2.pop('stepcnt')
                    df2.rename(columns={'time.1': 'time',}, inplace=True)
                    #df2['time'] = pd.to_datetime(df2['time'])
                    #reference_time = pd.to_datetime(df2['time'][0]) 
                    try:
                        addtime = list(df['time'][df.index[-1]])[-1]
                    except:
                        addtime = df['time'][list(df.index)[-1]]
                    time_list= []
                    for time in df2['time']:
                        t = time.replace(",", ".")
                        time_list.append((float(t)*3600)+addtime)
                    series = pd.Series(time_list)
                    df2['time'] = series
                    #df2['time'] = np.round((df2['time'] - reference_time).dt.total_seconds(),2)+addtime
                    df2.drop(0)
                    df = pd.concat([df, df2])
                df.set_index('time', inplace = True, drop=False)
                df = df.astype(np.float64)
                completeDict = df.to_dict()
                return completeDict
        else:
            return datadict
       
def serperaterawDataZSW(data:dict):
        '''
        Serperates and saves the raw data saved in a dictionary to lists of each in future used parameter.

        Passed Parameters
        ------------------
        data(Dictioanry): Dictionary containing all information from the read in csv result file from arbin.
        ------------------

        Returns
        ------------------
        testTime(list of int): List containing increasing testtime(int) read-in from csv-datafile.
        cycleIndex(list of int): List containing the corresponging cycleindex to each datapoint read-in from csv-file.
        I(list of float): List containing all Currentvalues(float) read-in from csv-datafile.
        V(list of float): List containing all Voltagevalues(float) read-in from csv-datafile.
        chargeEnergy(list of float): List containing all Chargeenergyvalues(loat) read-in from csv-datafile.
        dischargeEnergy(list of float): List containing all Dischargeenergyvalues(float) read-in from csv-datafile.
        ------------------
        '''
        testTime=list(data['time'].values())
        #datapoint = list(data['Data_Point'].values())                
        stepIndex = list(data['stepindex'].values())
        I = list(data['current'].values())
        V = list(data['voltage'].values())
        #cell_temp1 = list(data['cell_temp1'].values())
        #cell_temp2 = list(data['cell_temp2'].values())
        #chargeCapacity = list(data['Charge_Capacity(Ah)'].values())
        #dischargeCapacity = list(data['Discharge_Capacity(Ah)'].values())
        return  testTime,I,V, stepIndex

def check_cyclelist(cyclelist, testTime, voltage, current):
    new_cyclelist = []
    for i in range(len(cyclelist)):
        cycletesttime = []
        if i == 0:
            minlimit = 0
        else:
            minlimit = cyclelist[i-1]
        correction = testTime[minlimit]
        for l in range(minlimit,cyclelist[i]):
            cycletesttime.append(testTime[l]-correction)
        currentlist = []
        voltlist = []
        cycleend = cyclelist[i]
        cyclestart = cyclelist[i-1]
        if i == 0:
            cycleend = cyclelist[0]
            cyclestart = 0
        for l in range(cyclestart,cycleend):
            currentlist.append(current[l])
            voltlist.append(voltage[l])
        if len(voltlist) >=4:
            new_cyclelist.append(cyclelist[i])
    return new_cyclelist

def getandsafeCycledataZSW(cycleList:list,testTime:list,current:list, voltage: list,analysisDict:dict, cyclenamelist:list):
        '''
        Slices and saves data for each cycle.
        
        
        Passed Parameters
        ------------------
        cycleList(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        testTime(list of int): List containing increasing testtime(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(float) read-in from csv-datafile.
        voltage(list of float): List containing all Voltagevalues(float) read-in from csv-datafile.
        chargeEnergy(list of float): List containing all Chargeenergyvalues(loat) read-in from csv-datafile.
        dischargeEnergy(list of float): List containing all Dischargeenergyvalues(float) read-in from csv-datafile.
        analysisDict(Dictioanry): Global Dictionary to save all raw and split up Data
        cyclenamelist(list of int): List containing the limiting points of the charge/discharge cycles also including CV step.
        ------------------

        Returns
        ------------------
       
        ------------------
        
        '''
        for i in range(len(cycleList)):
            cycletesttime = []
            if i == 0:
                minlimit = 0
            else:
                minlimit = cycleList[i-1]
            correction = testTime[minlimit]
            for l in range(minlimit,cycleList[i]):
                cycletesttime.append(testTime[l]-correction)
            currentlist = []
            voltlist = []
            cycleend = cycleList[i]
            cyclestart = cycleList[i-1]
            if i == 0:
                cycleend = cycleList[0]
                cyclestart = 0
            for l in range(cyclestart,cycleend):
                currentlist.append(current[l])
                voltlist.append(voltage[l])
            capacity = np.trapz(currentlist,cycletesttime) * (1000/3600)
            capacityrounded = np.around(capacity,4)
            try:
                if currentlist[0] != currentlist[1]:
                    del currentlist[0]
                    del voltlist[0]
                    del cycletesttime[0]
                if currentlist[-1] != currentlist[-2]:
                    del currentlist[-1]
                    del voltlist[-1]
                    del cycletesttime[-1]
            except:
                pass
            analysisDict['split'][cyclenamelist[i]]['t'] = np.array(cycletesttime)
            analysisDict['split'][cyclenamelist[i]]['I'] = np.array(currentlist) 
            analysisDict['split'][cyclenamelist[i]]['V'] = np.array(voltlist)
            analysisDict['split'][cyclenamelist[i]]['C'] = [capacityrounded]
            #analysisDict['split']['Capacities(mAh)'].append(capacityrounded)

def saveDatatoDic(testTime:list,current:list, voltage: list, analysisDict:dict):
        '''
        Saves raw and within calcualted data of one test as numpy arrays to the analysisDict.

        Passed Parameters
        ------------------
        testTime(list of int): List containing increasing testtime(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(float) read-in from csv-datafile.
        voltage(list of float): List containing all Voltagevalues(float) read-in from csv-datafile.
        chargeEnergy(list of float): List containing all Chargeenergyvalues(loat) read-in from csv-datafile.
        dischargeEnergy(list of float): List containing all Dischargeenergyvalues(float) read-in from csv-datafile.
        analysisDict(Dictioanry): Global Dictionary to save all raw and split up Data
        cycleindexlist(list of int): List contaoning all Cycleinideces read in from csv file.
        ------------------
                
        Returns
        ------------------

        ------------------
        '''
        #analysisDict['raw']['Cycle_Index'] = np.array(cycleindexlist)
        #analysisDict['raw']['Step_Index'] = np.array(stepindex)
        analysisDict['raw']['Testtime(s)'] = np.array(testTime)
        analysisDict['raw']['Voltage(V)'] = np.array(voltage)
        analysisDict['raw']['Current(A)'] = np.array(current)
        #analysisDict['raw']['cell_temp1'] = np.array(cell_temp1)
        #analysisDict['raw']['cell_temp2'] = np.array(cell_temp2)

def changeDatatype(datadict:dict):
        '''
        Changes data to lists and floats.

        Passed Parameters
        ------------------
        datadict(nested dict): Containing all data from a hdf5file.
        ------------------
                
        Returns
        ------------------

        ------------------
        '''
        for key in datadict:
            if key == 'raw':
                for cycle in datadict[key]:
                    o = list(datadict[key][cycle])
                    datadict[key][cycle] = o
                    for i in range(len(datadict[key][cycle])):
                        datadict[key][cycle][i] = datadict[key][cycle][i].astype(np.float64)
            if key == 'split':
                for cycle in (datadict[key]):
                    for param in (datadict[key][cycle]):
                        try:
                            l = list(datadict[key][cycle][param])
                            datadict[key][cycle][param] = l
                            for i in range(len(datadict[key][cycle][param])):
                                datadict[key][cycle][param][i] = datadict[key][cycle][param][i].astype(np.float64)
                        except:
                            o = list(datadict[key][cycle])
                            datadict[key][cycle] = o
                            for i in range(len(datadict[key][cycle])):
                                datadict[key][cycle][i] = datadict[key][cycle][i].astype(np.float64)
    
def save_dict_to_hdf5(dic, filename, hdf5path):
        '''
        Saves a dictioanry to a hdf5file.

        This function was copied from stackoverflow.
        '''
        path = os.getcwd()
        os.chdir(hdf5path)
        with h5py.File(filename, 'a') as h5file:
            recursively_save_dict_contents_to_group(h5file, '/', dic)
        os.chdir(path)

def recursively_save_dict_contents_to_group(h5file, path, dic):
            '''
            Saves dictionary content to groups.

            This function was copied from stackoverflow.
            '''
            # argument type checking
            if not isinstance(dic, dict):
                raise ValueError("must provide a dictionary")        
            if not isinstance(path, str):
                raise ValueError("path must be a string")
            if not isinstance(h5file, h5py._hl.files.File):
                raise ValueError("must be an open h5py file")
            # save items to the hdf5 file
            for key, item in dic.items():
                #print(key,item)
                key = str(key)
                if isinstance(item, list):
                    item = np.array(item)
                    #print(item)
                if not isinstance(key, str):
                    raise ValueError("dict keys must be strings to save to hdf5")
                # save strings, numpy.int64, and numpy.float64 types
                if isinstance(item, (np.int64, np.float64, str, float, float, float,int)):
                    #print( 'here' )
                    h5file[path + key] = item
                    #print(h5file[path + key])
                    #print(item)
                    if not h5file[path + key].value == item:
                        raise ValueError('The data representation in the HDF5 file does not match the original dict.')
                # save numpy arrays
                elif isinstance(item, np.ndarray):            
                    try:
                        h5file[path + key] = item
                    except:
                        item = np.array(item).astype('|S32')      # S32 defines you length of reserved diskspace and max number of letters
                        h5file[path + key] = item
                    #if not np.array_equal(h5file[path + key].value, item):
                    #   raise ValueError('The data representation in the HDF5 file does not match the original dict.')
                # save dictionaries
                elif isinstance(item, dict):
                    recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
                # other types cannot be saved and will result in an error
                else:
                    #print(item)
                    raise ValueError('Cannot save %s type.' % type(item))



def getCyclelimitListEOLPulsZSW(stepindex:list,current:list, voltlist:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        stepindexlist = []
        currentlist = []

        for l in range(len(current)-1):
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l]==9 and stepindex[l-1]!=9 
                or voltlist[l] <4.2 and voltlist[l-1]==4.2 and current[l]<1 and stepindex[l]>7 
                or voltlist[l] >4 and current[l]<-1 and current[l-1]>-1 
                or voltlist[l-1]<=2.9 and current[l-1]<-10 and current[l]>-1 and stepindex[l]>9
                or voltlist[l] >=3 and current[l-1]<1 and current[l]>2
                ):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        # 29 Steps seem right
        return cyclelimitlist

pathZSWEOL = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Reference\EOL_Referenzformierung'

os.chdir(pathZSWEOL)
hdf5path = pathZSWEOL
folder = os.listdir()
namecriteria = 'Zelle'
counter = 0

for content in folder:
    origpath = os.getcwd()
    os.chdir(content)
    print(content)
    hdf5name = str(content) + '.hdf5'
    #print(hdf5name)
    data = getChannelDataZSW()
    os.chdir(origpath)
    testTime,I,V, stepindex = serperaterawDataZSW(data)
    cyclelist = getCyclelimitListEOLPulsZSW(stepindex,I, V)
    cyclelist = check_cyclelist(cyclelist, testTime, V, I)
    analysisDict={}
    cyclenamelist = []
    analysisDict['raw'] = {}
    analysisDict['split'] = {}
    for i in range(1,len(cyclelist)+1):
        analysisDict['split']['Step_'+str(i)] = {} 
        cyclenamelist.append('Step_'+str(i))
    print(cyclelist)
    getandsafeCycledataZSW(cyclelist,testTime,I,V,analysisDict,cyclenamelist )
    #for i in range(1,len(analysisDict['split'].keys())):
    #    t = analysisDict['split']['Step_'+str(i)]['t']
    #    v = analysisDict['split']['Step_'+str(i)]['V']
    #    plt.plot(t,v)
    #    plt.show()
    saveDatatoDic(testTime,I,V,analysisDict)
    '''print(analysisDict['raw'].keys())
    t = analysisDict['raw']['Testtime(s)']
    v = analysisDict['raw']['Voltage(V)']'''
    changeDatatype(analysisDict)  # Error of ValueError: [TypeError('cannot convert dictionary update sequence element #0 to a sequence') is because of datatype np.arra in beginning
    save_dict_to_hdf5(analysisDict,hdf5name, hdf5path)
    analysisDict.clear()
    cyclenamelist.clear()
    os.chdir(origpath)



def getCyclelimitListFormPulsZSW(stepindex:list,current:list, voltlist:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        stepindexlist = []
        currentlist = []

        for l in range(len(current)-1):
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l] == 5 and stepindex[l-1]!=5 or currentlist[l]<-1 and currentlist[l-1]>-1 or currentlist[l-1]<-1 and abs(currentlist[l]-currentlist[l-1])>=1):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist


def getCyclelimitListLongtermZSW(stepindex:list,current:list, voltlist:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        stepindexlist = []
        currentlist = []

        for l in range(len(current)-1):
            stepindexlist.append(stepindex[l])
            currentlist.append(current[l])
            if (stepindex[l]==6 and stepindex[l-1]!=6 
                or voltlist[l] >=2.9 and current[l-1]<-1 and current[l]>-0.9 
                or voltlist[l] >3 and current[l-1]<1 and current[l]>1
                or voltlist[l] <4.2 and voltlist[l-1]>=4.2 and current[l]<1 and stepindex[l]>7 
                or voltlist[l] >4 and current[l]<-1 and current[l-1]>-1                                   
                ):
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist

def getCyclelimitListFormRefZSW(current:list, voltlist:list):
        '''
        Returns a list containing cyclelimits based on changes in stepindices and current.

        Passed Parameters
        ------------------
        cycleindex(list of int):  List containing all Cycleindices(int) read-in from csv-datafile.
        stepindex(list of int): List containing all Stepindices(int) read-in from csv-datafile.
        current(list of float): List containing all Currentvalues(loat) read-in from csv-datafile.
        ------------------

        Returns
        ------------------
        cycleborderlist(list of int): List containing Limitindices(including Indice to cycle)
        ------------------
        '''
        cyclelimitlist = []
        for l in range(len(current)-1):
            if current[l-1] > -1 and current[l]<-24 or voltlist[l]==2.9 and current[l]>-24 and current[l-1]<-24 and voltlist[l]>=2.9 or voltlist[l]>= 2.9 and current[l]> -1 and current[l-1]<-1 :
                cyclelimitlist.append(l-1)
        cyclelimitlist.append(len(current))   
        return cyclelimitlist


# formation  

#pathdata = r"C:\Users\Leon Fischer\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Pulsed\Pulsed_Formation"
#pathdata = r"C:\Users\Leon Fischer\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Pulsed\End of Line Test"
#pathdata = r'C:\Users\Leon Fischer\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Pulsed\Zyklisierung'
pathdata = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Reference\Formierung'
hdf5path = pathdata

os.chdir(pathdata)
folder = os.listdir()
namecriteria = 'Zelle'
counter = 0

"""
>>> print(len(regularcycles[0])) 
308
>>> print(len(regularcycles[1])) 
360
>>> print(len(regularcycles[2])) 
886
>>> print(len(regularcycles[3])) 
362
>>> print(len(regularcycles[4])) 
740
"""
for content in folder:
        if namecriteria in content:
            #print(content)
            origpath = os.getcwd()
            os.chdir(content)
            hdf5name = content.split('.')[0]+'.hdf5'
            print(hdf5name)
            data = getChannelDataZSW()
            os.chdir(origpath)
            testTime,I,V= serperaterawDataZSW(data)
            #plt.plot(testTime, V)
            #plt.show()
            cyclelist = getCyclelimitListFormRefZSW(I, V)
            #starts with rest
            analysisDict={}
            cyclenamelist = []
            analysisDict['raw'] = {}
            analysisDict['split'] = {}
            for i in range(1,len(cyclelist)+1):
                analysisDict['split']['Step_'+str(i)] = {} 
                cyclenamelist.append('Step_'+str(i))
            print(len(cyclelist))
            getandsafeCycledataZSW(cyclelist,testTime,I,V,analysisDict,cyclenamelist )
            #for i in range(1,len(analysisDict['split'].keys())+1):
            #    print('Step_'+str(i))
            #    t = analysisDict['split']['Step_'+str(i)]['t']
            #    v = analysisDict['split']['Step_'+str(i)]['V']
            #    curr =  analysisDict['split']['Step_'+str(i)]['I']
            #    plt.scatter(t,v)
                #plt.scatter(t,curr)
            #    plt.show()
            saveDatatoDic(testTime,I,V,analysisDict,cyclelist)
            #print(analysisDict['raw'].keys())
            t = analysisDict['raw']['Testtime(s)']
            v = analysisDict['raw']['Voltage(V)']
            #curr = analysisDict['raw']['Current(A)']
            labels=["1","2","3","4","5"]
            #plt.plot(t,v, label=labels[counter])
            #counter = counter +1
            #plt.legend()
            changeDatatype(analysisDict)  # Error of ValueError: [TypeError('cannot convert dictionary update sequence element #0 to a sequence') is because of datatype np.arra in beginning
            save_dict_to_hdf5(analysisDict,hdf5name, hdf5path)
            analysisDict.clear()
            cyclenamelist.clear()
            os.chdir(origpath)


            
