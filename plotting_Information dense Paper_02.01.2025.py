import os
import pandas as pd
import numpy as np
import hdfdict
import matplotlib.pyplot as plt
import scienceplots
import re
import copy
import matplotlib.image as image
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from matplotlib.ticker import StrMethodFormatter
from scipy.interpolate import interp1d
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
import matplotlib.lines as mlines
from scipy.stats import linregress
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter

def atof(text):
        try:
            retval = float(text)
        except ValueError:
            retval = text
        return retval

def natural_keys(text):
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

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


def loadHDF5( hdf5name:str):
        hdf5dict = dict(hdfdict.load(hdf5name, mode='r+',lazy=False))
        changeDatatype(hdf5dict)
        return hdf5dict



pathpulse = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\PulsedC58"
pathCCref = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\CCCV_C20"
#pathCCref = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C20 CC - eff C20 pulsed ref'
pathC57ref1 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 2 C58 CC20rest batch 1'
pathC57ref2= r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 1 C58 5s pulse batch 2'
pathC5pulse = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\C5 pulse'
pathC83CC = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\C823 CCCV'
pathdoublepulse = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL Double Pulse batch 1'
pathdoublepulse2 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL Double Pulse batch 2'
path5spulse = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 5s pulse batch 2'
pathC58CC20 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse CC20 rest batch 1'
plotfolder = r'C:\Users\Leon\Desktop\PhD\KIT\My papers\Formation Strategies'


def showPulseandCCCV():
    
    def getCCCVdata():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        voltlists = []
        timelists = []
        currentlists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)

        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                    voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                    timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                    currentlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                    currentlist = data['split']["Step_"+str(form)]['I']
                    timelist = data['split']["Step_"+str(form)]['t']
                    voltlist = data['split']["Step_"+str(form)]['V']                
                voltlists.append(voltlist)
                timelists.append(timelist)
                currentlists.append(currentlist)       
        return voltlists, timelists, currentlists
    
    def getpulsedata():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        voltlists = []
        timelists = []
        currentlists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)

        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                    voltlist = list(data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)'])
                    timelist = list(data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)'])
                    currentlist = list(data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)'])
                except:
                    currentlist = list(data['split']["Step_"+str(form)]['I'])
                    timelist = list(data['split']["Step_"+str(form)]['t'])
                    voltlist = list(data['split']["Step_"+str(form)]['V'])                
                voltlists.append(voltlist)
                timelists.append(timelist)
                currentlists.append(currentlist)       
        return voltlists, timelists, currentlists

    voltref, timeref, curref = getCCCVdata()
    voltpulse, timepulse, currpulse = getpulsedata()

    fig,ax = plt.subplots(1, 1, figsize=(5, 3), sharex=True)
    plt.style.use(['science','nature','no-latex'])

    ax.plot(timeref[0], curref[0], color = 'r', linewidth = 2)
    ax.plot(timepulse[0], currpulse[0], color = 'b', linewidth = 2)
    ax.set_ylabel('Current', fontsize=8)
    ax.set_xlabel('Time', fontsize=8)
    ax.set_xlim([0, 10])
    ax.set_ylim([-0.00000000001, 0.0007])
    textsize = 12
    ax.text(1, 0.0005, r'$t_{on}$', fontsize = textsize)
    ax.text(3.5, 0.0005, r'$t_{off}$', fontsize = textsize)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.legend(['Constant Current CC', 'Pulsed'])
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('PulseExplanation.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def getEOLcoulombeffPlotForm():

    # idee auf der X-achse die mediane Formierzeit anzeigen und entsprechen ordnen

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefnew():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                try:
                    refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                    refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                    refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                    refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                    formtimelist.append(data['split']["Step_1"]["t"])
                except:
                    print(i)

        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                    except:
                        print(i)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                    except:
                        print(i)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))  
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    if capacities[5] == 0:
                        del capacities[5]
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                    except:
                        print(i)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))) 
                        formtimelist.append(data['split']["Step_1"]["t"])  
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0])))   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa


    refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtime, eolcapalist = getCCref()
    allrefcoulomb = [refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3]
    allrefformtime = []
    for i in range(len(formtime)):
        allrefformtime.append(formtime[i][-1])

    refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8, formtimeC823, eolcapalistC823 = getCCrefC83()
    allrefcoulombC8 = [refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8]
    allrefformtimeC8 = []
    for i in range(len(formtimeC823)):
        allrefformtimeC8.append(formtimeC823[i][-1])

    pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimepulse, eolcapalistpulse = getpulseform()
    allpulsecoulomb = [pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3]
    allpulseformtime = []
    for i in range(len(formtimepulse)):
        allpulseformtime.append(formtimepulse[i][-1])

    pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5, formtimepulseC5, eolcapalistC5 = getpulseformC5()
    allpulsecoulombC5 = [pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5]
    allpulseformtimeC5 = []
    for i in range(len(formtimepulseC5)):
        allpulseformtimeC5.append(formtimepulseC5[i][-1])

    pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s, formtimepulseC575s, eolcapalistC575s = getpulseformC575s()
    allpulsecoulombC5s = [pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s]
    allpulseformtimeC575s = [] 
    for i in range(len(formtimepulseC575s)):
        allpulseformtimeC575s.append(formtimepulseC575s[i][-1])

    pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse, formtimedoublepulse, eolcapalistdoubleulse = getpulseformdoublepulse()
    allpulsecoulombdoublepulse = [pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse]
    allpulseformtimedoublepulse = []
    for i in range(len(formtimedoublepulse)):
        allpulseformtimedoublepulse.append(formtimedoublepulse[i][-1])

    pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20, formtimepulseC57CC20, eolcapalistC57CC20 = getpulseformC57CC20()
    allpulsecoulombC57CC20 = [pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20]
    allpulseformtimeC57CC20 = []
    for i in range(len(formtimepulseC57CC20)):
        allpulseformtimeC57CC20.append(formtimepulseC57CC20[i][-1])

    print(allrefformtime)


    meanformtime = []

    decround = 1

    meanformtime.append(np.round(np.mean(allrefformtime)/3600, decround))

    meanformtime.append(np.round(np.mean(allpulseformtimedoublepulse)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allrefformtimeC8)/3600, decround))
    
    meanformtime.append(np.round(np.mean(allpulseformtimeC5)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC57CC20)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC575s)/3600, decround))
    

    print(meanformtime)


    xlabels = []
    for i in range(len(meanformtime)):
        xlabels.append(str(meanformtime[i]))
    print(xlabels)

    #t_stat, p_val = ttest_ind(refcoulomb4Form, pulsecoulomb4Form)
    #print('t-statistic:', t_stat)
    #print('p-value:', p_val)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1, figsize=(6, 5))#, sharex=True)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)

    # Add Capacheck 1
    #alphavalue = 0.25
    #CCrefCC=ax.boxplot(allrefcoulomb[1], positions = [0], patch_artist=True, boxprops=dict(facecolor='r', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    #pulseCC = ax.boxplot(allpulsecoulomb[1], positions = [0.5], patch_artist=True, boxprops=dict(facecolor='b', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #CCrefC83CC=ax.boxplot(allrefcoulombC8[1], positions = [1], patch_artist=True, boxprops=dict(facecolor='darkred', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulseC5sCC = ax.boxplot(allpulsecoulombC5s[1], positions = [1.5], patch_artist=True, boxprops=dict(facecolor='goldenrod', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulseC57CC20CC = ax.boxplot(allpulsecoulombC57CC20[1], positions = [2.0], patch_artist=True, boxprops=dict(facecolor='dimgrey', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulsedoublepulsCC = ax.boxplot(allpulsecoulombdoublepulse[1], positions = [2.5], patch_artist=True, boxprops=dict(facecolor='darkgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulseC5CC = ax.boxplot(allpulsecoulombC5[1], positions = [3], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))

    x = [0]
    y = [0]
    for o in range(len(x)):
        CCref=ax.boxplot(allrefcoulomb[o], positions = [0], patch_artist=True, boxprops=dict(facecolor='r'))

        pulse = ax.boxplot(allpulsecoulomb[o], positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
        pulsedoublepuls = ax.boxplot(allpulsecoulombdoublepulse[o], positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen')) 
        CCrefC83=ax.boxplot(allrefcoulombC8[o], positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))
    
        pulseC5 = ax.boxplot(allpulsecoulombC5[o], positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
        pulseC57CC20 = ax.boxplot(allpulsecoulombC57CC20[o], positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
        pulseC5s = ax.boxplot(allpulsecoulombC5s[o], positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))
        
         
        
    legendcolors = [CCref["boxes"][0], pulse["boxes"][0],pulsedoublepuls["boxes"][0],   CCrefC83["boxes"][0], pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    #legendtext = ['N(CC20) = ' + str(len(refcoulomb4Form)),  'N(PulsC5.7) = ' + str(len(pulsecoulomb4Form)),'N(CC8.23) = ' + str(len(refcoulomb4FormC8)),
    #                   'N(PulsC5.7 5s) = ' + str(len(pulsecoulomb4FormC5s)), 'N(PulsC5.7 CC20) = ' + str(len(pulsecoulomb4FormC57CC20)),
    #                    'N(PulsC5.7double) = ' + str(len(pulsecoulomb4Formdoublepulse)), 'N(PulsC5) = ' + str(len(pulsecoulomb4FormC5))]
    legendtext = ['Reference CC C/20', 'Pulsed C/5.7' ,'Pulsed double length',  'CC C/8.23', 'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s']
    legend = plt.legend(legendcolors, legendtext, bbox_to_anchor=(0.14, 0, 0.6, 0.24), loc="upper left",
                 borderaxespad=0, fancybox=True, shadow=True,frameon=True, ncol = 2)

    ax.set_xlim([-0.25,4.25])
    #ax.set_ylim([0.76,0.86])

    ax.set_xticks([0, 1, 1.5, 2,  3, 3.5, 4]) 
    ax.set_xticklabels(xlabels)
        
    ax.set_ylabel('Coulombic efficiency [1]', fontsize=8)
    ax.set_xlabel('Mean charge formation time [h]', fontsize=8)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('CoulombEfficienymorePulsenFormation5cells.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def getEOLcoulombeffPlotFormBrokenAxis():

    # idee auf der X-achse die mediane Formierzeit anzeigen und entsprechen ordnen

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                    except:
                        print(i)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                    except:
                        print(i)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)  
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    if capacities[5] == 0:
                        del capacities[5]
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                    except:
                        print(i)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100) 
                        formtimelist.append(data['split']["Step_1"]["t"])  
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa


    refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtime, eolcapalist = getCCref()
    allrefcoulomb = [refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3]
    allrefformtime = []
    for i in range(len(formtime)):
        allrefformtime.append(formtime[i][-1])

    refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8, formtimeC823, eolcapalistC823 = getCCrefC83()
    allrefcoulombC8 = [refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8]
    allrefformtimeC8 = []
    for i in range(len(formtimeC823)):
        allrefformtimeC8.append(formtimeC823[i][-1])

    pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimepulse, eolcapalistpulse = getpulseform()
    allpulsecoulomb = [pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3]
    allpulseformtime = []
    for i in range(len(formtimepulse)):
        allpulseformtime.append(formtimepulse[i][-1])

    pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5, formtimepulseC5, eolcapalistC5 = getpulseformC5()
    allpulsecoulombC5 = [pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5]
    allpulseformtimeC5 = []
    for i in range(len(formtimepulseC5)):
        allpulseformtimeC5.append(formtimepulseC5[i][-1])

    pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s, formtimepulseC575s, eolcapalistC575s = getpulseformC575s()
    allpulsecoulombC5s = [pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s]
    allpulseformtimeC575s = [] 
    for i in range(len(formtimepulseC575s)):
        allpulseformtimeC575s.append(formtimepulseC575s[i][-1])

    pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse, formtimedoublepulse, eolcapalistdoubleulse = getpulseformdoublepulse()
    allpulsecoulombdoublepulse = [pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse]
    allpulseformtimedoublepulse = []
    for i in range(len(formtimedoublepulse)):
        allpulseformtimedoublepulse.append(formtimedoublepulse[i][-1])

    pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20, formtimepulseC57CC20, eolcapalistC57CC20 = getpulseformC57CC20()
    allpulsecoulombC57CC20 = [pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20]
    allpulseformtimeC57CC20 = []
    for i in range(len(formtimepulseC57CC20)):
        allpulseformtimeC57CC20.append(formtimepulseC57CC20[i][-1])

    print(allrefformtime)


    meanformtime = []
    std = []

    decround = 1

    meanformtime.append(np.round(np.mean(allrefformtime)/3600, decround))

    meanformtime.append(np.round(np.mean(allpulseformtimedoublepulse)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allrefformtimeC8)/3600, decround))
    
    meanformtime.append(np.round(np.mean(allpulseformtimeC5)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC57CC20)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC575s)/3600, decround))
    
    std.append(np.round(np.std(allrefformtime)/3600, decround))

    std.append(np.round(np.std(allpulseformtimedoublepulse)/3600, decround))
    std.append(np.round(np.std(allpulseformtime)/3600, decround))
    std.append(np.round(np.std(allrefformtimeC8)/3600, decround))

    std.append(np.round(np.std(allpulseformtimeC5)/3600, decround))
    std.append(np.round(np.std(allpulseformtimeC57CC20)/3600, decround))
    std.append(np.round(np.std(allpulseformtimeC575s)/3600, decround))

    print(meanformtime)
    print(std)

    xlabels = []
    for i in range(len(meanformtime)):
        xlabels.append(str(meanformtime[i]))#+ r'$\pm$'+'\n'+str(std[i]))
    print(xlabels)

    #t_stat, p_val = ttest_ind(refcoulomb4Form, pulsecoulomb4Form)
    #print('t-statistic:', t_stat)
    #print('p-value:', p_val)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 2,gridspec_kw={'width_ratios': [1, 4]}, sharey=True)
    plt.style.use(['science','nature','no-latex'])
    plt.subplots_adjust(wspace=0.05)
    #fig.suptitle('Coulombic efficiency', fontsize=14)

    # Add Capacheck 1
    #alphavalue = 0.25
    #CCrefCC=ax.boxplot(allrefcoulomb[1], positions = [0], patch_artist=True, boxprops=dict(facecolor='r', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    #pulseCC = ax.boxplot(allpulsecoulomb[1], positions = [0.5], patch_artist=True, boxprops=dict(facecolor='b', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #CCrefC83CC=ax.boxplot(allrefcoulombC8[1], positions = [1], patch_artist=True, boxprops=dict(facecolor='darkred', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulseC5sCC = ax.boxplot(allpulsecoulombC5s[1], positions = [1.5], patch_artist=True, boxprops=dict(facecolor='goldenrod', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulseC57CC20CC = ax.boxplot(allpulsecoulombC57CC20[1], positions = [2.0], patch_artist=True, boxprops=dict(facecolor='dimgrey', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulsedoublepulsCC = ax.boxplot(allpulsecoulombdoublepulse[1], positions = [2.5], patch_artist=True, boxprops=dict(facecolor='darkgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))
    #pulseC5CC = ax.boxplot(allpulsecoulombC5[1], positions = [3], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(color = 'k',alpha = alphavalue))

    x = [0]
    y = [0]
    for o in range(len(x)):
        CCref=ax[0].boxplot(allrefcoulomb[o], positions = [24], patch_artist=True, boxprops=dict(facecolor='r'))

        pulse = ax[1].boxplot(allpulsecoulomb[o], positions = [9.7], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
        pulsedoublepuls = ax[1].boxplot(allpulsecoulombdoublepulse[o], positions = [9.4], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen'))
        CCrefC83=ax[1].boxplot(allrefcoulombC8[o], positions = [9.2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))
    
        pulseC5 = ax[1].boxplot(allpulsecoulombC5[o], positions = [8.4], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
        pulseC57CC20 = ax[1].boxplot(allpulsecoulombC57CC20[o], positions = [8.2], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
        pulseC5s = ax[1].boxplot(allpulsecoulombC5s[o], positions = [8.1], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))
        
         
        
    legendcolors = [CCref["boxes"][0], pulse["boxes"][0], pulsedoublepuls["boxes"][0],  CCrefC83["boxes"][0], pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    #legendtext = ['N(CC20) = ' + str(len(refcoulomb4Form)),  'N(PulsC5.7) = ' + str(len(pulsecoulomb4Form)),'N(CC8.23) = ' + str(len(refcoulomb4FormC8)),
    #                   'N(PulsC5.7 5s) = ' + str(len(pulsecoulomb4FormC5s)), 'N(PulsC5.7 CC20) = ' + str(len(pulsecoulomb4FormC57CC20)),
    #                    'N(PulsC5.7double) = ' + str(len(pulsecoulomb4Formdoublepulse)), 'N(PulsC5) = ' + str(len(pulsecoulomb4FormC5))]
    legendtext = ['Reference CC C/20','Pulsed C/5.7' , 'Pulsed double length',  'CC C/8.2', 'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s']
    legend = plt.legend(legendcolors, legendtext, fontsize=4, bbox_to_anchor=(-0.25, 0.3), loc="upper left",
                 borderaxespad=0, fancybox=True, shadow=True,frameon=True, ncol = 2)
    #ax.set_ylim([0.76,0.86])

    #ax.set_xticks([0, 1, 1.5, 2,  3, 3.5, 4]) 
    ax[0].set_xticklabels([xlabels[0]])
    ax[1].set_xticklabels(xlabels[1:])

    d = .02  # Größe der diagonalen Linien in Bruchteilen der Achsengröße
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)        # Diagonale Bruchlinie rechts
    #kwargs = dict(transform=ax[1].transAxes, color='k', clip_on=False)
    ax[1].plot((1.1-d, 1.1+d), (-d, +d), **kwargs)          # Diagonale Bruchlinie links

    #ax[0].spines['bottom'].set_visible(False)
    #ax[1].spines['bottom'].set_visible(False)
    ax[1].invert_xaxis()
    ax[1].yaxis.set_visible(False)

    ax[0].tick_params(axis='x', labelsize=8, rotation = 45)
    ax[1].tick_params(axis='x', labelsize=8, rotation = 45)

    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[1].spines['top'].set_visible(False)

    ax[0].xaxis.set_ticks_position('bottom')
    ax[1].xaxis.set_ticks_position('bottom')

    ax[0].minorticks_off()
    ax[1].minorticks_off()

        
    ax[0].set_ylabel('Coulombic efficiency [%]', fontsize=8)
    #ax[1].set_xlabel('Mean charge formation time [h]', fontsize=8)
    fig.text(0.5, -0.04, 'Mean charge formation time [h]', ha='center', fontsize=8)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('CoulombEfficienymorePulsenFormation5cellsBrokenAxis.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def getEOLcoulombeffPlotFormBrokenAxisAbsoluteValues():

    # idee auf der X-achse die mediane Formierzeit anzeigen und entsprechen ordnen

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)  
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    if capacities[5] == 0:
                        del capacities[5]
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100) 
                        formtimelist.append(data['split']["Step_1"]["t"])  
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa


    refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtime, eolcapalist = getCCref()
    allrefcoulomb = [refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3]
    allrefformtime = []
    for i in range(len(formtime)):
        allrefformtime.append(formtime[i][-1])
    formchargeCC = []
    formdischargeCC = []
    for i in range(len(eolcapalist)):
        formdischargeCC.append(abs(eolcapalist[i][2]+eolcapalist[i][3]))
        formchargeCC.append(eolcapalist[i][0])

    refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8, formtimeC823, eolcapalistC823 = getCCrefC83()
    allrefcoulombC8 = [refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8]
    allrefformtimeC8 = []
    for i in range(len(formtimeC823)):
        allrefformtimeC8.append(formtimeC823[i][-1])
    formchargeCC8 = []
    formdischargeCC8 = []
    for i in range(len(eolcapalistC823)):
        formdischargeCC8.append(abs(eolcapalistC823[i][2]+eolcapalistC823[i][3]))
        formchargeCC8.append(eolcapalistC823[i][0])

    pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimepulse, eolcapalistpulse = getpulseform()
    allpulsecoulomb = [pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3]
    allpulseformtime = []
    for i in range(len(formtimepulse)):
        allpulseformtime.append(formtimepulse[i][-1])
    formchargepulse = []
    formdischargepulse = []
    for i in range(len(eolcapalistpulse)):
        formdischargepulse.append(abs(eolcapalistpulse[i][2]+eolcapalistpulse[i][3]))
        formchargepulse.append(eolcapalistpulse[i][0])


    pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5, formtimepulseC5, eolcapalistC5 = getpulseformC5()
    allpulsecoulombC5 = [pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5]
    allpulseformtimeC5 = []
    for i in range(len(formtimepulseC5)):
        allpulseformtimeC5.append(formtimepulseC5[i][-1])
    formchargepulseC5 = []
    formdischargepulseC5 = []
    for i in range(len(eolcapalistC5)):
        formdischargepulseC5.append(abs(eolcapalistC5[i][2]+eolcapalistC5[i][3]))
        formchargepulseC5.append(eolcapalistC5[i][0])

    pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s, formtimepulseC575s, eolcapalistC575s = getpulseformC575s()
    allpulsecoulombC5s = [pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s]
    allpulseformtimeC575s = [] 
    for i in range(len(formtimepulseC575s)):
        allpulseformtimeC575s.append(formtimepulseC575s[i][-1])
    formchargepulseC5s = []
    formdischargepulseC5s = []
    for i in range(len(eolcapalistC575s)):
        formdischargepulseC5s.append(abs(eolcapalistC575s[i][2]+eolcapalistC575s[i][3]))
        formchargepulseC5s.append(eolcapalistC575s[i][0])

    pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse, formtimedoublepulse, eolcapalistdoubleulse = getpulseformdoublepulse()
    allpulsecoulombdoublepulse = [pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse]
    allpulseformtimedoublepulse = []
    for i in range(len(formtimedoublepulse)):
        allpulseformtimedoublepulse.append(formtimedoublepulse[i][-1])
    formchargepulsedouble = []
    formdischargepulsedouble = []
    for i in range(len(eolcapalistdoubleulse)):
        formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][2]+eolcapalistdoubleulse[i][3]))
        formchargepulsedouble.append(eolcapalistdoubleulse[i][0])

    pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20, formtimepulseC57CC20, eolcapalistC57CC20 = getpulseformC57CC20()
    allpulsecoulombC57CC20 = [pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20]
    allpulseformtimeC57CC20 = []
    for i in range(len(formtimepulseC57CC20)):
        allpulseformtimeC57CC20.append(formtimepulseC57CC20[i][-1])
    formchargepulseC57CC20 = []
    formdischargepulseC57CC20 = []
    for i in range(len(eolcapalistC57CC20)):
        formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][2]+eolcapalistC57CC20[i][3]))
        formchargepulseC57CC20.append(eolcapalistC57CC20[i][0])

    print(allrefformtime)


    meanformtime = []
    std = []

    decround = 1

    meanformtime.append(np.round(np.mean(allrefformtime)/3600, decround))

    meanformtime.append(np.round(np.mean(allpulseformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimedoublepulse)/3600, decround))
    meanformtime.append(np.round(np.mean(allrefformtimeC8)/3600, decround))
    
    meanformtime.append(np.round(np.mean(allpulseformtimeC5)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC57CC20)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC575s)/3600, decround))
    
    std.append(np.round(np.std(allrefformtime)/3600, decround))

    std.append(np.round(np.std(allpulseformtimedoublepulse)/3600, decround))
    std.append(np.round(np.std(allpulseformtime)/3600, decround))
    std.append(np.round(np.std(allrefformtimeC8)/3600, decround))

    std.append(np.round(np.std(allpulseformtimeC5)/3600, decround))
    std.append(np.round(np.std(allpulseformtimeC57CC20)/3600, decround))
    std.append(np.round(np.std(allpulseformtimeC575s)/3600, decround))

    print(meanformtime)
    print(std)
   
    #t_stat, p_val = ttest_ind(refcoulomb4Form, pulsecoulomb4Form)
    #print('t-statistic:', t_stat)
    #print('p-value:', p_val)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(2, 2,gridspec_kw={'width_ratios': [1, 4]})
    plt.style.use(['science','nature','no-latex'])
    plt.subplots_adjust(hspace=0.4, wspace=0.04)

    xlabels = []
    for i in range(len(meanformtime)):
        xlabels.append(str(meanformtime[i]))#+ r'$\pm$'+'\n'+str(std[i]))

    linew = 0.75
    capsize = 2.5

    ax[1,0].set_xlim([23.5, 25.2])
    xerrs = np.atleast_2d([(24-np.round(min(allrefformtime)/3600, decround),np.round(max(allrefformtime)/3600, decround)-24)]).T
    ax[1,0].errorbar([24], [84.59099211865941], xerr=xerrs, lw= linew, capsize= capsize, color = 'r')

    ax[1,1].set_xlim([7, 10.2])
    xerrs = np.atleast_2d([(9.7-np.round(min(allpulseformtime)/3600, decround),np.round(max(allpulseformtime)/3600, decround)-9.7)]).T
    ax[1,1].errorbar([9.7], [83.1725850122729], xerr=xerrs, lw= linew, capsize= capsize, color = 'cornflowerblue')

    xerrs = np.atleast_2d([(9.4-np.round(min(allpulseformtimedoublepulse)/3600, decround),np.round(max(allpulseformtimedoublepulse)/3600, decround)-9.4)]).T
    ax[1,1].errorbar([9.4], [82.75780594080686], xerr=xerrs, lw= linew, capsize= capsize, color = 'cadetblue')

    xerrs = np.atleast_2d([(9.2-np.round(min(allrefformtimeC8)/3600, decround),np.round(max(allrefformtimeC8)/3600, decround)-9.2)]).T
    ax[1,1].errorbar([9.2], [81.88764905897239], xerr=xerrs, lw= linew, capsize= capsize, color = 'steelblue')

    xerrs = np.atleast_2d([(8.4-np.round(min(allpulseformtimeC5)/3600, decround),np.round(max(allpulseformtimeC5)/3600, decround)-8.4)]).T
    ax[1,1].errorbar([8.4], [81.4205783415679], xerr=xerrs, lw= linew, capsize= capsize, color = 'forestgreen')

    xerrs = np.atleast_2d([(8.2-np.round(min(allpulseformtimeC57CC20)/3600, decround),np.round(max(allpulseformtimeC57CC20)/3600, decround)-8.2)]).T
    ax[1,1].errorbar([8.2], [82.55144887217685], xerr=xerrs, lw= linew, capsize= capsize, color = 'darkgreen')

    xerrs = np.atleast_2d([(8.1-np.round(min(allpulseformtimeC575s)/3600, decround),np.round(max(allpulseformtimeC575s)/3600, decround)-8.1)]).T
    ax[1,1].errorbar([8.1], [82.01773336430657], xerr=xerrs, lw= linew, capsize= capsize, color = 'darkolivegreen')    

    chargecaparef = ax[0,0].boxplot(formchargeCC, positions = [24],  patch_artist=True, boxprops=dict(facecolor='r'))

    chargecapapulse = ax[0,1].boxplot(formchargepulse, positions = [9.7], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
    chargecapapulsedoublepuls = ax[0,1].boxplot(formchargepulsedouble, positions = [9.4], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen'))
    chargecapaCCrefC83=ax[0,1].boxplot(formchargeCC8, positions = [9.2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))

    chargecapapulseC5 = ax[0,1].boxplot(formchargepulseC5, positions = [8.4], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
    chargecapapulseC57CC20 = ax[0,1].boxplot(formchargepulseC57CC20, positions = [8.2], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
    chargecapapulseC5s = ax[0,1].boxplot(formchargepulseC5s, positions = [8.1], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))

    alphavalue = 0.3
    CCrefdc=ax[0,0].boxplot(formdischargeCC, positions = [24], patch_artist=True, boxprops=dict(facecolor='r', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulsedc = ax[0,1].boxplot(formdischargepulse, positions = [9.7], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulsedoublepulsdc = ax[0,1].boxplot(formdischargepulsedouble, positions = [9.4], patch_artist=True, boxprops=dict(facecolor='cadetblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    CCrefC83dc=ax[0,1].boxplot(formdischargeCC8, positions = [9.2], patch_artist=True, boxprops=dict(facecolor='steelblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulseC5dc = ax[0,1].boxplot(formdischargepulseC5, positions = [8.4], patch_artist=True, boxprops=dict(facecolor='forestgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC57CC20dc = ax[0,1].boxplot(formdischargepulseC57CC20, positions = [8.2], patch_artist=True, boxprops=dict(facecolor='darkgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC5sdc = ax[0,1].boxplot(formdischargepulseC5s, positions = [8.1], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    

    widthref = 0.5
    CCref=ax[1,0].boxplot(allrefcoulomb[0], positions = [24], widths=widthref, patch_artist=True, boxprops=dict(facecolor='r'))

    pulse = ax[1,1].boxplot(allpulsecoulomb[0], positions = [9.7], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
    pulsedoublepuls = ax[1,1].boxplot(allpulsecoulombdoublepulse[0], positions = [9.4], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen'))
    CCrefC83=ax[1,1].boxplot(allrefcoulombC8[0], positions = [9.2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))

    pulseC5 = ax[1,1].boxplot(allpulsecoulombC5[0], positions = [8.4], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
    pulseC57CC20 = ax[1,1].boxplot(allpulsecoulombC57CC20[0], positions = [8.2], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
    pulseC5s = ax[1,1].boxplot(allpulsecoulombC5s[0], positions = [8.1], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))

    # calculate box middle for errobar
    box_path = pulseC5s['boxes'][0].get_path()
    vertices = box_path.vertices
    Q1 = vertices[1, 1]  # y-coordinate of lower box edge
    Q3 = vertices[2, 1]  # y-coordinate of upper box edge
    boxmiddle = Q1 + (Q3-Q1)/2
    print(boxmiddle)

        
    legendcolors = [CCref["boxes"][0], pulse["boxes"][0],pulsedoublepuls["boxes"][0],   CCrefC83["boxes"][0], pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    #legendtext = ['N(CC20) = ' + str(len(refcoulomb4Form)),  'N(PulsC5.7) = ' + str(len(pulsecoulomb4Form)),'N(CC8.23) = ' + str(len(refcoulomb4FormC8)),
    #                   'N(PulsC5.7 5s) = ' + str(len(pulsecoulomb4FormC5s)), 'N(PulsC5.7 CC20) = ' + str(len(pulsecoulomb4FormC57CC20)),
    #                    'N(PulsC5.7double) = ' + str(len(pulsecoulomb4Formdoublepulse)), 'N(PulsC5) = ' + str(len(pulsecoulomb4FormC5))]
    legendtext = ['Reference CC C/20','Pulsed C/5.7' , 'Pulsed double length',  'CC C/8.2', 'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s']
    #legend = plt.legend(legendcolors, legendtext, fontsize=3, bbox_to_anchor=(-0.2, 0.45), loc="upper left",
    #             borderaxespad=0, fancybox=True, shadow=True,frameon=True, ncol = 2)
    legend = plt.legend(legendcolors, legendtext, fontsize=4, bbox_to_anchor=(-0.2, 1.37), loc="upper left",
                 borderaxespad=0, ncol = 3)

    ax[1,0].set_xticks([23.5,24,25]) 
    ax[1,0].set_xticklabels(["23.5", xlabels[0], "25.0"])
    xlabels.append("7.0")
    xlabels.insert(1, "10.2")
    ax[1,1].set_xticks([10.2, 9.7, 9.4, 9.2, 8.4, 8.2, 8.1, 7.0])
    ax[1,1].set_xticklabels(xlabels[1:])
    ax[0,0].set_xticklabels([])
    ax[0,1].set_xticklabels([])

    ax[0,0].text(-0.1, 1.3, 'a)', transform=ax[0,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
    ax[1,0].text(-0.1, 1.3, 'b)', transform=ax[1,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

    d = .02  # Größe der diagonalen Linien in Bruchteilen der Achsengröße
    kwargs = dict(transform=ax[1,0].transAxes, color='k', clip_on=False)
    ax[1,0].plot((1-d, 1+d), (-d, +d), **kwargs)        # Diagonale Bruchlinie rechts
    #kwargs = dict(transform=ax[1].transAxes, color='k', clip_on=False)
    ax[1,1].plot((1.1-d, 1.1+d), (-d, +d), **kwargs)          # Diagonale Bruchlinie links

    d = .02  # Größe der diagonalen Linien in Bruchteilen der Achsengröße
    kwargs = dict(transform=ax[0,0].transAxes, color='k', clip_on=False)
    ax[0,0].plot((1-d, 1+d), (-d, +d), **kwargs)        # Diagonale Bruchlinie rechts
    #kwargs = dict(transform=ax[1].transAxes, color='k', clip_on=False)
    ax[0,1].plot((1.1-d, 1.1+d), (-d, +d), **kwargs)          # Diagonale Bruchlinie links

    #ax[0].spines['bottom'].set_visible(False)
    #ax[1].spines['bottom'].set_visible(False)
    ax[1,1].invert_xaxis()
    ax[0,1].invert_xaxis()
    ax[1,1].yaxis.set_visible(False)
    ax[0,1].yaxis.set_visible(False)
    

    ax[1,0].tick_params(axis='x', labelsize=8, rotation = 45)
    ax[1,1].tick_params(axis='x', labelsize=8, rotation = 45)

    ax[1,0].spines['right'].set_visible(False)
    ax[1,1].spines['left'].set_visible(False)
    
    ax[0,0].spines['right'].set_visible(False)
    ax[0,1].spines['left'].set_visible(False)
    

    ax[1,0].xaxis.set_ticks_position('bottom')
    ax[1,1].xaxis.set_ticks_position('bottom')

    ax[0,0].minorticks_off()
    ax[0,1].minorticks_off()
    ax[1,0].minorticks_off()
    ax[1,1].minorticks_off()

    ax[0,0].set_ylim([2.4, 4.1])
    ax[0,1].set_ylim([2.4, 4.1])
    ax[1,0].set_ylim([76, 86])
    ax[1,1].set_ylim([76, 86])


    ax[1,0].set_ylabel('coulombic efficiency [%]', fontsize=8)
    ax[0,0].set_ylabel('capacity [mAh]', fontsize=8)
    #ax[1].set_xlabel('Mean charge formation time [h]', fontsize=8)
    fig.text(0.5, -0.04, 'mean wall time for formation charge [h]', ha='center', fontsize=8)

    fig.text(0.6, 0.8, 'charge', ha='center', fontsize=4)
    fig.text(0.6, 0.675, 'discharge', ha='center', fontsize=4)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('CoulombEfficienymorePulsenFormation5cellsBrokenAxisAbsolute.png', dpi=800)
    os.chdir(origpath)
    plt.show()



def getEOLFormChargeDischarge():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)  
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    if capacities[5] == 0:
                        del capacities[5]
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100) 
                        formtimelist.append(data['split']["Step_1"]["t"])  
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtime, eolcapalist = getCCref()
    refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8, formtimeC823, eolcapalistC823 = getCCrefC83()
    pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimepulse, eolcapalistpulse = getpulseform()
    pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5, formtimepulseC5, eolcapalistC5 = getpulseformC5()
    pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s, formtimepulseC575s, eolcapalistC575s = getpulseformC575s()
    pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse, formtimedoublepulse, eolcapalistdoubleulse = getpulseformdoublepulse()
    pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20, formtimepulseC57CC20, eolcapalistC57CC20 = getpulseformC57CC20()
    formchargeCC = []
    formdischargeCC = []
    for i in range(len(eolcapalist)):
        formdischargeCC.append(abs(eolcapalist[i][2]+eolcapalist[i][3]))
        formchargeCC.append(eolcapalist[i][0])
    allrefformtime = []
    for i in range(len(formtime)):
        allrefformtime.append(formtime[i][-1])

    
    formchargeCC8 = []
    formdischargeCC8 = []
    for i in range(len(eolcapalistC823)):
        formdischargeCC8.append(abs(eolcapalistC823[i][2]+eolcapalistC823[i][3]))
        formchargeCC8.append(eolcapalistC823[i][0])
    allrefformtimeC8 = []
    for i in range(len(formtimeC823)):
        allrefformtimeC8.append(formtimeC823[i][-1])

    
    formchargepulse = []
    formdischargepulse = []
    for i in range(len(eolcapalistpulse)):
        formdischargepulse.append(abs(eolcapalistpulse[i][2]+eolcapalistpulse[i][3]))
        formchargepulse.append(eolcapalistpulse[i][0])
    allpulseformtime = []
    for i in range(len(formtimepulse)):
        allpulseformtime.append(formtimepulse[i][-1])

    
    formchargepulseC5 = []
    formdischargepulseC5 = []
    for i in range(len(eolcapalistC5)):
        formdischargepulseC5.append(abs(eolcapalistC5[i][2]+eolcapalistC5[i][3]))
        formchargepulseC5.append(eolcapalistC5[i][0])
    allpulseformtimeC5 = []
    for i in range(len(formtimepulseC5)):
        allpulseformtimeC5.append(formtimepulseC5[i][-1])

    
    formchargepulseC5s = []
    formdischargepulseC5s = []
    for i in range(len(eolcapalistC575s)):
        formdischargepulseC5s.append(abs(eolcapalistC575s[i][2]+eolcapalistC575s[i][3]))
        formchargepulseC5s.append(eolcapalistC575s[i][0])
    allpulseformtimeC575s = [] 
    for i in range(len(formtimepulseC575s)):
        allpulseformtimeC575s.append(formtimepulseC575s[i][-1])

    
    formchargepulsedouble = []
    formdischargepulsedouble = []
    for i in range(len(eolcapalistdoubleulse)):
        formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][2]+eolcapalistdoubleulse[i][3]))
        formchargepulsedouble.append(eolcapalistdoubleulse[i][0])
    allpulseformtimedoublepulse = []
    for i in range(len(formtimedoublepulse)):
        allpulseformtimedoublepulse.append(formtimedoublepulse[i][-1])

    
    formchargepulseC57CC20 = []
    formdischargepulseC57CC20 = []
    for i in range(len(eolcapalistC57CC20)):
        formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][2]+eolcapalistC57CC20[i][3]))
        formchargepulseC57CC20.append(eolcapalistC57CC20[i][0])
    allpulseformtimeC57CC20 = []
    for i in range(len(formtimepulseC57CC20)):
        allpulseformtimeC57CC20.append(formtimepulseC57CC20[i][-1])

    print(allrefformtime)
    meanformtime = []
    decround = 1
    meanformtime.append(np.round(np.mean(allrefformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimedoublepulse)/3600, decround))
    meanformtime.append(np.round(np.mean(allrefformtimeC8)/3600, decround))
    
    meanformtime.append(np.round(np.mean(allpulseformtimeC5)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC57CC20)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC575s)/3600, decround))
    

    print(meanformtime)

    xlabels = []
    for i in range(len(meanformtime)):
        xlabels.append(str(meanformtime[i]))
    print(xlabels)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1)#, sharex=True, sharey=True)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)


    CCref=ax.boxplot(formchargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r'), medianprops=dict(color = 'limegreen'))

    pulse = ax.boxplot(formchargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
    pulsedoublepuls = ax.boxplot(formchargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen'))
    CCrefC83=ax.boxplot(formchargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))
    
    pulseC5 = ax.boxplot(formchargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
    pulseC57CC20 = ax.boxplot(formchargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
    pulseC5s = ax.boxplot(formchargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))
   
    
    alphavalue = 0.3
    CCrefdc=ax.boxplot(formdischargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    
    pulsedc = ax.boxplot(formdischargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulsedoublepulsdc = ax.boxplot(formdischargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    CCrefC83dc=ax.boxplot(formdischargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulseC5dc = ax.boxplot(formdischargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC57CC20dc = ax.boxplot(formdischargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC5sdc = ax.boxplot(formdischargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
                
        
    legendcolors = [CCref["boxes"][0],pulse["boxes"][0], pulsedoublepuls["boxes"][0],  CCrefC83["boxes"][0], pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    #legendtext = ['N(CC20) = ' + str(len(refcoulomb4Form)),  'N(PulsC5.7) = ' + str(len(pulsecoulomb4Form)),'N(CC8.23) = ' + str(len(refcoulomb4FormC8)),
    #                   'N(PulsC5.7 5s) = ' + str(len(pulsecoulomb4FormC5s)), 'N(PulsC5.7 CC20) = ' + str(len(pulsecoulomb4FormC57CC20)),
    #                    'N(PulsC5.7double) = ' + str(len(pulsecoulomb4Formdoublepulse)), 'N(PulsC5) = ' + str(len(pulsecoulomb4FormC5))]
    legendtext = ['Reference CC C/20', 'Pulsed C/5.7' ,'Pulsed double length', 'CC C/8.2', 'Pulsed C/5',  'Mixed C/5.7+C/20', 'Pulsed 5s' ]
    legend = plt.legend(legendcolors, legendtext, bbox_to_anchor=(1.01, 0.99),
                 borderaxespad=0, ncol = 1)

    ax.set_xlim([-0.25,4.25])

    ax.set_xticks([0, 1, 1.5, 2,  3, 3.5, 4]) 
    ax.set_xticklabels(xlabels)

    fig.text(0.6, 0.675, 'charge', ha='center', fontsize=6)
    fig.text(0.6, 0.4, 'discharge', ha='center', fontsize=6)        

    ax.set_xlabel('mean wall time for formation charge [h]', fontsize=8)
    ax.set_ylabel('capacity [mAh]', fontsize=8)
    
    origpath = os.getcwd()
    os.chdir(plotfolder)
    fig.savefig('FormationChargeDischarge5cells.png', dpi=800 )
    os.chdir(origpath)
    plt.show()

def getEOLCapaCheck1ChargeDischarge():


    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)  
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    if capacities[5] == 0:
                        del capacities[5]
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100) 
                        formtimelist.append(data['split']["Step_1"]["t"])  
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa


    refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtime, eolcapalist = getCCref()
    refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8, formtimeC823, eolcapalistC823 = getCCrefC83()
    pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimepulse, eolcapalistpulse = getpulseform()
    pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5, formtimepulseC5, eolcapalistC5 = getpulseformC5()
    pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s, formtimepulseC575s, eolcapalistC575s = getpulseformC575s()
    pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse, formtimedoublepulse, eolcapalistdoubleulse = getpulseformdoublepulse()
    pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20, formtimepulseC57CC20, eolcapalistC57CC20 = getpulseformC57CC20()


    formchargeCC = []
    formdischargeCC = []
    for i in range(len(eolcapalist)):
        if len(eolcapalist[i]) == 37:
            formdischargeCC.append(abs(eolcapalist[i][8]))
            formchargeCC.append(eolcapalist[i][5]+ eolcapalist[i][6])
        else:
            formdischargeCC.append(abs(eolcapalist[i][9]))
            formchargeCC.append(eolcapalist[i][6]+ eolcapalist[i][7])
    allrefformtime = []
    for i in range(len(formtime)):
        allrefformtime.append(formtime[i][-1])

    for i in range(len(eolcapalist)):
        if len(eolcapalist[i]) == 38:
            print(str(eolcapalist[i][6])+ '  '+str( eolcapalist[i][7])) 
        else:
            print(str(eolcapalist[i][5])+ '  '+str( eolcapalist[i][6])) 
            
    formchargeCC8 = []
    formdischargeCC8 = []
    for i in range(len(eolcapalistC823)):
        if len(eolcapalistC823[i]) == 37:
            formdischargeCC8.append(abs(eolcapalistC823[i][8]))
            formchargeCC8.append(eolcapalistC823[i][5]+ eolcapalistC823[i][6])
        else:
            formdischargeCC8.append(abs(eolcapalistC823[i][9]))
            formchargeCC8.append(eolcapalistC823[i][6]+ eolcapalistC823[i][7])
    allrefformtimeC8 = []
    for i in range(len(formtimeC823)):
        allrefformtimeC8.append(formtimeC823[i][-1])

    formchargepulse = []
    formdischargepulse = []
    for i in range(len(eolcapalistpulse)):
        if len(eolcapalistpulse[i]) == 37:
            formdischargepulse.append(abs(eolcapalistpulse[i][8]))
            formchargepulse.append(eolcapalistpulse[i][5]+ eolcapalistpulse[i][6])
        else:
            formdischargepulse.append(abs(eolcapalistpulse[i][9]))
            formchargepulse.append(eolcapalistpulse[i][6]+ eolcapalistpulse[i][7])
    allpulseformtime = []
    for i in range(len(formtimepulse)):
        allpulseformtime.append(formtimepulse[i][-1])

    formchargepulseC5 = []
    formdischargepulseC5 = []
    for i in range(len(eolcapalistC5)):
        if len(eolcapalistC5[i]) == 37:
            formdischargepulseC5.append(abs(eolcapalistC5[i][8]))
            formchargepulseC5.append(eolcapalistC5[i][5]+ eolcapalistC5[i][6])
        else:
            formdischargepulseC5.append(abs(eolcapalistC5[i][9]))
            formchargepulseC5.append(eolcapalistC5[i][6]+ eolcapalistC5[i][7])
    allpulseformtimeC5 = []
    for i in range(len(formtimepulseC5)):
        allpulseformtimeC5.append(formtimepulseC5[i][-1])

    formchargepulseC5s = []
    formdischargepulseC5s = []
    for i in range(len(eolcapalistC575s)):
        if len(eolcapalistC575s[i]) == 37:
            formdischargepulseC5s.append(abs(eolcapalistC575s[i][8]))
            formchargepulseC5s.append(eolcapalistC575s[i][5]+ eolcapalistC575s[i][6])
        else:
            formdischargepulseC5s.append(abs(eolcapalistC575s[i][9]))
            formchargepulseC5s.append(eolcapalistC575s[i][6]+ eolcapalistC575s[i][7])
    allpulseformtimeC575s = [] 
    for i in range(len(formtimepulseC575s)):
        allpulseformtimeC575s.append(formtimepulseC575s[i][-1])


    formchargepulsedouble = []
    formdischargepulsedouble = []
    for i in range(len(eolcapalistdoubleulse)):
        if len(eolcapalistdoubleulse[i]) == 37:
            formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][8]))
            formchargepulsedouble.append(eolcapalistdoubleulse[i][5]+ eolcapalistdoubleulse[i][6])
        else:
            formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][9]))
            formchargepulsedouble.append(eolcapalistdoubleulse[i][6]+ eolcapalistdoubleulse[i][7])
    allpulseformtimedoublepulse = []
    for i in range(len(formtimedoublepulse)):
        allpulseformtimedoublepulse.append(formtimedoublepulse[i][-1])


    formchargepulseC57CC20 = []
    formdischargepulseC57CC20 = []
    for i in range(len(eolcapalistC57CC20)):
        if len(eolcapalistC57CC20[i]) == 37:
            formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][8]))
            formchargepulseC57CC20.append(eolcapalistC57CC20[i][5]+ eolcapalistC57CC20[i][6])
        else:
            formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][9]))
            formchargepulseC57CC20.append(eolcapalistC57CC20[i][6]+ eolcapalistC57CC20[i][7])
    allpulseformtimeC57CC20 = []
    for i in range(len(formtimepulseC57CC20)):
        allpulseformtimeC57CC20.append(formtimepulseC57CC20[i][-1])

    print(allrefformtime)


    meanformtime = []

    decround = 1

    meanformtime.append(np.round(np.mean(allrefformtime)/3600, decround))

    meanformtime.append(np.round(np.mean(allpulseformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimedoublepulse)/3600, decround))
    meanformtime.append(np.round(np.mean(allrefformtimeC8)/3600, decround))
    
    meanformtime.append(np.round(np.mean(allpulseformtimeC5)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC57CC20)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC575s)/3600, decround))
    

    print(meanformtime)

    xlabels = []
    for i in range(len(meanformtime)):
        xlabels.append(str(meanformtime[i]))
    print(xlabels)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1)#, sharex=True, sharey=True)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)


    CCref=ax.boxplot(formchargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r'), medianprops=dict(color = 'limegreen'))

    pulse = ax.boxplot(formchargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
    pulsedoublepuls = ax.boxplot(formchargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen'))
    CCrefC83=ax.boxplot(formchargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))
    
    pulseC5 = ax.boxplot(formchargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
    pulseC57CC20 = ax.boxplot(formchargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
    pulseC5s = ax.boxplot(formchargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))
       
    alphavalue = 0.3
    CCrefdc=ax.boxplot(formdischargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulsedc = ax.boxplot(formdischargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue)) 
    pulsedoublepulsdc = ax.boxplot(formdischargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))  
    CCrefC83dc=ax.boxplot(formdischargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulseC5dc = ax.boxplot(formdischargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC57CC20dc = ax.boxplot(formdischargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC5sdc = ax.boxplot(formdischargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
                
    legendcolors = [CCref["boxes"][0],pulse["boxes"][0], pulsedoublepuls["boxes"][0],  CCrefC83["boxes"][0], pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    #legendtext = ['N(CC20) = ' + str(len(refcoulomb4Form)),  'N(PulsC5.7) = ' + str(len(pulsecoulomb4Form)),'N(CC8.23) = ' + str(len(refcoulomb4FormC8)),
    #                   'N(PulsC5.7 5s) = ' + str(len(pulsecoulomb4FormC5s)), 'N(PulsC5.7 CC20) = ' + str(len(pulsecoulomb4FormC57CC20)),
    #                    'N(PulsC5.7double) = ' + str(len(pulsecoulomb4Formdoublepulse)), 'N(PulsC5) = ' + str(len(pulsecoulomb4FormC5))]
    legendtext = ['Reference CC C/20','Pulsed C/5.7' ,  'Pulsed double length','CC C/8.2', 'Pulsed C/5',  'Mixed C/5.7+C/20', 'Pulsed 5s' ]
    legend = plt.legend(legendcolors, legendtext, bbox_to_anchor=(1.01, 0.99),
                 borderaxespad=0, ncol = 1)

    ax.set_xlim([-0.25,4.25])

    ax.set_xticks([0, 1, 1.5, 2,  3, 3.5, 4]) 
    ax.set_xticklabels(xlabels)

    fig.text(0.6, 0.4, 'charge', ha='center', fontsize=6)
    fig.text(0.6, 0.25, 'discharge', ha='center', fontsize=6)        

    ax.set_xlabel('mean wall time for formation charge [h]', fontsize=8)
    ax.set_ylabel('capacity [mAh]', fontsize=8)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    fig.savefig('CapaCheck1ChargeDischarge5cells.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def getEOLCapaCheck2ChargeDischarge():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)  
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    if capacities[5] == 0:
                        del capacities[5]
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100) 
                        formtimelist.append(data['split']["Step_1"]["t"])  
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa


    refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtime, eolcapalist = getCCref()
    refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8, formtimeC823, eolcapalistC823 = getCCrefC83()
    pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimepulse, eolcapalistpulse = getpulseform()
    pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5, formtimepulseC5, eolcapalistC5 = getpulseformC5()
    pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s, formtimepulseC575s, eolcapalistC575s = getpulseformC575s()
    pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse, formtimedoublepulse, eolcapalistdoubleulse = getpulseformdoublepulse()
    pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20, formtimepulseC57CC20, eolcapalistC57CC20 = getpulseformC57CC20()

    formchargeCC = []
    formdischargeCC = []
    for i in range(len(eolcapalist)):
        if len(eolcapalist[i]) == 37:
            formdischargeCC.append(abs(eolcapalist[i][27]))
            formchargeCC.append(eolcapalist[i][24]+ eolcapalist[i][25])
        else:
            formdischargeCC.append(abs(eolcapalist[i][28]))
            formchargeCC.append(eolcapalist[i][25]+ eolcapalist[i][26])
    allrefformtime = []
    for i in range(len(formtime)):
        allrefformtime.append(formtime[i][-1])


    formchargeCC8 = []
    formdischargeCC8 = []
    for i in range(len(eolcapalistC823)):
        if len(eolcapalistC823[i]) == 37:
            formdischargeCC8.append(abs(eolcapalistC823[i][27]))
            formchargeCC8.append(eolcapalistC823[i][24]+ eolcapalistC823[i][25])
        else:
            formdischargeCC8.append(abs(eolcapalistC823[i][28]))
            formchargeCC8.append(eolcapalistC823[i][25]+ eolcapalistC823[i][26])
    allrefformtimeC8 = []
    for i in range(len(formtimeC823)):
        allrefformtimeC8.append(formtimeC823[i][-1])

    formchargepulse = []
    formdischargepulse = []
    for i in range(len(eolcapalistpulse)):
        if len(eolcapalistpulse[i]) == 37:
            formdischargepulse.append(abs(eolcapalistpulse[i][27]))
            formchargepulse.append(eolcapalistpulse[i][24]+ eolcapalistpulse[i][25])
        else:
            formdischargepulse.append(abs(eolcapalistpulse[i][28]))
            formchargepulse.append(eolcapalistpulse[i][25]+ eolcapalistpulse[i][26])
    allpulseformtime = []
    for i in range(len(formtimepulse)):
        allpulseformtime.append(formtimepulse[i][-1])


    formchargepulseC5 = []
    formdischargepulseC5 = []
    for i in range(len(eolcapalistC5)):
        if len(eolcapalistC5[i]) == 37:
            formdischargepulseC5.append(abs(eolcapalistC5[i][27]))
            formchargepulseC5.append(eolcapalistC5[i][24]+ eolcapalistC5[i][25])
        else:
            formdischargepulseC5.append(abs(eolcapalistC5[i][28]))
            formchargepulseC5.append(eolcapalistC5[i][25]+ eolcapalistC5[i][26])
    allpulseformtimeC5 = []
    for i in range(len(formtimepulseC5)):
        allpulseformtimeC5.append(formtimepulseC5[i][-1])


    formchargepulseC5s = []
    formdischargepulseC5s = []
    for i in range(len(eolcapalistC575s)):
        if len(eolcapalistC575s[i]) == 37:
            formdischargepulseC5s.append(abs(eolcapalistC575s[i][27]))
            formchargepulseC5s.append(eolcapalistC575s[i][24]+ eolcapalistC575s[i][25])
        else:
            formdischargepulseC5s.append(abs(eolcapalistC575s[i][28]))
            formchargepulseC5s.append(eolcapalistC575s[i][25]+ eolcapalistC575s[i][26])
    allpulseformtimeC575s = [] 
    for i in range(len(formtimepulseC575s)):
        allpulseformtimeC575s.append(formtimepulseC575s[i][-1])


    formchargepulsedouble = []
    formdischargepulsedouble = []
    for i in range(len(eolcapalistdoubleulse)):
        if len(eolcapalistdoubleulse[i]) == 37:
            formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][27]))
            formchargepulsedouble.append(eolcapalistdoubleulse[i][24]+ eolcapalistdoubleulse[i][25])
        else:
            formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][28]))
            formchargepulsedouble.append(eolcapalistdoubleulse[i][25]+ eolcapalistdoubleulse[i][26])
    allpulseformtimedoublepulse = []
    for i in range(len(formtimedoublepulse)):
        allpulseformtimedoublepulse.append(formtimedoublepulse[i][-1])


    formchargepulseC57CC20 = []
    formdischargepulseC57CC20 = []
    for i in range(len(eolcapalistC57CC20)):
        if len(eolcapalistC57CC20[i]) == 37:
            formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][27]))
            formchargepulseC57CC20.append(eolcapalistC57CC20[i][24]+ eolcapalistC57CC20[i][25])
        else:
            formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][28]))
            formchargepulseC57CC20.append(eolcapalistC57CC20[i][25]+ eolcapalistC57CC20[i][26])
    allpulseformtimeC57CC20 = []
    for i in range(len(formtimepulseC57CC20)):
        allpulseformtimeC57CC20.append(formtimepulseC57CC20[i][-1])

    print(allrefformtime)


    meanformtime = []

    decround = 1

    meanformtime.append(np.round(np.mean(allrefformtime)/3600, decround))

    meanformtime.append(np.round(np.mean(allpulseformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimedoublepulse)/3600, decround))
    meanformtime.append(np.round(np.mean(allrefformtimeC8)/3600, decround))
    
    meanformtime.append(np.round(np.mean(allpulseformtimeC5)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC57CC20)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC575s)/3600, decround))
    

    print(meanformtime)
    xlabels = []
    for i in range(len(meanformtime)):
        xlabels.append(str(meanformtime[i]))
    print(xlabels)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1)#, sharex=True, sharey=True)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)


    CCref=ax.boxplot(formchargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r'), medianprops=dict(color = 'limegreen'))

    pulse = ax.boxplot(formchargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
    pulsedoublepuls = ax.boxplot(formchargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen'))
    CCrefC83=ax.boxplot(formchargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))
    
    pulseC5 = ax.boxplot(formchargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
    pulseC57CC20 = ax.boxplot(formchargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
    pulseC5s = ax.boxplot(formchargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))
    
    ax.set_xlim([-0.25,4.25])
        

    
    alphavalue = 0.3
    CCrefdc=ax.boxplot(formdischargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulsedc = ax.boxplot(formdischargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))   
    pulsedoublepulsdc = ax.boxplot(formdischargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    CCrefC83dc=ax.boxplot(formdischargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulseC5dc = ax.boxplot(formdischargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC57CC20dc = ax.boxplot(formdischargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC5sdc = ax.boxplot(formdischargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
                
    legendcolors = [CCref["boxes"][0], pulse["boxes"][0],pulsedoublepuls["boxes"][0],   CCrefC83["boxes"][0], pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    legendtext = ['Reference CC C/20', 'Pulsed C/5.7' ,'Pulsed double length', 'CC C/8.2', 'Pulsed C/5',  'Mixed C/5.7+C/20', 'Pulsed 5s' ]
    legend = plt.legend(legendcolors, legendtext, bbox_to_anchor=(1.01, 0.99),
                 borderaxespad=0, ncol = 1)
    
    fig.text(0.25, 0.55, 'charge', ha='center', fontsize=6)
    fig.text(0.25, 0.625, 'discharge', ha='center', fontsize=6)        

    ax.set_xlabel('mean wall time for formation charge [h]', fontsize=8)
    ax.set_ylabel('capacity [mAh]', fontsize=8)

    ax.set_xticks([0,  1, 1.5, 2, 3, 3.5, 4]) 
    ax.set_xticklabels(xlabels)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    fig.savefig('CapaCheck2ChargeDischarge5cells.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def getEOLCapaCheck3ChargeDischarge():


    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        refcoulomb4Form = []
        refcoulomb4Capa1 = []
        refcoulomb4Capa2 = []
        refcoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        refcoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        refcoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        refcoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        refcoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        refcoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtimelist, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)  
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    if capacities[5] == 0:
                        del capacities[5]
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100) 
                        formtimelist.append(data['split']["Step_1"]["t"])  
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        pulsecoulomb4Form = []
        pulsecoulomb4Capa1 = []
        pulsecoulomb4Capa2 = []
        pulsecoulomb4Capa3 = []
        formtimelist = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> -0.000164:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[32]/(capacities[29]+capacities[30])))
                        pulsecoulomb4Capa2.append(abs(capacities[27]/(capacities[24]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[8]/(capacities[5]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
                else:
                    try:
                        pulsecoulomb4Capa3.append(abs(capacities[33]/(capacities[30]+capacities[31])))
                        pulsecoulomb4Capa2.append(abs(capacities[28]/(capacities[26]+capacities[25])))
                        pulsecoulomb4Capa1.append(abs(capacities[9]/(capacities[7]+capacities[6])))
                        pulsecoulomb4Form.append(abs((capacities[2]+capacities[3]+ capacities[4])/(capacities[0]))*100)   
                        formtimelist.append(data['split']["Step_1"]["t"])
                    except:
                        print(i)
        return pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimelist, pulsecapa


    refcoulomb4Form, refcoulomb4Capa1, refcoulomb4Capa2, refcoulomb4Capa3, formtime, eolcapalist = getCCref()
    refcoulomb4FormC8, refcoulomb4Capa1C8, refcoulomb4Capa2C8, refcoulomb4Capa3C8, formtimeC823, eolcapalistC823 = getCCrefC83()
    pulsecoulomb4Form, pulsecoulomb4Capa1, pulsecoulomb4Capa2, pulsecoulomb4Capa3, formtimepulse, eolcapalistpulse = getpulseform()
    pulsecoulomb4FormC5, pulsecoulomb4Capa1C5, pulsecoulomb4Capa2C5, pulsecoulomb4Capa3C5, formtimepulseC5, eolcapalistC5 = getpulseformC5()
    pulsecoulomb4FormC5s, pulsecoulomb4Capa1C5s, pulsecoulomb4Capa2C5s, pulsecoulomb4Capa3C5s, formtimepulseC575s, eolcapalistC575s = getpulseformC575s()
    pulsecoulomb4Formdoublepulse, pulsecoulomb4Capa1doublepulse, pulsecoulomb4Capa2doublepulse, pulsecoulomb4Capa3doublepulse, formtimedoublepulse, eolcapalistdoubleulse = getpulseformdoublepulse()
    pulsecoulomb4FormC57CC20, pulsecoulomb4Capa1C57CC20, pulsecoulomb4Capa2C57CC20, pulsecoulomb4Capa3C57CC20, formtimepulseC57CC20, eolcapalistC57CC20 = getpulseformC57CC20()

    formchargeCC = []
    formdischargeCC = []
    for i in range(len(eolcapalist)):
        if len(eolcapalist[i]) == 37:
            formdischargeCC.append(abs(eolcapalist[i][32]))
            formchargeCC.append(eolcapalist[i][29]+ eolcapalist[i][30])
        else:
            formdischargeCC.append(abs(eolcapalist[i][33]))
            formchargeCC.append(eolcapalist[i][30]+ eolcapalist[i][31])
    allrefformtime = []
    for i in range(len(formtime)):
        allrefformtime.append(formtime[i][-1])


    formchargeCC8 = []
    formdischargeCC8 = []
    for i in range(len(eolcapalistC823)):
        if len(eolcapalistC823[i]) == 37:
            formdischargeCC8.append(abs(eolcapalistC823[i][32]))
            formchargeCC8.append(eolcapalistC823[i][29]+ eolcapalistC823[i][30])
        else:
            formdischargeCC8.append(abs(eolcapalistC823[i][33]))
            formchargeCC8.append(eolcapalistC823[i][30]+ eolcapalistC823[i][31])
    allrefformtimeC8 = []
    for i in range(len(formtimeC823)):
        allrefformtimeC8.append(formtimeC823[i][-1])

    formchargepulse = []
    formdischargepulse = []
    for i in range(len(eolcapalistpulse)):
        if len(eolcapalistpulse[i]) == 37:
            formdischargepulse.append(abs(eolcapalistpulse[i][32]))
            formchargepulse.append(eolcapalistpulse[i][29]+ eolcapalistpulse[i][30])
        else:
            formdischargepulse.append(abs(eolcapalistpulse[i][33]))
            formchargepulse.append(eolcapalistpulse[i][30]+ eolcapalistpulse[i][31])
    allpulseformtime = []
    for i in range(len(formtimepulse)):
        allpulseformtime.append(formtimepulse[i][-1])


    formchargepulseC5 = []
    formdischargepulseC5 = []
    for i in range(len(eolcapalistC5)):
        if len(eolcapalistC5[i]) == 37:
            formdischargepulseC5.append(abs(eolcapalistC5[i][32]))
            formchargepulseC5.append(eolcapalistC5[i][29]+ eolcapalistC5[i][30])
        else:
            formdischargepulseC5.append(abs(eolcapalistC5[i][33]))
            formchargepulseC5.append(eolcapalistC5[i][30]+ eolcapalistC5[i][31])
    allpulseformtimeC5 = []
    for i in range(len(formtimepulseC5)):
        allpulseformtimeC5.append(formtimepulseC5[i][-1])


    formchargepulseC5s = []
    formdischargepulseC5s = []
    for i in range(len(eolcapalistC575s)):
        if len(eolcapalistC575s[i]) == 37:
            formdischargepulseC5s.append(abs(eolcapalistC575s[i][32]))
            formchargepulseC5s.append(eolcapalistC575s[i][29]+ eolcapalistC575s[i][30])
        else:
            formdischargepulseC5s.append(abs(eolcapalistC575s[i][33]))
            formchargepulseC5s.append(eolcapalistC575s[i][30]+ eolcapalistC575s[i][31])
    allpulseformtimeC575s = [] 
    for i in range(len(formtimepulseC575s)):
        allpulseformtimeC575s.append(formtimepulseC575s[i][-1])


    formchargepulsedouble = []
    formdischargepulsedouble = []
    for i in range(len(eolcapalistdoubleulse)):
        if len(eolcapalistdoubleulse[i]) == 37:
            formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][32]))
            formchargepulsedouble.append(eolcapalistdoubleulse[i][29]+ eolcapalistdoubleulse[i][30])
        else:
            formdischargepulsedouble.append(abs(eolcapalistdoubleulse[i][33]))
            formchargepulsedouble.append(eolcapalistdoubleulse[i][30]+ eolcapalistdoubleulse[i][31])
    allpulseformtimedoublepulse = []
    for i in range(len(formtimedoublepulse)):
        allpulseformtimedoublepulse.append(formtimedoublepulse[i][-1])


    formchargepulseC57CC20 = []
    formdischargepulseC57CC20 = []
    for i in range(len(eolcapalistC57CC20)):
        if len(eolcapalistC57CC20[i]) == 37:
            formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][32]))
            formchargepulseC57CC20.append(eolcapalistC57CC20[i][29]+ eolcapalistC57CC20[i][30])
        else:
            formdischargepulseC57CC20.append(abs(eolcapalistC57CC20[i][33]))
            formchargepulseC57CC20.append(eolcapalistC57CC20[i][30]+ eolcapalistC57CC20[i][31])
    allpulseformtimeC57CC20 = []
    for i in range(len(formtimepulseC57CC20)):
        allpulseformtimeC57CC20.append(formtimepulseC57CC20[i][-1])

    print(allrefformtime)


    meanformtime = []

    decround = 1

    meanformtime.append(np.round(np.mean(allrefformtime)/3600, decround))

    meanformtime.append(np.round(np.mean(allpulseformtime)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimedoublepulse)/3600, decround))
    meanformtime.append(np.round(np.mean(allrefformtimeC8)/3600, decround))
    
    meanformtime.append(np.round(np.mean(allpulseformtimeC5)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC57CC20)/3600, decround))
    meanformtime.append(np.round(np.mean(allpulseformtimeC575s)/3600, decround))
    

    print(meanformtime)


    xlabels = []
    for i in range(len(meanformtime)):
        xlabels.append(str(meanformtime[i]))
    print(xlabels)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1)#, sharex=True, sharey=True)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)


    CCref=ax.boxplot(formchargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r'), medianprops=dict(color = 'limegreen'))

    pulse = ax.boxplot(formchargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue'), medianprops=dict(color = 'limegreen'))
    pulsedoublepuls = ax.boxplot(formchargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue'), medianprops=dict(color = 'limegreen'))
    CCrefC83=ax.boxplot(formchargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue'), medianprops=dict(color = 'limegreen'))
    
    pulseC5 = ax.boxplot(formchargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen'), medianprops=dict(color = 'limegreen'))
    pulseC57CC20 = ax.boxplot(formchargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen'), medianprops=dict(color = 'limegreen'))
    pulseC5s = ax.boxplot(formchargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen'), medianprops=dict(color = 'limegreen'))
    

    ax.set_xlim([-0.25,4.25])
        
    alphavalue = 0.3
    CCrefdc=ax.boxplot(formdischargeCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
    pulsedc = ax.boxplot(formdischargepulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))   
    pulsedoublepulsdc = ax.boxplot(formdischargepulsedouble, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    CCrefC83dc=ax.boxplot(formdischargeCC8, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
   
    pulseC5dc = ax.boxplot(formdischargepulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC57CC20dc = ax.boxplot(formdischargepulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    pulseC5sdc = ax.boxplot(formdischargepulseC5s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', alpha= alphavalue), whiskerprops=dict(alpha=alphavalue), medianprops=dict(color = 'k',alpha = alphavalue), capprops=dict(alpha = alphavalue))
    
                
    legendcolors = [CCref["boxes"][0],pulse["boxes"][0], pulsedoublepuls["boxes"][0],   CCrefC83["boxes"][0],pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    legendtext = ['Reference CC C/20', 'Pulsed C/5.7' ,'Pulsed double length',  'CC C/8.2', 'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s' ]
    legend = plt.legend(legendcolors, legendtext, bbox_to_anchor=(1.01, 0.99),
                 borderaxespad=0, ncol = 1)
    
    fig.text(0.25, 0.575, 'charge', ha='center', fontsize=6)
    fig.text(0.25, 0.5, 'discharge', ha='center', fontsize=6)        

    ax.set_xlabel('mean wall time for formation charge [h]', fontsize=8)
    ax.set_ylabel('capacity [mAh]', fontsize=8)

    ax.set_xticks([0,  1, 1.5, 2, 3, 3.5, 4]) 
    ax.set_xticklabels(xlabels)
    
    origpath = os.getcwd()
    os.chdir(plotfolder)
    fig.savefig('CapaCheck3ChargeDischarge5cells.png', dpi=800)
    os.chdir(origpath)
    plt.show()


def getEOLselfdischargeVoltageDecay():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        CCcapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname
    
    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname

    selfdischarge, selfdccapa, voltlistref, timelistref, currentlistC2ref, voltlistC2ref, cellnamesref = getCCref()

    selfdischargepulse, selfdccapapulse, voltlistpulse, timelistpulse, currentlistC2pulse, voltlistC2pulse, cellnamespulse = getpulseform()

    selfdischargeC82, selfdccapaC82, voltlistC82, timelistC82, currentlistCC2C82, voltlistCC2C82, cellnamesC82 = getCCrefC83()

    selfdischargepulseC5, selfdccapapulseC5, voltlistpulseC5, timelistpulseC5, currentlistC2pulseC5, voltlistC2pulseC5, cellnamespulseC5 = getpulseformC5()

    selfdischargepulseC575s, selfdccapapulseC575s, voltlistpulseC575s, timelistpulseC575s, currentlistC2pulseC575s, voltlistC2pulseC575s, cellnamespulseC575s = getpulseformC575s()

    selfdischargepulsedouble, selfdccapapulsedouble, voltlistpulsedouble, timelistpulsedouble, currentlistC2pulsedouble, voltlistC2pulsedouble, cellnamespulsedouble = getpulseformdoublepulse()

    selfdischargepulseC57CC20, selfdccapapulseC57CC20, voltlistpulseC57CC20, timelistpulseC57CC20, currentlistC2pulseC57CC20, voltlistC2pulseC57CC20, cellnamespulseC57CC20 = getpulseformC57CC20()

    def change_time_sqrt(timelist, voltlist):
        fittimepulse = []
        fitvoltpulse = []
        for p in range(len(timelist)):
            temptime = []
            tempvolt = []
            for time in range(len(timelist[p])):
                if timelist[p][time]/(60*60) > 24 :
                    temptime.append(np.sqrt(timelist[p][time]/(60*60)))
                    tempvolt.append(voltlist[p][time])
                #timelist[p][time] = timelist[p][time]/(60*60)
            fittimepulse.append(temptime)
            fitvoltpulse.append(tempvolt)
        return fittimepulse, fitvoltpulse

    fittimeref, fitvoltref = change_time_sqrt(timelistref, voltlistref)
    fittimeC82, fitvoltC82 = change_time_sqrt(timelistC82, voltlistC82)

    fittimepulse, fitvoltpulse = change_time_sqrt(timelistpulse, voltlistpulse)
    fittimepulseC5, fitvoltpulseC5 = change_time_sqrt(timelistpulseC5, voltlistpulseC5)
    fittimepulseC575s, fitvoltpulseC575s = change_time_sqrt(timelistpulseC575s, voltlistpulseC575s)
    fittimepulsedouble, fitvoltpulsedouble = change_time_sqrt(timelistpulsedouble, voltlistpulsedouble)
    fittimepulseC57CC20, fitvoltpulseC57CC20 = change_time_sqrt(timelistpulseC57CC20, voltlistpulseC57CC20)

    def linear_sqrt_function(x, a, b):
        return a * x + b
    
    def linear_fit(timelist, voltlist):
        voltdecayref = []
        rsqrt = []
        for r in range(len(timelist)):
            params, covariance = curve_fit(linear_sqrt_function, timelist[r], voltlist[r])
            a, b = params
            residuals = []
            totals = []
            for u in range(len(voltlist)):
                residuals.append((voltlist[r][u]- linear_sqrt_function(timelist[r][u], a, b))**2)
                totals.append((voltlist[r][u]-np.mean(voltlist[r]))**2)
            ss_res = np.sum(residuals)
            ss_tot = np.sum(totals)
            r_sq = np.round(1 - (ss_res / ss_tot), 5)
            voltdecayref.append(a *1000)
            rsqrt.append(r_sq)
        return voltdecayref, rsqrt

    voltdecayref, rsqrtref = linear_fit(fittimeref, fitvoltref)
    voltdecayC82, rsqrtC82 = linear_fit(fittimeC82, fitvoltC82)

    voltdecaypulse, rsqrtpulse = linear_fit(fittimepulse, fitvoltpulse)
    voltdecaypulseC5, rsqrtpulseC5 = linear_fit(fittimepulseC5, fitvoltpulseC5)
    voltdecaypulseC575s, rsqrtpulseC575s = linear_fit(fittimepulseC575s, fitvoltpulseC575s)
    voltdecaypulsedouble, rsqrtpulsedouble = linear_fit(fittimepulsedouble, fitvoltpulsedouble)
    voltdecaypulseC57CC20, rsqrtpulseC57CC20 = linear_fit(fittimepulseC57CC20, fitvoltpulseC57CC20)

    fig,ax = plt.subplots(1, 1, figsize=(4, 2), sharex=True)
    plt.style.use(['science','nature','no-latex'])
    ax.plot(list(np.linspace(1,len(rsqrtref),len(rsqrtref))), rsqrtref, color = 'r')
    ax.plot(list(np.linspace(1,len(rsqrtpulse),len(rsqrtpulse))), rsqrtpulse, color = 'cornflowerblue')
    ax.plot(list(np.linspace(1,len(rsqrtpulsedouble),len(rsqrtpulsedouble))), rsqrtpulsedouble, color = 'cadetblue')
    ax.plot(list(np.linspace(1,len(rsqrtC82),len(rsqrtC82))), rsqrtC82, color = 'steelblue')
        
    ax.plot(list(np.linspace(1,len(rsqrtpulseC57CC20),len(rsqrtpulseC57CC20))), rsqrtpulseC57CC20, color = 'darkgreen')
    ax.plot(list(np.linspace(1,len(rsqrtpulseC575s),len(rsqrtpulseC575s))), rsqrtpulseC575s, color = 'darkolivegreen')
    ax.plot(list(np.linspace(1,len(rsqrtpulseC5),len(rsqrtpulseC5))), rsqrtpulseC5, color = 'forestgreen')

    ax.scatter(list(np.linspace(1,len(rsqrtref),len(rsqrtref))), rsqrtref, color = 'r')
    ax.scatter(list(np.linspace(1,len(rsqrtpulse),len(rsqrtpulse))), rsqrtpulse, color = 'cornflowerblue')
    ax.scatter(list(np.linspace(1,len(rsqrtpulsedouble),len(rsqrtpulsedouble))), rsqrtpulsedouble, color = 'cadetblue')
    ax.scatter(list(np.linspace(1,len(rsqrtC82),len(rsqrtC82))), rsqrtC82, color = 'steelblue')
        
    ax.scatter(list(np.linspace(1,len(rsqrtpulseC57CC20),len(rsqrtpulseC57CC20))), rsqrtpulseC57CC20, color = 'darkgreen')
    ax.scatter(list(np.linspace(1,len(rsqrtpulseC575s),len(rsqrtpulseC575s))), rsqrtpulseC575s, color = 'darkolivegreen')
    ax.scatter(list(np.linspace(1,len(rsqrtpulseC5),len(rsqrtpulseC5))), rsqrtpulseC5, color = 'forestgreen')

    #decimals = 2  # Change this value as needed
    #ax.yaxis.set_major_formatter(StrMethodFormatter(f'{{x:.{decimals}f}}'))

    ax.set_xlabel('cell number', fontsize=8)
    ax.set_ylabel(r'R$^2$ of linear fit', fontsize=8)

    plt.legend(['Reference CC C/20','Pulsed C/5.7' , 'Pulsed double length',  'CC C/8.23',  'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s' ],
                fontsize = 4, bbox_to_anchor=(0.6, 0.25), borderaxespad=0,frameon=True, ncol = 2 )
    ax.set_xticks([1,  2, 3, 4, 5]) 
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('SelfdischargelinearfitsqrtT-V.png', dpi=800)
    os.chdir(origpath)
    plt.show()


    fig,ax = plt.subplots(1, 1, figsize=(4, 2), sharex=True)
    plt.style.use(['science','nature','no-latex'])
    size = 5
    ax.scatter(list(np.linspace(1,len(voltdecayref),len(voltdecayref))), voltdecayref, s = size, color = 'r')
    
    ax.scatter(list(np.linspace(1,len(voltdecaypulse),len(voltdecaypulse))), voltdecaypulse, s = size, color = 'cornflowerblue')
    ax.scatter(list(np.linspace(1,len(voltdecaypulsedouble),len(voltdecaypulsedouble))), voltdecaypulsedouble, s = size, color = 'cadetblue')
    ax.scatter(list(np.linspace(1,len(voltdecayC82),len(voltdecayC82))), voltdecayC82, s = size, color = 'steelblue')
    
    ax.scatter(list(np.linspace(1,len(voltdecaypulseC57CC20),len(voltdecaypulseC57CC20))), voltdecaypulseC57CC20, s = size, color = 'darkgreen')
    ax.scatter(list(np.linspace(1,len(voltdecaypulseC575s),len(voltdecaypulseC575s))), voltdecaypulseC575s, s = size, color = 'darkolivegreen')
    ax.scatter(list(np.linspace(1,len(voltdecaypulseC5),len(voltdecaypulseC5))), voltdecaypulseC5, s = size, color = 'forestgreen')
    
    ax.plot(list(np.linspace(1,len(voltdecayref),len(voltdecayref))), voltdecayref, color = 'r')
    
    ax.plot(list(np.linspace(1,len(voltdecaypulse),len(voltdecaypulse))), voltdecaypulse, color = 'cornflowerblue')
    ax.plot(list(np.linspace(1,len(voltdecaypulsedouble),len(voltdecaypulsedouble))), voltdecaypulsedouble, color = 'cadetblue')
    ax.plot(list(np.linspace(1,len(voltdecayC82),len(voltdecayC82))), voltdecayC82, color = 'steelblue')
    
    
    ax.plot(list(np.linspace(1,len(voltdecaypulseC57CC20),len(voltdecaypulseC57CC20))), voltdecaypulseC57CC20, color = 'darkgreen')
    ax.plot(list(np.linspace(1,len(voltdecaypulseC575s),len(voltdecaypulseC575s))), voltdecaypulseC575s, color = 'darkolivegreen')
    ax.plot(list(np.linspace(1,len(voltdecaypulseC5),len(voltdecaypulseC5))), voltdecaypulseC5, color = 'forestgreen')
    
    ax.set_xlabel('cell number', fontsize=8)
    ax.set_ylabel(r'voltage decay [mV/$\sqrt{h}$]', fontsize=8)

    #decimals = 4  # Change this value as needed
    #ax.yaxis.set_major_formatter(StrMethodFormatter(f'{{x:.{decimals}f}}'))
    plt.ylim([-2, -4.6])
    plt.xticks([1,  2, 3, 4, 5]) 
    plt.legend(['Reference CC C/20','Pulsed C/5.7' , 'Pulsed double length',  'CC C/8.23',  'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s' ],
                fontsize = 3.5, bbox_to_anchor=(0.7, 0.2), frameon=True, ncol = 3 )
      
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('SelfdischargeVoltageDecay5cells.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def EOLIRCheck():

    area = np.around(0.7*0.7*np.pi,4)
    
    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        irvalueCC = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvalueCC.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvalueCC

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvalueC83 = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])

                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvalueC83.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvalueC83


    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        irvaluePulse = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluePulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 

        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluePulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area)  

        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluePulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvaluePulse

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluepulseC5 = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])

                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluepulseC5.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])) 
        return irvaluepulseC5

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluepulseC575s = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluepulseC575s.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvaluepulseC575s
    
    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluedoublepulse = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluedoublepulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluedoublepulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvaluedoublepulse

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvalueC57CC20 = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvalueC57CC20.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvalueC57CC20

    irvaluesCC = getCCref()

    irvaluesPulse= getpulseform()

    irvaluesC83 = getCCrefC83()

    irvaluesPulseC575s = getpulseformC575s()

    irvaluesPulseC57CC20 = getpulseformC57CC20()

    irvaluesdoublePulse = getpulseformdoublepulse()
    
    irvaluesPulseC5 = getpulseformC5()
    
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1,1)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)
    fig.tight_layout(w_pad=1.0)
    fontsizelabelaxis = 4
    medianwidth = 0.5
    markersize = 2
    marker = 'D'
    boxlinewidth = 0.5


    CCref=ax.boxplot(irvaluesCC, positions = [0], patch_artist=True, boxprops=dict(facecolor='r', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    
    pulse = ax.boxplot(irvaluesPulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulsedoublepuls = ax.boxplot(irvaluesdoublePulse, positions = [1.5], patch_artist=True, boxprops=dict(facecolor='cadetblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    CCrefC83=ax.boxplot(irvaluesC83, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    
    pulseC5 = ax.boxplot(irvaluesPulseC5, positions = [3], patch_artist=True, boxprops=dict(facecolor='forestgreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulseC57CC20 = ax.boxplot(irvaluesPulseC57CC20, positions = [3.5], patch_artist=True, boxprops=dict(facecolor='darkgreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulseC5s = ax.boxplot(irvaluesPulseC575s, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    

    #ax[x[o],y[o]].set_xlim([-0.5,2])
    #ax[x[o],y[o]].set_ylim([0.5,1.2])
    ax.set_xticks([])  

    legendcolors = [CCref["boxes"][0],pulse["boxes"][0],  pulsedoublepuls["boxes"][0], CCrefC83["boxes"][0], 
                         pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0]]
    legendtext = ['Reference CC C/20','Pulsed C/5.7' ,  'Pulsed double length', 'CC C/8.2',  'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s']

    legend = ax.legend(legendcolors, legendtext, fontsize = 3, bbox_to_anchor=(0.44, 0.95),
                             frameon=True, borderaxespad=0, ncol = 2 )
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    #ax.set_ylim([0,35])
    ax.set_ylabel(r'Area specific impedance [Ohm / cm$^2$]', fontsize=8)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('IRcheck5cells.png', dpi=1000)
    os.chdir(origpath)
    plt.show()

def EOLIRCheckBrokenAxis():

    area = np.around(0.7*0.7*np.pi,4)
    
    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        irvalueCC = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvalueCC.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvalueCC

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvalueC83 = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])

                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvalueC83.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvalueC83


    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        irvaluePulse = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluePulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 

        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluePulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area)  

        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluePulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvaluePulse

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluepulseC5 = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])

                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluepulseC5.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])) 
        return irvaluepulseC5

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluepulseC575s = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluepulseC575s.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvaluepulseC575s
    
    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluedoublepulse = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluedoublepulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvaluedoublepulse.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvaluedoublepulse

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvalueC57CC20 = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                if len(cycles) == 38:
                    IRcycles = 16
                else:
                    IRcycles = 15
                try:
                        voltlist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(IRcycles)]["Step_"+str(IRcycles)+'_Testtime(s)']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]["Step_"+str(IRcycles)+'_Voltage(V)']
                except:
                        voltlist = data['split']["Step_"+str(IRcycles)]['V']
                        voltlistrest = data['split']["Step_"+str(IRcycles-1)]['V']
                        timelist = data['split']["Step_"+str(IRcycles)]['t']
                        currlist = data['split']["Step_"+str(IRcycles)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
                irvalueC57CC20.append(abs((voltlistrest[-1]- voltlist[-1])/currlist[2])/area) 
        return irvalueC57CC20

    irvaluesCC = getCCref()

    irvaluesPulse= getpulseform()

    irvaluesC83 = getCCrefC83()

    irvaluesPulseC575s = getpulseformC575s()

    irvaluesPulseC57CC20 = getpulseformC57CC20()

    irvaluesdoublePulse = getpulseformdoublepulse()
    
    irvaluesPulseC5 = getpulseformC5()
    
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 2,gridspec_kw={'width_ratios': [1, 4]}, sharey=True)
    plt.style.use(['science','nature','no-latex'])
    plt.subplots_adjust(wspace=0.05)

    fontsizelabelaxis = 4
    medianwidth = 0.5
    markersize = 2
    marker = 'D'
    boxlinewidth = 0.5

    CCref=ax[0].boxplot(irvaluesCC, positions = [24], patch_artist=True, boxprops=dict(facecolor='r', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    
    pulse = ax[1].boxplot(irvaluesPulse, positions = [9.7], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulsedoublepuls = ax[1].boxplot(irvaluesdoublePulse, positions = [9.4], patch_artist=True, boxprops=dict(facecolor='cadetblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    CCrefC83=ax[1].boxplot(irvaluesC83, positions = [9.2], patch_artist=True, boxprops=dict(facecolor='steelblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    
    pulseC5 = ax[1].boxplot(irvaluesPulseC5, positions = [8.4], patch_artist=True, boxprops=dict(facecolor='forestgreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulseC57CC20 = ax[1].boxplot(irvaluesPulseC57CC20, positions = [8.2], patch_artist=True, boxprops=dict(facecolor='darkgreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulseC5s = ax[1].boxplot(irvaluesPulseC575s, positions = [8.1], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth,color = 'limegreen'),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    

    legendcolors = [CCref["boxes"][0],pulse["boxes"][0], pulsedoublepuls["boxes"][0],   CCrefC83["boxes"][0], pulseC5["boxes"][0], pulseC57CC20["boxes"][0], pulseC5s["boxes"][0] ]
    #legendtext = ['N(CC20) = ' + str(len(refcoulomb4Form)),  'N(PulsC5.7) = ' + str(len(pulsecoulomb4Form)),'N(CC8.23) = ' + str(len(refcoulomb4FormC8)),
    #                   'N(PulsC5.7 5s) = ' + str(len(pulsecoulomb4FormC5s)), 'N(PulsC5.7 CC20) = ' + str(len(pulsecoulomb4FormC57CC20)),
    #                    'N(PulsC5.7double) = ' + str(len(pulsecoulomb4Formdoublepulse)), 'N(PulsC5) = ' + str(len(pulsecoulomb4FormC5))]
    legendtext = ['Reference CC C/20', 'Pulsed C/5.7' ,'Pulsed double length',  'CC C/8.2', 'Pulsed C/5', 'Mixed C/5.7+C/20', 'Pulsed 5s']
    legend = plt.legend(legendcolors, legendtext, fontsize=4, bbox_to_anchor=(-0.25, 0.9), loc="upper left",
                 borderaxespad=0, fancybox=True, shadow=True,frameon=True, ncol = 2)
    xlabels = ['24.0', '9.7', '9.4', '9.2', '8.4', '8.2', '8.1']
    #ax.set_ylim([0.76,0.86])

    #ax.set_xticks([0, 1, 1.5, 2,  3, 3.5, 4]) 
    ax[0].set_xticklabels([xlabels[0]])
    ax[1].set_xticklabels(xlabels[1:])

    d = .02  # Größe der diagonalen Linien in Bruchteilen der Achsengröße
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)        # Diagonale Bruchlinie rechts
    #kwargs = dict(transform=ax[1].transAxes, color='k', clip_on=False)
    ax[1].plot((1.1-d, 1.1+d), (-d, +d), **kwargs)          # Diagonale Bruchlinie links

    #ax[0].spines['bottom'].set_visible(False)
    #ax[1].spines['bottom'].set_visible(False)
    ax[1].invert_xaxis()
    ax[1].yaxis.set_visible(False)

    ax[0].tick_params(axis='x', labelsize=8, rotation=45)
    ax[1].tick_params(axis='x', labelsize=8, rotation=45)

    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[1].spines['top'].set_visible(False)

    ax[0].xaxis.set_ticks_position('bottom')
    ax[1].xaxis.set_ticks_position('bottom')

    ax[0].minorticks_off()
    ax[1].minorticks_off()

        
    ax[0].set_ylabel(r'area specific impedance [Ohm / cm$^2$]', fontsize=8)
    #ax[1].set_xlabel('Mean charge formation time [h]', fontsize=8)
    fig.text(0.5, -0.04, 'mean wall time for formation charge [h]', ha='center', fontsize=8)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('IRcheck5cellsbrokenaxis.png', dpi=1000)
    os.chdir(origpath)
    plt.show()

def ExplainSchema():

    pathpulse = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\PulsedC58"

    hdf5path = pathpulse
    namecriteria = '_PulsedC58Temp'

    os.chdir(hdf5path)
    originpath = os.getcwd()
    content = os.listdir()
    content.sort(key=natural_keys)
    for i in range(len(content)):
        if namecriteria in content[i]:
            print(content[i])
            data = loadHDF5(content[i])
            currentlist = data['raw']["Current(A)"]
            timelist = data['raw']["Testtime(s)"]
            voltlist = data['raw']["Voltage(V)"]
            break
    
    pathpulse = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\PulsedC58\LongtermC58"
    hdf5path = pathpulse
    namecriteria = '_PulsedC58Temp'
    os.chdir(hdf5path)
    originpath = os.getcwd()
    content = os.listdir()
    content.sort(key=natural_keys)
    for i in range(len(content)):
        if namecriteria in content[i]:
            print(content[i])
            data = loadHDF5(content[i])
            currentlistlong = data['raw']["Current(A)"]
            timelistlong = data['raw']["Testtime(s)"]
            voltlistlong = data['raw']["Voltage(V)"]
            break
    
    for i in range(len(timelistlong)):
        timelistlong[i] = (timelistlong[i]+timelist[-1])/(3600*24)
    
    for i in range(len(timelist)):
        timelist[i] = (timelist[i])/(3600*24)

    for i in range(len(currentlist)):
        currentlist[i] = (currentlist[i])*1000
    
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 5,gridspec_kw={'width_ratios': [2, 2, 2, 1, 1]}, sharey=True)
    plt.style.use(['science','nature','no-latex'])
    plt.subplots_adjust(wspace=0.05)
   
    slicevalueform = 179000
    slicevalueself1 = 188350
    sliceself= 194500
    slicevalue = 51800

    ax[0].plot(timelist[:slicevalueform], voltlist[:slicevalueform], color = 'r')

    ax[1].plot(timelist[slicevalueform:slicevalueself1], voltlist[slicevalueform:slicevalueself1], color = 'b')

    ax[2].plot(timelist[sliceself:], voltlist[sliceself:], color = 'b')

    ax[3].plot(timelistlong[:slicevalue], voltlistlong[:slicevalue], color = 'g')

    ax[4].plot(timelistlong[-slicevalue:], voltlistlong[-slicevalue:], color = 'g')

    ax_inset = inset_axes(ax[0], width="60%", height="60%", bbox_to_anchor=(0.28, -0.06, 1.0, 0.6), bbox_transform=ax[0].transAxes)

    x1 = 9998
    x2 = 10021
    labelsize = 6
    labelpadval = 0.5
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append(timelist[l]*(3600*24)-timelist[x1]*(3600*24))
    ax_inset.plot(pulsetime, voltlist[x1:x2], color = 'r')
    ax_inset.set_ylabel('Voltage [V]', fontsize=labelsize, labelpad=labelpadval)
    ax_inset.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval)

    ax_inset2 = inset_axes(ax[1], width="60%", height="60%", bbox_to_anchor=(1.48, -0.06, 1.0, 0.6), bbox_transform=ax[0].transAxes)

    x1 = 9998
    x2 = 10021
    labelsize = 6
    labelpadval = 0.5
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append(timelist[l]*(3600*24)-timelist[x1]*(3600*24))
    ax_inset2.plot(pulsetime, currentlist[x1:x2], color = 'dodgerblue', linewidth = 2)
    ax_inset2.set_ylabel('current [mA]', fontsize=labelsize, labelpad=labelpadval)
    ax_inset2.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval)
    ax_inset2.set_ylim([0,1])
    #ax_inset2.set_yticks([0,0.2,0.4,0.6,0.8,1])

    ax_inset.minorticks_off()  # Disable minor ticks on primary axis
    ax_inset2.minorticks_off()  #

    ax_inset.tick_params(axis='both', which='minor', bottom=False, left=False)
    ax_inset2.tick_params(axis='both', which='minor', bottom=False, left=False)

    ax_inset.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_inset2.tick_params(axis='both', which='major', labelsize=labelsize)
    #ax.indicate_inset_zoom(ax_inset, edgecolor="black")

    d = .02  # Größe der diagonalen Linien in Bruchteilen der Achsengröße
    kwargs = dict(transform=ax[1].transAxes, color='k', clip_on=False)
    ax[1].plot((1-d, 1+d), (-d, +d), **kwargs)        # Diagonale Bruchlinie rechts
    #kwargs = dict(transform=ax[1].transAxes, color='k', clip_on=False)
    ax[2].plot((1.05-d, 1.05+d), (-d, +d), **kwargs)          # Diagonale Bruchlinie links

    kwargs2 = dict(transform=ax[3].transAxes, color='k', clip_on=False)
    ax[3].plot((1-d, 1+d), (-d, +d), **kwargs2)        # Diagonale Bruchlinie rechts
    #kwargs = dict(transform=ax[1].transAxes, color='k', clip_on=False)
    ax[4].plot((1.1-d, 1.1+d), (-d, +d), **kwargs2)          # Diagonale Bruchlinie links

    ax[2].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #ax[0].spines['bottom'].set_visible(False)
    #ax[1].spines['bottom'].set_visible(False)
    ax[1].yaxis.set_visible(False)
    ax[2].yaxis.set_visible(False)
    ax[3].yaxis.set_visible(False)
    ax[4].yaxis.set_visible(False)

    ax[0].tick_params(axis='x', labelsize=8, rotation=45)
    ax[1].tick_params(axis='x', labelsize=8, rotation=45)
    ax[2].tick_params(axis='x', labelsize=8, rotation=45)
    ax[3].tick_params(axis='x', labelsize=8, rotation=45)
    ax[4].tick_params(axis='x', labelsize=8, rotation=45)

    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[2].spines['left'].set_visible(False)
    ax[2].spines['right'].set_visible(False)
    ax[3].spines['left'].set_visible(False)
    ax[3].spines['right'].set_visible(False)
    ax[4].spines['left'].set_visible(False)
    ax[4].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[3].spines['top'].set_visible(False)
    ax[4].spines['top'].set_visible(False)

    ax[0].xaxis.set_ticks_position('bottom')
    ax[1].xaxis.set_ticks_position('bottom')
    ax[2].xaxis.set_ticks_position('bottom')
    ax[3].xaxis.set_ticks_position('bottom')
    ax[4].xaxis.set_ticks_position('bottom')
       
    ax[0].set_ylabel('Voltage [V]', fontsize=6)
    fig.text(0.2, -0.03, 'Formation', ha='center', fontsize=6, color = 'r')
    fig.text(0.5, -0.03, 'End-of-Line Test', ha='center', fontsize=6, color = 'b')
    fig.text(0.8, -0.03, 'longterm cycling', ha='center', fontsize=6, color = 'g')
    fig.text(0.5, -0.08, 'time [d]', ha='center', fontsize=8)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('ExplainSchema.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def ExplainSchema2():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)

        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
        for i in range(len(timelist)):
            timelist[i] = timelist[i]/ 3600
        for i in range(len(currlist)):
            currlist[i] = currlist[i]*1000
        return voltlist, timelist, currlist

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)

        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])

                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']
                        #print(currlist[0])
                        #plt.scatter(timelist, voltlist)
                        #plt.show()
                        #print(voltlist[0]-voltlist[1])
        for i in range(len(timelist)):
            timelist[i] = timelist[i]/ 3600
        for i in range(len(currlist)):
            currlist[i] = currlist[i]*1000
        return voltlist, timelist, currlist


    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']

        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I'] 

        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']
        for i in range(len(timelist)):
            timelist[i] = timelist[i]/ 3600
        for i in range(len(currlist)):
            currlist[i] = currlist[i]*1000
        return voltlist, timelist, currlist

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])

                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']
        for i in range(len(timelist)):
            timelist[i] = timelist[i]/ 3600
        for i in range(len(currlist)):
            currlist[i] = currlist[i]*1000
        return voltlist, timelist, currlist

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluepulseC575s = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(4,len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']
                break
        for i in range(len(timelist)):
            timelist[i] = timelist[i]/ 3600
        for i in range(len(currlist)):
            currlist[i] = currlist[i]*1000
        return voltlist, timelist, currlist
    
    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvaluedoublepulse = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']

        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']
        for i in range(len(timelist)):
            timelist[i] = timelist[i]/ 3600
        for i in range(len(currlist)):
            currlist[i] = currlist[i]*1000
        return voltlist, timelist, currlist

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        irvalueC57CC20 = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                form = 1
                try:
                        voltlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Voltage(V)']
                        timelist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Testtime(s)']
                        currlist = data['split']["Step_"+str(form)]["Step_"+str(form)+'_Current(A)']
                except:
                        voltlist = data['split']["Step_"+str(form)]['V']
                        timelist = data['split']["Step_"+str(form)]['t']
                        currlist = data['split']["Step_"+str(form)]['I']
        for i in range(len(timelist)):
            timelist[i] = timelist[i]/ 3600
        for i in range(len(currlist)):
            currlist[i] = currlist[i]*1000
        return voltlist, timelist, currlist



    voltref, timeref, curref = getCCref()
    voltC82, timeC82, curC82 = getCCrefC83()

    voltpulseC57, timepulseC57, currpulseC57 = getpulseform()
    voltpulseC5, timepulseC5, currpulseC5 = getpulseformC5()
    voltpulsedouble, timepulsedouble, currpulsedouble = getpulseformdoublepulse()
    voltpulseC575s, timepulseC575s, currpulseC575s = getpulseformC575s()
    voltpulseC57CC20, timepulseC57CC20, currpulseC57CC20 = getpulseformC57CC20()

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(2, 3)
    plt.style.use(['science','nature','no-latex'])
    plt.subplots_adjust(wspace=0.25, hspace = 0.25)
    color = ['cadetblue', 'cornflowerblue', 'steelblue', 'forestgreen', 'darkgreen', 'darkolivegreen']
   
    ax[0,0].plot(timeref, voltref, color = 'r', label = 'CC C/20')
    ax[0,0].plot(timeC82, voltC82, color = color[2], label = 'CC C/8.2')
    ax[0,0].legend(fontsize = 2, bbox_to_anchor=(1.0, 0.05), loc="lower right",)

    ax[0,1].plot(timepulseC57, voltpulseC57, color = color[1])

    ax[0,2].plot(timepulsedouble, voltpulsedouble, color = color[0])

    ax[1,0].plot(timepulseC5, voltpulseC5, color = color[3])

    ax[1,1].plot(timepulseC57CC20, voltpulseC57CC20, color = color[4])

    ax[1,2].plot(timepulseC575s, voltpulseC575s, color = color[5])

    ax_inset = inset_axes(ax[0,0], width="60%", height="60%", bbox_to_anchor=(0.02, -0.02, 0.8, 0.8), bbox_transform=ax[0,0].transAxes)
    ax2 = ax_inset.twinx()
    x1 = 9998
    x2 = 10001
    labelsize = 4
    labelpadval = 0.4
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timeref[l]*3600)-(timeref[x1]*3600))
    ax2.set_ylim([-0.05,1])
    ax2.plot(pulsetime, curref[x1:x2], color = 'k', alpha = 0.6)
    ax2.set_ylabel('current [mA]', fontsize=labelsize, labelpad=labelpadval)
    ax2.tick_params(axis='y', labelsize=labelsize-1, pad = labelpadval)  
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timeref[l]*3600)-(timeref[x1]*3600))
    ax_inset.plot(pulsetime, list(savgol_filter(voltref[x1:x2], 3, 1)) , color = 'r')
    ax_inset.set_ylabel('Voltage [V]', fontsize=labelsize, labelpad=labelpadval-1.4, color='r')
    ax_inset.set_yticklabels([])
    ax_inset.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval-1.4)
    ax_inset.tick_params(axis='y', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.tick_params(axis='x', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.yaxis.tick_left()  # Only show ticks on the left for the primary y-axis
    ax2.yaxis.tick_right() 

    ax_inset = inset_axes(ax[0,2], width="60%", height="60%", bbox_to_anchor=(0.02, -0.02, 0.8, 0.8), bbox_transform=ax[0,2].transAxes)
    ax2 = ax_inset.twinx()
    x1 = 9998
    x2 = 10050
    labelsize = 4
    labelpadval = 0.4
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulsedouble[l]*3600)-(timepulsedouble[x1]*3600))
    ax2.set_ylim([-0.05,1])
    ax2.plot(pulsetime, currpulsedouble[x1:x2], color = 'k', alpha = 0.6)
    ax2.set_ylabel('current [mA]', fontsize=labelsize, labelpad=labelpadval)
    ax2.tick_params(axis='y', labelsize=labelsize-1, pad = labelpadval)  
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulsedouble[l]*3600)-(timepulsedouble[x1]*3600))
    ax_inset.plot(pulsetime, voltpulsedouble[x1:x2], color = color[0])
    ax_inset.set_ylabel('Voltage [V]', fontsize=labelsize, labelpad=labelpadval-1.4, color=color[0])
    ax_inset.set_yticklabels([])
    ax_inset.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval-1.4)
    ax_inset.tick_params(axis='y', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.tick_params(axis='x', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.yaxis.tick_left()  # Only show ticks on the left for the primary y-axis
    ax2.yaxis.tick_right()  

    ax_inset = inset_axes(ax[0,1], width="60%", height="60%", bbox_to_anchor=(0.02, -0.02, 0.8, 0.8), bbox_transform=ax[0,1].transAxes)
    ax2 = ax_inset.twinx()
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC57[l]*3600)-(timepulseC57[x1]*3600))
    ax2.set_ylim([-0.05,1])
    ax2.plot(pulsetime, currpulseC57[x1:x2], color = 'k', alpha = 0.6)
    ax2.set_ylabel('current [mA]', fontsize=labelsize, labelpad=labelpadval)
    ax2.tick_params(axis='y', labelsize=labelsize-1, pad = labelpadval)  
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC57[l]*3600)-(timepulseC57[x1]*3600))
    ax_inset.plot(pulsetime, voltpulseC57[x1:x2], color = color[1])
    ax_inset.set_ylabel('Voltage [V]', fontsize=labelsize, labelpad=labelpadval-1.4, color=color[1])
    ax_inset.set_yticklabels([])
    ax_inset.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval-1.4)
    ax_inset.tick_params(axis='y', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.tick_params(axis='x', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.yaxis.tick_left()  # Only show ticks on the left for the primary y-axis
    ax2.yaxis.tick_right()  

    ax_inset = inset_axes(ax[1,0], width="60%", height="60%", bbox_to_anchor=(0.02, -0.02, 0.8, 0.8), bbox_transform=ax[1,0].transAxes)
    ax2 = ax_inset.twinx()
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC5[l]*3600)-(timepulseC5[x1]*3600))
    ax2.set_ylim([-0.05,1])
    ax2.plot(pulsetime, currpulseC5[x1:x2], color = 'k', alpha = 0.6)
    ax2.set_ylabel('current [mA]', fontsize=labelsize, labelpad=labelpadval)
    ax2.tick_params(axis='y', labelsize=labelsize-1, pad = labelpadval)  
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC5[l]*3600)-(timepulseC5[x1]*3600))
    ax_inset.plot(pulsetime, voltpulseC5[x1:x2], color = color[3])
    ax_inset.set_ylabel('Voltage [V]', fontsize=labelsize, labelpad=labelpadval-1.4, color=color[3])
    ax_inset.set_yticklabels([])
    ax_inset.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval-1.4)
    ax_inset.tick_params(axis='y', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.tick_params(axis='x', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.yaxis.tick_left()  # Only show ticks on the left for the primary y-axis
    ax2.yaxis.tick_right()

    ax_inset = inset_axes(ax[1,1], width="60%", height="60%", bbox_to_anchor=(0.02, -0.02, 0.8, 0.8), bbox_transform=ax[1,1].transAxes)
    ax2 = ax_inset.twinx()
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC57CC20[l]*3600)-(timepulseC57CC20[x1]*3600))
    ax2.set_ylim([-0.05,1])
    ax2.plot(pulsetime, currpulseC57CC20[x1:x2], color = 'k', alpha = 0.6)
    ax2.set_ylabel('current [mA]', fontsize=labelsize, labelpad=labelpadval)
    ax2.tick_params(axis='y', labelsize=labelsize-1, pad = labelpadval)  
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC57CC20[l]*3600)-(timepulseC57CC20[x1]*3600))
    ax_inset.plot(pulsetime, voltpulseC57CC20[x1:x2], color = color[4])
    ax_inset.set_ylabel('Voltage [V]', fontsize=labelsize, labelpad=labelpadval-1.4, color=color[4])
    ax_inset.set_yticklabels([])
    ax_inset.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval-1.4)
    ax_inset.tick_params(axis='y', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.tick_params(axis='x', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.yaxis.tick_left()  # Only show ticks on the left for the primary y-axis
    ax2.yaxis.tick_right()

    ax_inset = inset_axes(ax[1,2], width="60%", height="60%", bbox_to_anchor=(0.02, -0.02, 0.8, 0.8), bbox_transform=ax[1,2].transAxes)
    ax2 = ax_inset.twinx()
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC575s[l]*3600)-(timepulseC575s[x1]*3600))
    ax2.set_ylim([-0.05,1])
    ax2.plot(pulsetime, currpulseC575s[x1:x2], color = 'k', alpha = 0.6)
    ax2.set_ylabel('current [mA]', fontsize=labelsize, labelpad=labelpadval)
    ax2.tick_params(axis='y', labelsize=labelsize-1, pad = labelpadval)  
    pulsetime = []
    for l in range(x1,x2):
        pulsetime.append((timepulseC5[l]*3600)-(timepulseC5[x1]*3600))
    ax_inset.plot(pulsetime, voltpulseC575s[x1:x2], color = color[5])
    ax_inset.set_ylabel('Voltage [V]', fontsize=labelsize, labelpad=labelpadval-1.4, color=color[5])
    ax_inset.set_yticklabels([])
    ax_inset.set_xlabel('time [s]', fontsize=labelsize, labelpad=labelpadval-1.4)
    ax_inset.tick_params(axis='y', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.tick_params(axis='x', labelsize=labelsize-1, pad=labelpadval)
    ax_inset.yaxis.tick_left()  # Only show ticks on the left for the primary y-axis
    ax2.yaxis.tick_right() 

    legendtext = [ 'CC C/20 + CC C/8.2','Pulsed C/5.7' ,'Pulsed double length', 'Pulsed C/5',  'Mixed C/5.7+C/20', 'Pulsed 5s' ]
    #ax[1, 0].legend(lines, legendtext, fontsize=3, bbox_to_anchor=(-0.025, 1.175), loc="upper left", ncol=6, handletextpad=0.25)
    ax[0,0].set_title(legendtext[0], fontsize=4, pad=3, loc='center')
    ax[0,1].set_title(legendtext[1], fontsize=4, pad=3, loc='center')
    ax[0,2].set_title(legendtext[2], fontsize=4, pad=3, loc='center')
    ax[1,0].set_title(legendtext[3], fontsize=4, pad=3, loc='center')
    ax[1,1].set_title(legendtext[4], fontsize=4, pad=3, loc='center')
    ax[1,2].set_title(legendtext[5], fontsize=4, pad=3, loc='center')
    fig.text(0.04, 0.5, 'Voltage [V]', ha='center', va='center', rotation='vertical', fontsize=6)
    ax[1,1].set_xlabel('Time [h]', fontsize=6)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('ExplainSchema2.png', dpi=800)
    os.chdir(origpath)
    plt.show()


def EOLselfdischargeVoltcourse_all():
    
    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        CCcapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa
    
    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    selfdischarge, selfdccapa, voltlistref, timelistref, currentlistC2ref, voltlistC2ref, cellnamesref, capalistCC20 = getCCref()

    selfdischargepulse, selfdccapapulse, voltlistpulse, timelistpulse, currentlistC2pulse, voltlistC2pulse, cellnamespulse, capalistC57 = getpulseform()

    selfdischargeC82, selfdccapaC82, voltlistC82, timelistC82, currentlistCC2C82, voltlistCC2C82, cellnamesC82, capalistC82 = getCCrefC83()

    selfdischargepulseC5, selfdccapapulseC5, voltlistpulseC5, timelistpulseC5, currentlistC2pulseC5, voltlistC2pulseC5, cellnamespulseC5, capalistC5 = getpulseformC5()

    selfdischargepulseC575s, selfdccapapulseC575s, voltlistpulseC575s, timelistpulseC575s, currentlistC2pulseC575s, voltlistC2pulseC575s, cellnamespulseC575s, capalistC575s = getpulseformC575s()

    selfdischargepulsedouble, selfdccapapulsedouble, voltlistpulsedouble, timelistpulsedouble, currentlistC2pulsedouble, voltlistC2pulsedouble, cellnamespulsedouble, capalistdouble = getpulseformdoublepulse()

    selfdischargepulseC57CC20, selfdccapapulseC57CC20, voltlistpulseC57CC20, timelistpulseC57CC20, currentlistC2pulseC57CC20, voltlistC2pulseC57CC20, cellnamespulseC57CC20, capalistC57CC20 = getpulseformC57CC20()


    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)
    fig.tight_layout()
    fontsizelabelaxis = 6

    #a)

    def change_time_h(timelist):
        newtime = []
        for p in range(len(timelist)):
            temptime = []
            for time in range(len(timelist[p])):
                temptime.append(timelist[p][time]/(60*60*24))
            newtime.append(temptime)
        return newtime

    newtimeref = change_time_h(timelistref)
    newtimeC82 = change_time_h(timelistC82)
    newtimepulse = change_time_h(timelistpulse)
    newtimepulseC5 = change_time_h(timelistpulseC5)
    newtimepulseC575s = change_time_h(timelistpulseC575s)
    newtimepulsedouble = change_time_h(timelistpulsedouble)
    newtimepulseC57CC20 = change_time_h(timelistpulseC57CC20)    
    
    ax.plot(newtimeref[-1], voltlistref[-1], label = 'Reference CC C/20', color = "r")
    ax.plot(newtimepulse[-1], voltlistpulse[-1], label = 'Pulsed C/5.7', color = "cadetblue")
    ax.plot(newtimepulsedouble[1], voltlistpulsedouble[1], label = 'Pulsed C/5.7 double length', color = "cornflowerblue")
    ax.plot(newtimeC82[-1], voltlistC82[-1], label = 'CC C/8.23', color = "steelblue")
    
    ax.plot(newtimepulseC5[-1], voltlistpulseC5[-1], label = 'Pulsed C/5', color = "forestgreen")
    ax.plot(newtimepulseC57CC20[-2], voltlistpulseC57CC20[-2], label = 'Mixed C/5.7 CC C/20', color = "darkgreen")
    ax.plot(newtimepulseC575s[-1], voltlistpulseC575s[-1], label = 'Pulse 5s', color = "darkolivegreen")
     
    ax.legend(bbox_to_anchor=(0.98, 0.95),borderaxespad=0, ncol = 2, fontsize = 3)
    
    ax.set_ylabel('voltage [V]', fontsize=fontsizelabelaxis)
    ax.set_xlabel('time [day]', fontsize=fontsizelabelaxis)
    ax.set_xticks([1,2,3,4,5])
    ax.yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax.xaxis.set_tick_params(labelsize=fontsizelabelaxis)

    ax.axhspan(4.11, 4.192, xmin=0.020, xmax=0.225 , color='#008000', alpha=0.4) # red
    ax.text(0.02, 0.98, "Voltage relaxation", transform=ax.transAxes, fontsize=4,
        verticalalignment='top',weight='bold', color = 'k')
    ax.axhspan(4.11, 4.192, xmin=0.0225, xmax=0.99 , color='m', alpha=0.2) # magenta
    ax.text(0.3, 0.98, "Self-discharge", transform=ax.transAxes, fontsize=4,
        verticalalignment='top', weight='bold',color = 'k')

    '''
    ############################################# Course too similar to regular Voltage vs time
    ax2 = ax.twinx()
    dVdtref = []
    for i in range(len(voltlistref)):
        tempdvdt=[]
        for l in range(1, len(voltlistref[i])):
            tempdvdt.append(voltlistref[i][l]-voltlistref[i][l-1] / (timelistref[i][l]-timelistref[i][l-1]))
        dVdtref.append(tempdvdt)
    
    dVdtpulseC57 = []
    for i in range(len(voltlistpulse)):
        tempdvdt=[]
        for l in range(1, len(voltlistpulse[i])):
            tempdvdt.append(voltlistpulse[i][l]-voltlistpulse[i][l-1]/ (timelistpulse[i][l]-timelistpulse[i][l-1]))
        dVdtpulseC57.append(tempdvdt)

    
    dVdtrefC82 = []
    for i in range(len(voltlistC82)):
        tempdvdt=[]
        for l in range(1, len(voltlistC82[i])):
            tempdvdt.append(voltlistC82[i][l]-voltlistC82[i][l-1] / (timelistC82[i][l]-timelistC82[i][l-1]))
        dVdtrefC82.append(tempdvdt)


    dVdtpulsedouble = []
    for i in range(len(voltlistpulsedouble)):
        tempdvdt=[]
        for l in range(1, len(voltlistpulsedouble[i])):
            tempdvdt.append(voltlistpulsedouble[i][l]-voltlistpulsedouble[i][l-1]/ (timelistpulsedouble[i][l]-timelistpulsedouble[i][l-1]))
        dVdtpulsedouble.append(tempdvdt)

    
    dVdtpulseC5 = []
    for i in range(len(voltlistpulseC5)):
        tempdvdt=[]
        for l in range(1, len(voltlistpulseC5[i])):
            tempdvdt.append(voltlistpulseC5[i][l]-voltlistpulseC5[i][l-1]/ (timelistpulseC5[i][l]-timelistpulseC5[i][l-1]))
        dVdtpulseC5.append(tempdvdt)

    dVdtpulseC57CC20 = []
    for i in range(len(voltlistpulseC57CC20)):
        tempdvdt=[]
        for l in range(1, len(voltlistpulseC57CC20[i])):
            tempdvdt.append(voltlistpulseC57CC20[i][l]-voltlistpulseC57CC20[i][l-1] / (timelistpulseC57CC20[i][l]-timelistpulseC57CC20[i][l-1]))
        dVdtpulseC57CC20.append(tempdvdt)

    dVdtpulseC5s = []
    for i in range(len(voltlistpulseC575s)):
        tempdvdt=[]
        for l in range(1, len(voltlistpulseC575s[i])):
            tempdvdt.append(voltlistpulseC575s[i][l]-voltlistpulseC575s[i][l-1]/ (timelistpulseC575s[i][l]-timelistpulseC575s[i][l-1]))
        dVdtpulseC5s.append(tempdvdt)

    cell = 0
    numb_coeffdiff = 349
    poly_degreediff = 1
    ax2.plot(newtimeref[cell][1:], list(savgol_filter(dVdtref[cell], numb_coeffdiff, poly_degreediff)),  color = "r")
    ax2.plot(newtimepulse[cell][1:], list(savgol_filter(dVdtpulseC57[cell], numb_coeffdiff, poly_degreediff)), color = "cadetblue")
    ax2.plot(newtimepulsedouble[cell][1:], list(savgol_filter(dVdtpulsedouble[cell], numb_coeffdiff, poly_degreediff)),  color = "cornflowerblue")
    ax2.plot(newtimeC82[cell][1:], list(savgol_filter(dVdtrefC82[cell], numb_coeffdiff, poly_degreediff)),  color = "steelblue")
    
    ax2.plot(newtimepulseC5[cell][1:], list(savgol_filter(dVdtpulseC5[cell], numb_coeffdiff, poly_degreediff)),  color = "forestgreen")
    ax2.plot(newtimepulseC57CC20[cell][1:], list(savgol_filter(dVdtpulseC57CC20[cell], numb_coeffdiff, poly_degreediff)),  color = "darkgreen")
    ax2.plot(newtimepulseC575s[cell][1:], list(savgol_filter(dVdtpulseC5s[cell], numb_coeffdiff, poly_degreediff)),  color = "darkolivegreen")
    '''
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('SelfdischargeCoursedVdt5cells.png', dpi=800)
    os.chdir(origpath)
    plt.show()


def getEOLselfdischarge4erPlotrevirrevLoss():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        CCcapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, CCcapa

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCcapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, CCcapa

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa
    
    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        voltlists = []
        timelists = []
        cellsname = []
        currentlistscapa2 = []
        voltlistcapa2 = []
        selfdischarge = []
        selfdccapa= []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cellsname.append(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        sdcycle = 22
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        sdcycle = 23
                        capa2 = 6
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                currentliste = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Current(A)']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle + capa2)]["Step_"+str(sdcycle+ capa2)+'_Voltage(V)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                currentliste = data['split']["Step_"+str(sdcycle+ capa2)]['I']
                                voltlistCap2 = data['split']["Step_"+str(sdcycle+ capa2)]['V']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                currentlistscapa2.append(currentliste)
                                voltlistcapa2.append(voltlistCap2)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                #print(len(capacities))
                #print(capacities)
        return selfdischarge, selfdccapa, voltlists, timelists, currentlistscapa2, voltlistcapa2, cellsname, pulsecapa

    selfdischarge, selfdccapa, voltlistref, timelistref, currentlistC2ref, voltlistC2ref, cellnamesref, capalistCC20 = getCCref()

    selfdischargepulse, selfdccapapulse, voltlistpulse, timelistpulse, currentlistC2pulse, voltlistC2pulse, cellnamespulse, capalistC57 = getpulseform()

    selfdischargeC82, selfdccapaC82, voltlistC82, timelistC82, currentlistCC2C82, voltlistCC2C82, cellnamesC82, capalistC82 = getCCrefC83()

    selfdischargepulseC5, selfdccapapulseC5, voltlistpulseC5, timelistpulseC5, currentlistC2pulseC5, voltlistC2pulseC5, cellnamespulseC5, capalistC5 = getpulseformC5()

    selfdischargepulseC575s, selfdccapapulseC575s, voltlistpulseC575s, timelistpulseC575s, currentlistC2pulseC575s, voltlistC2pulseC575s, cellnamespulseC575s, capalistC575s = getpulseformC575s()

    selfdischargepulsedouble, selfdccapapulsedouble, voltlistpulsedouble, timelistpulsedouble, currentlistC2pulsedouble, voltlistC2pulsedouble, cellnamespulsedouble, capalistdouble = getpulseformdoublepulse()

    selfdischargepulseC57CC20, selfdccapapulseC57CC20, voltlistpulseC57CC20, timelistpulseC57CC20, currentlistC2pulseC57CC20, voltlistC2pulseC57CC20, cellnamespulseC57CC20, capalistC57CC20 = getpulseformC57CC20()


    # a) Voltage difference
    # b) reversible loss
    # c) irrev loss
    # d) self-dischagre current

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(2, 2)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)
    fig.tight_layout(w_pad=1.0)
    fontsizelabelaxis = 4

    # a) Voltage difference

    def voltage_diff(timelist, voltlist):
        fittimepulse = []
        fitvoltpulse = []
        selfdcvoltdiff = []
        for p in range(len(timelist)):
            temptime = []
            tempvolt = []
            for time in range(len(timelist[p])):
                if timelist[p][time]/(60*60) > 24 :
                    temptime.append(np.sqrt(timelist[p][time]/(60*60)))
                    tempvolt.append(voltlist[p][time])
                #timelist[p][time] = timelist[p][time]/(60*60)
            fittimepulse.append(temptime)
            fitvoltpulse.append(tempvolt)
            selfdcvoltdiff.append(tempvolt[-1] - tempvolt[0])
        return selfdcvoltdiff

    selfdcvoltdiffref = voltage_diff(timelistref, voltlistref)
    selfdcvoltdiffC82 = voltage_diff(timelistC82, voltlistC82)

    selfdcvoltdiffpulse = voltage_diff(timelistpulse, voltlistpulse)
    selfdcvoltdiffpulseC5 = voltage_diff(timelistpulseC5, voltlistpulseC5)
    selfdcvoltdiffpulseC575s = voltage_diff(timelistpulseC575s, voltlistpulseC575s)
    selfdcvoltdiffpulsedouble = voltage_diff(timelistpulsedouble, voltlistpulsedouble)
    selfdcvoltdiffpulseC57CC20= voltage_diff(timelistpulseC57CC20, voltlistpulseC57CC20)

    for i in range(len(cellnamespulse)):
        print(cellnamespulse[i])
        print(selfdcvoltdiffpulse[i])
    #    print(abs((capalistC57[i][2]+capalistC57[i][3])/capalistC57[i][0]))
    
    medianwidth = 0.5
    markersize = 2
    marker = 'D'
    boxlinewidth = 0.5
    
    CCref=ax[0,0].boxplot(selfdcvoltdiffref, positions = [0], patch_artist=True, boxprops=dict(facecolor='r', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    
    pulsedoublepuls = ax[0,0].boxplot(selfdcvoltdiffpulsedouble, positions = [5], patch_artist=True, boxprops=dict(facecolor='cadetblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    CCrefC83=ax[0,0].boxplot(selfdcvoltdiffC82, positions = [2], patch_artist=True, boxprops=dict(facecolor='steelblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulse = ax[0,0].boxplot(selfdcvoltdiffpulse, positions = [1], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    
    pulseC57CC20 = ax[0,0].boxplot(selfdcvoltdiffpulseC57CC20, positions = [4], patch_artist=True, boxprops=dict(facecolor='darkgreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulseC5s = ax[0,0].boxplot(selfdcvoltdiffpulseC575s, positions = [3], patch_artist=True, boxprops=dict(facecolor='darkolivegreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    pulseC5 = ax[0,0].boxplot(selfdcvoltdiffpulseC5, positions = [6], patch_artist=True, boxprops=dict(facecolor='forestgreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    
    #legendcolorsboxplot = [CCref["boxes"][0], pulse["boxes"][0], CCrefC83["boxes"][0],
    #                    pulseC5s["boxes"][0], pulseC57CC20["boxes"][0], pulsedoublepuls["boxes"][0],
    #                    pulseC5["boxes"][0]]
    
    #legendtext = ['N(CC20) = 8',  'N(PulsC5.7) = 29','N(CC8.23) = 10',
    #                   'N(PulsC5.7 5s) = 8', 'N(PulsC5.7 CC20) = 10',
    #                    'N(PulsC5.7double) = 17', 'N(PulsC5) = 12']
    #legend = ax[0,0].legend(legendcolorsboxplot, legendtext, fontsize = 2.5, bbox_to_anchor=(0.54, 0.55),
    #             borderaxespad=0, ncol = 1 )

    ax[0,0].set_ylabel('\u0394 Voltage [V]', fontsize=fontsizelabelaxis)
    
    ax[0,0].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,0].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,0].xaxis.set_ticklabels([])

    #b) irrev
    # Calculate difference of DC after self-discharge and CC2 
    revcapalossCC20 = []
    for i in range(len(capalistCC20)):
        if len(capalistCC20[i]) == 37:
            revcapalossCC20.append(capalistCC20[i][27]- capalistCC20[i][22])
        else:
            revcapalossCC20.append(capalistCC20[i][28]- capalistCC20[i][23])
    
    revcapalossCC82 = []
    for i in range(len(capalistC82)):
        if len(capalistC82[i]) == 37:
            revcapalossCC82.append(capalistC82[i][27]- capalistC82[i][22])
        else:
            revcapalossCC82.append(capalistC82[i][28]- capalistC82[i][23])
    
    revcapalossC57 = []
    for i in range(len(capalistC57)):
        if len(capalistC57[i]) == 37:
            revcapalossC57.append(capalistC57[i][27]- capalistC57[i][22])
        else:
            revcapalossC57.append(capalistC57[i][28]- capalistC57[i][23])
    
    revcapalossC5 = []
    for i in range(len(capalistC5)):
        if len(capalistC5[i]) == 37:
            revcapalossC5.append(capalistC5[i][27]- capalistC5[i][22])
        else:
            revcapalossC5.append(capalistC5[i][28]- capalistC5[i][23])

    revcapalossC575s = []
    for i in range(len(capalistC575s)):
        if len(capalistC575s[i]) == 37:
            revcapalossC575s.append(capalistC575s[i][27]- capalistC575s[i][22])
        else:
            revcapalossC575s.append(capalistC575s[i][28]- capalistC575s[i][23])

    revcapalossdouble = []
    for i in range(len(capalistdouble)):
        if len(capalistdouble[i]) == 37:
            revcapalossdouble.append(capalistdouble[i][27]- capalistdouble[i][22])
        else:
            revcapalossdouble.append(capalistdouble[i][28]- capalistdouble[i][23])

    revcapalossC57CC20 = []
    for i in range(len(capalistC57CC20)):
        if len(capalistC57CC20[i]) == 37:
            revcapalossC57CC20.append(capalistC57CC20[i][27]- capalistC57CC20[i][22])
        else:
            revcapalossC57CC20.append(capalistC57CC20[i][28]- capalistC57CC20[i][23])

    vioCC20 = ax[0,1].violinplot(revcapalossCC20, positions=[0],showmeans=False, showextrema=True, showmedians=True)
    vioC57 = ax[0,1].violinplot(revcapalossC57, positions=[1],showmeans=False, showextrema=True, showmedians=True)
    vioC82 = ax[0,1].violinplot(revcapalossCC82, positions=[2],showmeans=False, showextrema=True, showmedians=True)
    vioC575s = ax[0,1].violinplot(revcapalossC575s, positions=[3],showmeans=False, showextrema=True, showmedians=True)
    vioC57CC20 = ax[0,1].violinplot(revcapalossC57CC20, positions=[4],showmeans=False, showextrema=True, showmedians=True)
    viodouble = ax[0,1].violinplot(revcapalossdouble, positions=[5],showmeans=False, showextrema=True, showmedians=True)
    vioC5 = ax[0,1].violinplot(revcapalossC5, positions=[6],showmeans=False, showextrema=True, showmedians=True)

    for settings in vioCC20['bodies']:
        settings.set_facecolor('r')
        settings.set_alpha(1)
    for settings in vioC57['bodies']:
        settings.set_facecolor('b')
        settings.set_alpha(1)
    for settings in vioC82['bodies']:
        settings.set_facecolor('darkred')
        settings.set_alpha(1)
    for settings in vioC575s['bodies']:
        settings.set_facecolor('goldenrod')
        settings.set_alpha(1)
    for settings in vioC57CC20['bodies']:
        settings.set_facecolor('dimgrey')
        settings.set_alpha(1)
    for settings in viodouble['bodies']:
        settings.set_facecolor('darkgreen')
        settings.set_alpha(1)
    for settings in vioC5['bodies']:
        settings.set_facecolor('cornflowerblue')
        settings.set_alpha(1)

    templist = [vioCC20, vioC57, vioC82, vioC575s, vioC57CC20, viodouble, vioC5]
    for i in range(len(templist)):
        try:
            for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
                vp = templist[i][partname]
                vp.set_edgecolor("k")
                vp.set_linewidth(0.5)
            vp.set_edgecolor('lime')
        except:
            pass

    #legendcolors = [vioCC20["bodies"][0], vioC57["bodies"][0], vioC82["bodies"][0],
    #                    vioC575s["bodies"][0], vioC57CC20["bodies"][0], viodouble["bodies"][0],
    #                    vioC5["bodies"][0]]
    #legendtext = ['N(CC20) = '+str(len(revcapalossCC20)),  'N(PulsC5.7) = '+str(len(revcapalossC57)),'N(CC8.23) = '+str(len(revcapalossCC82)),
    #                   'N(PulsC5.7 5s) = '+str(len(revcapalossC575s)), 'N(PulsC5.7 CC20) = '+str(len(revcapalossC57CC20)),
    #                    'N(PulsC5.7double) = '+str(len(revcapalossdouble)), 'N(PulsC5) = '+str(len(revcapalossC5))]
    #legend = ax[0,1].legend(legendcolors, legendtext, fontsize = 2, bbox_to_anchor=(0.75, 0.95),
    #             borderaxespad=0, ncol =2 )
    
    ax[0,1].set_ylabel('reversible Loss [mAh]', fontsize=fontsizelabelaxis)
    
    ax[0,1].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,1].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,1].xaxis.set_ticklabels([])

    # c) irreversible capa loss CC1 - SD - rev Loss

    irrrevcapalossCC20 = []
    for i in range(len(capalistCC20)):
        if len(capalistCC20[i]) == 37:
                irrrevcapalossCC20.append(capalistCC20[i][27]- capalistCC20[i][8])#-revcapalossCC20[i])
        else:
                irrrevcapalossCC20.append(capalistCC20[i][28]- capalistCC20[i][9])#-revcapalossCC20[i])
    
    irrrevcapalossCC82 = []
    for i in range(len(capalistC82)):
        if len(capalistC82[i]) == 37:
                irrrevcapalossCC82.append(capalistC82[i][27]- capalistC82[i][8])#-revcapalossCC82[i])
        else:
                irrrevcapalossCC82.append(capalistC82[i][28]- capalistC82[i][9])#-revcapalossCC82[i])
           
    
    irrrevcapalossC57 = []
    for i in range(len(capalistC57)):
        if len(capalistC57[i]) == 37:
                irrrevcapalossC57.append(capalistC57[i][27]- capalistC57[i][8])#-revcapalossC57[i])
        else:
                irrrevcapalossC57.append(capalistC57[i][28]- capalistC57[i][9])#-revcapalossC57[i])
    
    irrrevcapalossC5 = []
    for i in range(len(capalistC5)):
        if len(capalistC5[i]) == 37:
                irrrevcapalossC5.append(capalistC5[i][27]- capalistC5[i][8])#-revcapalossC5[i])
        else:
                irrrevcapalossC5.append(capalistC5[i][28]- capalistC5[i][9])#-revcapalossC5[i])

    irrrevcapalossC575s = []
    for i in range(len(capalistC575s)):
        if len(capalistC575s[i]) == 37:
            if capalistC575s[i][8] == 0:
                irrev = capalistC575s[i][27]- capalistC575s[i][9]#-revcapalossC57[i]
            else:
                irrev = capalistC575s[i][27]- capalistC575s[i][8]#-revcapalossC57[i]

            irrrevcapalossC575s.append(irrev)
        else:
                irrrevcapalossC575s.append(capalistC575s[i][28]- capalistC575s[i][9])#-revcapalossC575s[i])

    irrrevcapalossdouble = []
    for i in range(len(capalistdouble)):
        if len(capalistdouble[i]) == 37:
                irrrevcapalossdouble.append(capalistdouble[i][27]- capalistdouble[i][8])#-revcapalossdouble[i])
        else:
                irrrevcapalossdouble.append(capalistdouble[i][28]- capalistdouble[i][9])#-revcapalossdouble[i])

    irrrevcapalossC57CC20 = []
    for i in range(len(capalistC57CC20)):
        if len(capalistC57CC20[i]) == 37:
                irrrevcapalossC57CC20.append(capalistC57CC20[i][27]- capalistC57CC20[i][8])#-revcapalossC57CC20[i])
        else:
                irrrevcapalossC57CC20.append(capalistC57CC20[i][28]- capalistC57CC20[i][9])#-revcapalossC57CC20[i])

    if len(irrrevcapalossCC20)>1:
        vioCC20irr = ax[1,0].violinplot(irrrevcapalossCC20, positions=[0],showmeans=False, showextrema=True, showmedians=True)
    else:
        vioCC20irr = ax[1,0].scatter(list([0]*len(irrrevcapalossCC20)), irrrevcapalossCC20, color = 'r')
    
    if len(irrrevcapalossC57)>1:
        vioC57irr = ax[1,0].violinplot(irrrevcapalossC57, positions=[1], showmeans=False, showextrema=True, showmedians=True)
    else:
        vioC57irr = ax[1,0].scatter(list([1]*len(irrrevcapalossC57)),irrrevcapalossC57, color = 'b')
    
    if len(irrrevcapalossCC82)>1:
        vioC82irr = ax[1,0].violinplot(irrrevcapalossCC82, positions=[2], showmeans=False, showextrema=True, showmedians=True)
    else:
        vioC82irr = ax[1,0].scatter(list([2]*len(irrrevcapalossCC82)),irrrevcapalossCC82, color = 'darkred')

    if len(irrrevcapalossC575s)>1:
        vioC575sirr = ax[1,0].violinplot(irrrevcapalossC575s, positions=[3], showmeans=False, showextrema=True, showmedians=True)
    else:
        vioC575sirr = ax[1,0].scatter(list([3]*len(irrrevcapalossC575s)),irrrevcapalossC575s, color = 'goldenrod')

    if len(irrrevcapalossC57CC20)>1:
        vioC57CC20irr = ax[1,0].violinplot(irrrevcapalossC57CC20, positions=[4], showmeans=False, showextrema=True, showmedians=True)
    else:
        vioC57CC20irr = ax[1,0].scatter(list([4]*len(irrrevcapalossC57CC20)),irrrevcapalossC57CC20, color = 'dimgrey')

    if len(irrrevcapalossdouble)>1:
        viodoubleirr = ax[1,0].violinplot(irrrevcapalossdouble, positions=[5], showmeans=False, showextrema=True, showmedians=True)
    else:
        viodoubleirr = ax[1,0].scatter(list([5]*len(irrrevcapalossdouble)),irrrevcapalossdouble, color = 'darkgreen')
    
    if len(irrrevcapalossC5)>1:
        vioC5irr = ax[1,0].violinplot(irrrevcapalossC5, positions=[6], showmeans=False, showextrema=True, showmedians=True)
    else:
        vioC5irr = ax[1,0].scatter(list([6]*len(irrrevcapalossC5)),irrrevcapalossC5, color = 'cornflowerblue')

    for settings in vioCC20irr['bodies']:
        settings.set_facecolor('r')
        settings.set_alpha(1)
    for settings in vioC57irr['bodies']:
        settings.set_facecolor('b')
        settings.set_alpha(1)
    for settings in vioC82irr['bodies']:
        settings.set_facecolor('darkred')
        settings.set_alpha(1)
    for settings in vioC575sirr['bodies']:
        settings.set_facecolor('goldenrod')
        settings.set_alpha(1)
    for settings in vioC57CC20irr['bodies']:
        settings.set_facecolor('dimgrey')
        settings.set_alpha(1)
    for settings in viodoubleirr['bodies']:
        settings.set_facecolor('darkgreen')
        settings.set_alpha(1)
    for settings in vioC5irr['bodies']:
        settings.set_facecolor('cornflowerblue')
        settings.set_alpha(1)

    templist = [vioCC20irr, vioC57irr, vioC82irr, vioC575sirr, vioC57CC20irr, viodoubleirr, vioC5irr]
    for i in range(len(templist)):
        try:
            for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
                vp = templist[i][partname]
                vp.set_edgecolor("k")
                vp.set_linewidth(0.5)
            vp.set_edgecolor('lime')
        except:
            pass

    legendcolors = [vioCC20["bodies"][0], vioC57["bodies"][0], vioC82["bodies"][0],
                        vioC575s["bodies"][0], vioC57CC20["bodies"][0], viodouble["bodies"][0],
                        vioC5["bodies"][0]]
    legendtext = ['N(CC20) = ' + str(len(vioCC20irr)),  'N(PulsC5.7) = ' + str(len(vioC57irr)),'N(CC8.23) = ' + str(len(vioC82irr)),
                       'N(PulsC5.7 5s) = ' + str(len(vioC575sirr)), 'N(PulsC5.7 CC20) = ' + str(len(vioC57CC20irr)),
                        'N(PulsC5.7double) = ' + str(len(viodoubleirr)), 'N(PulsC5) = ' + str(len(vioC5irr))]
    legend = ax[1,0].legend(legendcolors, legendtext, fontsize = 3.5, bbox_to_anchor=(2.35, -0.06),
                             frameon=True, borderaxespad=0, ncol = 4 )
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    ax[1,0].set_ylabel('irreversible Loss [mAh]', fontsize=fontsizelabelaxis)
    
    ax[1,0].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,0].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,0].xaxis.set_ticklabels([])

    # d) self-dischagre current of irreversible + reversible capacity loss / storage time

    irrcurrentCC20 = []
    for i in range(len(irrrevcapalossCC20)):
        irrcurrentCC20.append((irrrevcapalossCC20[i]+ revcapalossCC20[i])/(5*24)*1000)   
    
    irrcurrentC82 = []
    for i in range(len(irrrevcapalossCC82)):
        irrcurrentC82.append((irrrevcapalossCC82[i]+ revcapalossCC82[i])/(5*24)*1000) 

    irrcurrentC57 = []
    for i in range(len(irrrevcapalossC57)):
        irrcurrentC57.append((irrrevcapalossC57[i]+ irrrevcapalossC57[i])/(5*24)*1000) 
    
    irrcurrentC575s = []
    for i in range(len(irrrevcapalossC575s)):
        irrcurrentC575s.append((irrrevcapalossC575s[i]+ irrrevcapalossC575s[i])/(5*24)*1000) 
    
    irrcurrentdouble = []
    for i in range(len(irrrevcapalossdouble)):
        irrcurrentdouble.append((irrrevcapalossdouble[i]+ irrrevcapalossdouble[i])/(5*24)*1000) 
    
    irrcurrentC5 = []
    for i in range(len(irrrevcapalossC5)):
        irrcurrentC5.append((irrrevcapalossC5[i]+ irrrevcapalossC5[i])/(5*24)*1000) 
    
    irrcurrentC57CC20 = []
    for i in range(len(irrrevcapalossC57CC20)):
        irrcurrentC57CC20.append((irrrevcapalossC57CC20[i]+ irrrevcapalossC57CC20[i])/(5*24)*1000) 



    if len(irrcurrentCC20)>1:
        CCrefirr =ax[1,1].boxplot(irrcurrentCC20, positions = [0], patch_artist=True, boxprops=dict(facecolor='r', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    else:
        CCrefirr =ax[1,1].scatter(list([0]*len(irrrevcapalossC5)),irrcurrentCC20, color = 'r')

    if len(irrcurrentC57)>1:
        pulseirr = ax[1,1].boxplot(irrcurrentC57, positions = [1], patch_artist=True, boxprops=dict(facecolor='b', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    else:
        pulseirr =ax[1,1].scatter(list([1]*len(irrcurrentC57)),irrcurrentC57, color = 'b')
    
    if len(irrcurrentC82)>1:
        CCrefC83irr =ax[1,1].boxplot(irrcurrentC82, positions = [2], patch_artist=True, boxprops=dict(facecolor='darkred', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    else:
        CCrefC83irr =ax[1,1].scatter(list([2]*len(irrcurrentC82)),irrcurrentC82, color = 'darkred')

    if len(irrcurrentC575s)>1:
        pulseC5sirr = ax[1,1].boxplot(irrcurrentC575s, positions = [3], patch_artist=True, boxprops=dict(facecolor='goldenrod', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    else:
        pulseC5sirr =ax[1,1].scatter(list([3]*len(irrcurrentC575s)),irrcurrentC575s, color = 'goldenrod')

    if len(irrcurrentC57CC20)>1:
        pulseC57CC20irr = ax[1,1].boxplot(irrcurrentC57CC20, positions = [4], patch_artist=True, boxprops=dict(facecolor='dimgrey', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    else:
        pulseC57CC20irr =ax[1,1].scatter(list([4]*len(irrcurrentC57CC20)),irrcurrentC57CC20, color = 'dimgrey')
    
    if len(irrcurrentdouble)>1:
        pulsedoublepulsirr = ax[1,1].boxplot(irrcurrentdouble, positions = [5], patch_artist=True, boxprops=dict(facecolor='darkgreen', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    else:
        pulsedoublepulsirr =ax[1,1].scatter(list([5]*len(irrcurrentdouble)),irrcurrentdouble, color = 'darkgreen')
    
    if len(irrcurrentC5)>1:
        pulseC5irr = ax[1,1].boxplot(irrcurrentC5, positions = [6], patch_artist=True, boxprops=dict(facecolor='cornflowerblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    else:
        pulseC5irr =ax[1,1].scatter(list([6]*len(irrcurrentC5)),irrcurrentC5, color = 'cornflowerblue')
        
    #legendcolorsboxplot = [CCrefirr["boxes"][0], pulseirr["boxes"][0], CCrefC83irr["boxes"][0],
    #                    pulseC5sirr["boxes"][0], pulseC57CC20irr["boxes"][0], pulsedoublepulsirr["boxes"][0],
    #                    pulseC5irr["boxes"][0]]
    #legendtext = ['N(CC20) = '+str(len(irrcurrentCC20)),  'N(PulsC5.7) = '+str(len(irrcurrentC57)),'N(CC8.23) = '+str(len(irrcurrentC82)),
    #                   'N(PulsC5.7 5s) = '+str(len(irrcurrentC575s)), 'N(PulsC5.7 CC20) = '+str(len(irrcurrentC57CC20)),
    #                    'N(PulsC5.7double) = '+str(len(irrcurrentdouble)), 'N(PulsC5) = '+str(len(irrcurrentC5))]
    #legend = ax[1,1].legend(legendcolorsboxplot, legendtext, fontsize = 2.5, bbox_to_anchor=(0.45, 0.5),
    #             borderaxespad=0, ncol = 1 )

    ax[1,1].set_ylabel('self-discharge current [\u03BCA]', fontsize=fontsizelabelaxis)
    
    ax[1,1].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,1].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,1].xaxis.set_ticklabels([])

    xpos = -0.13
    ypos = 1.15
    ax[0, 0].text(xpos, ypos, 'a)', transform=ax[0, 0].transAxes, fontsize=6, fontweight='bold', va='top')
    ax[0, 1].text(xpos, ypos, 'b)', transform=ax[0, 1].transAxes, fontsize=6, fontweight='bold', va='top')
    ax[1, 0].text(xpos, ypos, 'c)', transform=ax[1, 0].transAxes, fontsize=6, fontweight='bold', va='top')
    ax[1, 1].text(xpos, ypos, 'd)', transform=ax[1, 1].transAxes, fontsize=6, fontweight='bold', va='top')

    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('SelfdischargeMultiPlot5cells.png', dpi=1000)
    os.chdir(origpath)
    plt.show()


# Not for SI but good to look at 
def getCurrentCVAfterForm():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        CCcapa = []
        CVcurrentlist = []
        CVtimelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                print(len(cycles))
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(data['split'][l]['t'],data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                    temptime.append(timelist)
                    tempcurrent.append(currentlist)
                CCcapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])

        return CCcapa, CVcurrentlist, CVtimelists

    def getCCrefC83():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCcapa = []
        CVcurrentlist = []
        CVtimelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)


        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent= []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])
                CCcapa.append(capacities)

        return CCcapa, CVcurrentlist, CVtimelists
    
    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        CVcurrentlist = []
        CVtimelists = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    temptime.append(timelist)
                    tempcurrent.append(currentlist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])

        return pulsecapa, CVcurrentlist, CVtimelists

    def getpulseformC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        CVcurrentlist = []
        CVtimelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)            
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities) 
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])    
        return pulsecapa, CVcurrentlist, CVtimelists

    def getpulseformC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        CVcurrentlist = []
        CVtimelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])
                
        return pulsecapa, CVcurrentlist, CVtimelists

    def getpulseformdoublepulse():
        hdf5path = pathdoublepulse
        hdf5path2 = pathdoublepulse2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        CVcurrentlist = []
        CVtimelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])

        os.chdir(hdf5path2)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])

        return pulsecapa, CVcurrentlist, CVtimelists

    def getpulseformC57CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        pulsecapa = []
        CVcurrentlist = []
        CVtimelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                tempcurrent = []
                temptime = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    tempcurrent.append(currentlist)
                    temptime.append(timelist)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                CVcurrentlist.append(tempcurrent[3])
                CVtimelists.append(temptime[3])
                
        return pulsecapa, CVcurrentlist, CVtimelists


    eolcapalist, cvcurrent, cvtime = getCCref()
    eolcapalistC823, cvcurrentC823, cvtimeC823 = getCCrefC83()
    eolcapalistpulse,  cvcurrentC57, cvtimeC57 = getpulseform()
    eolcapalistC5,  cvcurrentC5, cvtimeC5 = getpulseformC5()
    eolcapalistC575s, cvcurrentC575s, cvtimeC575s = getpulseformC575s()
    eolcapalistdoubleulse, cvcurrentdoublepulse, cvtimedoublepulse = getpulseformdoublepulse()
    eolcapalistC57CC20, cvcurrentC57CC20, cvtimeC57CC20 = getpulseformC57CC20()

    # cutoff value for CC20 and C57 is postiv 0.00167 and therefore it was charged for some time as well
    # cutoff vlaue for all other is 0 and nice

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(3,3)
    fig.delaxes(ax[2,1])
    fig.delaxes(ax[2,2])
    plt.style.use(['science','nature','no-latex'])


    for i in range(len(cvcurrent)):
        ax[0,0].plot(cvtime[i],cvcurrent[i], color = 'r')

    for i in range(len(cvcurrentC823)):
        ax[0,1].plot(cvtimeC823[i],cvcurrentC823[i], color = 'b')

    for i in range(len(cvcurrentC5)):
        ax[0,2].plot(cvtimeC5[i],cvcurrentC5[i], color = 'g')

    for i in range(len(cvcurrentC57)):
        ax[1,0].plot(cvtimeC57[i],cvcurrentC57[i], color = 'cornflowerblue')
    
    for i in range(len(cvcurrentC575s)):
        ax[1,1].plot(cvtimeC575s[i],cvcurrentC575s[i], color = 'darkred')
    
    for i in range(len(cvcurrentdoublepulse)):
        ax[1,2].plot(cvtimedoublepulse[i],cvcurrentdoublepulse[i], color = 'darkgreen')
    
    for i in range(len(cvcurrentC57CC20)):
        ax[2,0].plot(cvtimeC57CC20[i],cvcurrentC57CC20[i], color = 'dimgrey')

    ax[0,0].set_xlabel('time [h]', fontsize=8)
    ax[0,0].set_ylabel('Current [mAh]', fontsize=8)

    ax[0,0].set_ylim([-0.004,0.00025])

    origpath = os.getcwd()
    os.chdir(plotfolder)
    #fig.savefig('FormationDC-CV.png', dpi=800)
    os.chdir(origpath)
    plt.show()


# Maybe SI
def getFormCapavsVolt():

    def getCCrefCapaVolt():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        CCcapa = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        refvolt = []
        refcapa = []
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                c=[]
                for p in range(len(currentlist)):
                    if p == 0:
                        cap = 0
                    else:
                        cap = c[-1]+ currentlist[p]*(timelist[p]-timelist[p-1])* (1000/3600)
                    c.append(abs(cap))
                #print(content[i])
                #plt.plot(timelist,data['split'][l]['V'])
                #plt.show()
                #print(c)
                refcapa.append(c)
                refvolt.append(voltlist)            
        return refcapa, refvolt
    
    def getpulseCapaVolt():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)

        pulsevolt = []
        pulsecapa = []
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                if currentlist[0] != currentlist[1]:
                    del currentlist[0]
                    del timelist[0]
                if currentlist[-1] != currentlist[-2]:
                    del currentlist[-1]
                    del timelist[-1]
                c=[]
                for p in range(len(currentlist)):
                    if p == 0:
                        cap =0
                    else:
                        cap = c[-1]+ currentlist[p]*(timelist[p]-timelist[p-1])* (1000/3600)
                    c.append(abs(cap))
                #print(content[i])
                #plt.plot(timelist,data['split'][l]['V'])
                #plt.show()
                #print(c)
                pulsecapa.append(c)
                pulsevolt.append(voltlist)            
        return pulsecapa, pulsevolt

    refcapa, refvolt = getCCrefCapaVolt()

    pulsecapa, pulsevolt = getpulseCapaVolt()

    for i in range(len(refvolt)):
        plt.plot(refvolt[i],refcapa[i], color= 'r')
    for i in range(len(pulsevolt)):
        plt.plot(pulsevolt[i], pulsecapa[i], color = 'b')

    plt.style.use(['science','no-latex'])
    plt.xlabel('Voltage (V)')
    plt.ylabel('Capacity (mAh)')
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('FormationCapavsVolt.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def getEOLselfdischargeCapaVergleich():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        selfdischarge = []
        selfdccapa = []
        Capa1 = []
        Capa2 = []
        Capa3 = []
        voltlists = []
        timelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)

        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        Capa3.append(abs(capacities[32]))
                        Capa2.append(abs(capacities[27]))
                        Capa1.append(abs(capacities[8])) 
                        sdcycle = 22
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        Capa3.append(abs(capacities[33]))
                        Capa2.append(abs(capacities[28]))
                        Capa1.append(abs(capacities[9]))
                        sdcycle = 23
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                voltlists.append(voltlist)
                                timelists.append(timeliste)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
        return selfdischarge, selfdccapa, Capa1, Capa2, Capa3 , voltlists, timelists

    def getpulseform():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        selfdischargepulse = []
        selfdccapapulse = []
        pulseCapa1 = []
        pulseCapa2 = []
        pulseCapa3 = []
        voltlistspulse = []
        timelistspulse = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                #print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)

                if len(capacities) == 37:
                    try:
                        selfdccapapulse.append(abs(capacities[22]))
                        pulseCapa3.append(abs(capacities[32]))
                        pulseCapa2.append(abs(capacities[27]))
                        pulseCapa1.append(abs(capacities[8]))
                        sdcycle = 22
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0])) 
                    except:
                        print(i)
                else:
                    try:
                        selfdccapapulse.append(abs(capacities[23]))
                        pulseCapa3.append(abs(capacities[33]))
                        pulseCapa2.append(abs(capacities[28]))
                        pulseCapa1.append(abs(capacities[9]))
                        sdcycle = 23
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0])) 
                    except:
                        print(i)
        return selfdischargepulse, selfdccapapulse, pulseCapa1, pulseCapa2, pulseCapa3, voltlistspulse, timelistspulse


    selfdischarge, selfdccapa, Capa1, Capa2, Capa3, voltlistref, timelistref = getCCref()

    selfdischargepulse, selfdccapapulse, pulseCapa1, pulseCapa2, pulseCapa3, voltlistpulse, timelistpulse = getpulseform()

    CC = [selfdccapa ,selfdischarge ]
    pulsedata = [ selfdccapapulse, selfdischargepulse]

    """
    for p in range(len(voltlistref)):
        plt.scatter(timelistref[p],voltlistref[p], color = 'r')

    for p in range(len(voltlistpulse)):
        plt.scatter(timelistpulse[p],voltlistpulse[p], color = 'b')

    plt.show()
    """
    #ax[y[o]].set_xticks([])  
    #ax[y[o]].title.set_text(plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1, figsize=(4, 2), sharex=True)
    plt.style.use(['science','nature','no-latex'])
    #fig.suptitle('Coulombic efficiency', fontsize=14)

    boxwidth = 0.15
    colorcapa1 = 'b'
    colorcapa2 = 'g'
    colorcapa1pulse = 'b'
    colorcapa2pulse = 'g'
    whiskerlinewidth = 3
    CCref = ax.boxplot(CC[0], positions = [0], patch_artist=True, boxprops=dict(facecolor='r'))
    for element in CCref.keys():
        plt.setp(CCref[element], visible=False)
    boxcapa1 = ax.boxplot(Capa1, positions = [1],  patch_artist=True, widths=(boxwidth+0.05), boxprops=dict(facecolor=colorcapa1, alpha = 0.5))
    boxcapa2 = ax.boxplot(Capa2, positions = [1], patch_artist=True, widths=(boxwidth+0.1), boxprops=dict(facecolor=colorcapa2, alpha = 0.5))
    CCref = ax.boxplot(CC[0], positions = [1], patch_artist=True,widths=(boxwidth), boxprops=dict(facecolor='r', linewidth=1.5))

    boxcapa1pulse = ax.boxplot(pulseCapa1, positions = [1.5], patch_artist=True, widths=(boxwidth+0.05), boxprops=dict(facecolor=colorcapa1pulse, alpha = 0.5))
    boxcapa2pulse = ax.boxplot(pulseCapa2, positions = [1.5], patch_artist=True, widths=(boxwidth+0.1), boxprops=dict(facecolor=colorcapa2pulse, alpha = 0.5))
    pulse = ax.boxplot(pulsedata[0], positions = [1.5], patch_artist=True,widths=(boxwidth), boxprops=dict(facecolor='r', linewidth=1.5)) 
    ax.set_xlim([0.75,1.9])
    ax.legend([boxcapa1["boxes"][0], boxcapa2["boxes"][0], CCref["boxes"][0]], ['DC Capacity-Check-1', 'DC Capacity-Check-2','DC after self-discharge' ], loc='upper right')
    #ax[x[o],y[o]].set_ylimtitles[o])
    for median in CCref['medians']:
        median.set( linewidth = 1.5, color = 'k',)
    for median in pulse['medians']:
        median.set( linewidth = 1.5, color = 'k')
    for median in boxcapa1['medians']:
        median.set( linewidth = 1, color = colorcapa1, alpha = 0.7)
    for median in boxcapa2['medians']:
        median.set( linewidth = 1, color = colorcapa2, alpha = 0.7)
    for median in boxcapa1pulse['medians']:
        median.set( linewidth = 1, color = colorcapa1pulse, alpha = 0.7)
    for median in boxcapa2pulse['medians']:
        median.set( linewidth = 1, color = colorcapa2pulse, alpha = 0.7)

    for cap in CCref['caps']:
        cap.set(linewidth = 2)
    
    for flier in CCref['fliers']:
        flier.set(
            markerfacecolor  ="r",
                )
        
    for cap in boxcapa1['caps']:
        cap.set(color = colorcapa1,
                linewidth = 2)
    for flier in boxcapa1['fliers']:
        flier.set(
            markerfacecolor  =colorcapa1,
                )
    for whisker in boxcapa1['whiskers']:
        whisker.set(color =colorcapa1,
                linewidth = whiskerlinewidth,
                linestyle =":")
        
    for cap in boxcapa2['caps']:
        cap.set(color = colorcapa2,
                linewidth = 2)
    for flier in boxcapa2['fliers']:
        flier.set(
            markerfacecolor  =colorcapa2,
                )
    for whisker in boxcapa2['whiskers']:
        whisker.set(color =colorcapa2,
                linewidth = whiskerlinewidth,
                linestyle =":")
    
    for cap in pulse['caps']:
        cap.set(linewidth = 2)

    for flier in pulse['fliers']:
        flier.set(
            markerfacecolor  ="r",
                )
    for cap in boxcapa1pulse['caps']:
        cap.set(color = colorcapa1pulse,
                linewidth = 2)
    for flier in boxcapa1pulse['fliers']:
        flier.set(
            markerfacecolor  =colorcapa1pulse,
                )
    for whisker in boxcapa1pulse['whiskers']:
        whisker.set(color =colorcapa1pulse,
                linewidth = whiskerlinewidth,
                linestyle =":")
        
    for cap in boxcapa2pulse['caps']:
        cap.set(color = colorcapa2pulse,
                linewidth = 2)
    for flier in boxcapa2pulse['fliers']:
        flier.set(
            markerfacecolor  =colorcapa2pulse,
                )
    for whisker in boxcapa2pulse['whiskers']:
        whisker.set(color =colorcapa2pulse,
                linewidth = whiskerlinewidth,
                linestyle =":")

    ax.set_ylabel('Discharge capacity [mAh]', fontsize=8)
    ax.set_xticks((1,1.5))
    ax.set_xticklabels(['CC', 'Pulse'], color = 'k', fontsize = 8)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('SelfdischargeCapavsChecks.png', dpi=800)
    os.chdir(origpath)
    plt.show()


# vielleicht SI als 7er Subplot 
def getEOLselfdischargeVoltageCourse():

    def getCCref():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        selfdischarge = []
        selfdccapa = []
        Capa1 = []
        Capa2 = []
        Capa3 = []
        voltlists = []
        timelists = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)

        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]['V']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    if len(cycles) == 37 and l == "Step_22":
                        voltlists.append(voltlist)
                        timelists.append(timelist)
                    elif len(cycles) == 38 and l == "Step_23":
                        voltlists.append(voltlist)
                        timelists.append(timelist)
                    else:
                        pass
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                if len(capacities) == 37:
                    try:
                        selfdccapa.append(abs(capacities[22]))
                        Capa3.append(abs(capacities[32]))
                        Capa2.append(abs(capacities[27]))
                        Capa1.append(abs(capacities[8])) 
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
                else:
                    try:
                        selfdccapa.append(abs(capacities[23]))
                        Capa3.append(abs(capacities[33]))
                        Capa2.append(abs(capacities[28]))
                        Capa1.append(abs(capacities[9]))
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0]))
                    except:
                        print(i)
        return selfdischarge, selfdccapa, Capa1, Capa2, Capa3 , voltlists, timelists

    def getpulseform():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        selfdischargepulse = []
        selfdccapapulse = []
        pulseCapa1 = []
        pulseCapa2 = []
        pulseCapa3 = []
        voltlistspulse = []
        timelistspulse = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                #print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)

                if len(capacities) == 37:
                    try:
                        selfdccapapulse.append(abs(capacities[22]))
                        pulseCapa3.append(abs(capacities[32]))
                        pulseCapa2.append(abs(capacities[27]))
                        pulseCapa1.append(abs(capacities[8]))
                        sdcycle = 22
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0])) 
                    except:
                        print(i)
                else:
                    try:
                        selfdccapapulse.append(abs(capacities[23]))
                        pulseCapa3.append(abs(capacities[33]))
                        pulseCapa2.append(abs(capacities[28]))
                        pulseCapa1.append(abs(capacities[9]))
                        sdcycle = 23
                        try:
                                voltlist = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Voltage(V)']
                                timeliste = data['split']["Step_"+str(sdcycle)]["Step_"+str(sdcycle)+'_Testtime(s)']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                        except:
                                voltlist = data['split']["Step_"+str(sdcycle)]['V']
                                timeliste = data['split']["Step_"+str(sdcycle)]['t']
                                voltlistspulse.append(voltlist)
                                timelistspulse.append(timeliste)
                                #plt.plot(timeliste, data['split'][l]['V'])
                                #plt.show()
                        selfdischarge.append(abs(voltlist[-1]- voltlist[0])) 
                    except:
                        print(i)
        return selfdischargepulse, selfdccapapulse, pulseCapa1, pulseCapa2, pulseCapa3, voltlistspulse, timelistspulse


    selfdischarge, selfdccapa, Capa1, Capa2, Capa3, voltlistref, timelistref = getCCref()

    selfdischargepulse, selfdccapapulse, pulseCapa1, pulseCapa2, pulseCapa3, voltlistpulse, timelistpulse = getpulseform()

    CC = [selfdccapa ,selfdischarge ]
    pulsedata = [ selfdccapapulse, selfdischargepulse]

    for p in range(len(voltlistref)):
        for time in range(len(timelistref[p])):
            timelistref[p][time] = timelistref[p][time]/(60*60*24)


    for p in range(len(voltlistpulse)):
        for time in range(len(timelistpulse[p])):
            timelistpulse[p][time] = timelistpulse[p][time]/(60*60*24)

    """
    linewidth = 2
    alpha = 0.1

    durchschnittspositionen = [np.mean(kurve) for kurve in voltlistref]

    # Index der Kurve mit der durchschnittlichsten Position
    mittigste_kurve_index = np.argmin(np.abs(np.array(durchschnittspositionen) - np.median(durchschnittspositionen)))

    # Index der Kurve mit der höchsten durchschnittlichen Position
    hoechste_kurve_index = np.argmax(durchschnittspositionen)

    # Index der Kurve mit der niedrigsten durchschnittlichen Position
    niedrigste_kurve_index = np.argmin(durchschnittspositionen)

    # Plot der mittigsten Kurve
    plt.plot(voltlistref[mittigste_kurve_index], color='darkred', label='Mittigste Kurve', linewidth = linewidth)
    plt.fill_between(range(len(voltlistref[hoechste_kurve_index])), voltlistref[niedrigste_kurve_index], voltlistref[hoechste_kurve_index], color='r', alpha=alpha)


    durchschnittspositionen = [np.mean(kurve) for kurve in voltlistpulse]

    # Index der Kurve mit der durchschnittlichsten Position
    mittigste_kurve_index = np.argmin(np.abs(np.array(durchschnittspositionen) - np.median(durchschnittspositionen)))

    # Index der Kurve mit der höchsten durchschnittlichen Position
    hoechste_kurve_index = np.argmax(durchschnittspositionen)

    # Index der Kurve mit der niedrigsten durchschnittlichen Position
    niedrigste_kurve_index = np.argmin(durchschnittspositionen)

    # Plot der mittigsten Kurve
    plt.plot(voltlistpulse[mittigste_kurve_index], color='cornflowerblue', linewidth = linewidth)
    plt.fill_between(range(len(voltlistpulse[hoechste_kurve_index])), voltlistpulse[niedrigste_kurve_index], voltlistpulse[hoechste_kurve_index], color='b', alpha=alpha)

    """
    #ax[y[o]].set_xticks([])  
    #ax[y[o]].title.set_text(plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1, 1, figsize=(4, 2), sharex=True)
    plt.style.use(['science','nature','no-latex'])    
    ax.axhspan(4.11, 4.188, xmin=0.020, xmax=0.21 , color='#008000', alpha=0.4) # red
    ax.text(0.03, 0.98, "Relaxation", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'k')
    ax.axhspan(4.11, 4.188, xmin=0.021, xmax=0.99 , color='m', alpha=0.2) # red
    ax.text(0.5, 0.98, "Self-discharge", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'k')
    
    ax.scatter(timelistref[0], voltlistref[0], label = "CC", color = 'r')
    ax.scatter(timelistpulse[0], voltlistpulse[0], label = "Pulse", color = 'b')

    plt.legend(loc='upper right', fontsize = 8)
    ax.set_xlim([-0.1,5.1])
    ax.set_ylim([4.11, 4.188])
    ax.set_xlabel('Time [d]', fontsize=8)
    ax.set_ylabel('Voltage [V]', fontsize=8)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    #plt.savefig('SelfdischargeVoltageCourse.png', dpi=800)
    os.chdir(origpath)
    plt.show()

# falls nötig backup

# langzeit auswertungen

pathpulse = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\PulsedC58\LongtermC58"
pathCCref = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\CCCV_C20\CCCV_InForm_300Cycles"
#pathCCref = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20Call\EOL C20 CC - eff C20 pulsed ref\300Cycles"
pathC5pulse = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\C5 pulse\Longterm'
pathC83CC = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\C823 CCCV\Longterm'
pathC57ref1 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 1 C58 5s pulse batch 2\300 Cycles'
pathC57ref2= r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 2 C58 CC20rest batch 1\300 Cycles'
pathdoublepulse = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL Double Pulse batch 1\300Cycles'
pathdoublepulse2 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL Double Pulse batch 2\300Cycles'
path5spulse = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 5s pulse batch 2\300Cycles'
pathC58CC20 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse CC20 rest batch 1\300Cycles'
plotfolder = r'C:\Users\Leon\Desktop\PhD\KIT\My papers\Formation Strategies'

def getCycleCapaplot5cellsrealtiv():


    def getCCCapaCycles():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'ZSWInform'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa

    def getpulseCapaCycles():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        
        hdf5path = pathC57ref1
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)

        hdf5path = pathC57ref2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa
    
    def getCC83CapaCycles():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa
            
    def getpulseCapaCyclesC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa

    def getpulseCapaCyclesC58CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa

    def getdoublepulseCapaCycles():
        hdf5path = pathdoublepulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        
        hdf5path = pathdoublepulse2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa
    
    def getpulseCapaCyclesC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa


    cccvdccapa, legendscccv, ccallcapa = getCCCapaCycles()
    pulsedccapa, legendspulse, ccallcapapulse = getpulseCapaCycles()
    cccvC83dccapa, legendscccvC83, cc823allcapa = getCC83CapaCycles()
    pulsedC585scapa, legendspulseC585s, ccallcapapulseC585s = getpulseCapaCyclesC575s()
    pulsedC58CC20ccapa, legendspulseC58CC20, ccallcapaC58CC20 = getpulseCapaCyclesC58CC20()
    doublepulsedcapa, doublelegendspulse, ccallcapadouble = getdoublepulseCapaCycles()
    pulsedC5ccapa, legendspulseC5, ccallcapaC5  = getpulseCapaCyclesC5()

    '''
    # rot: Zelle 1 am ehesten noch 2
    del ccallcapa[0]

    
    # blau hier bleiben folgende Zellen
    #  3,  6, 10, 11, 19

    pulsebackup = copy.deepcopy(ccallcapapulse)
    temp = []
    temp.append(ccallcapapulse[2])
    temp.append(ccallcapapulse[5])
    temp.append(ccallcapapulse[9])
    temp.append(ccallcapapulse[10])
    temp.append(ccallcapapulse[18])
    ccallcapapulse = temp
    

    # dunkelrot: 1, 3, 5, 6
    del cc823allcapa[0]
    del cc823allcapa[1]
    del cc823allcapa[2]
    del cc823allcapa[2]

    # gelb: 1, 3, 4
    del ccallcapapulseC585s[0]
    del ccallcapapulseC585s[1]
    del ccallcapapulseC585s[1]

    # grau: 2, 5, 7, 8, 10
    del ccallcapaC58CC20[1]
    del ccallcapaC58CC20[3]
    del ccallcapaC58CC20[4]
    del ccallcapaC58CC20[4]
    del ccallcapaC58CC20[5]

    # grün: 2, 3, 4, 5,  8, 10, 12, 13,  14, 15, 16
    del ccallcapadouble[1]
    del ccallcapadouble[1]
    del ccallcapadouble[1]
    del ccallcapadouble[1]
    del ccallcapadouble[3]
    del ccallcapadouble[4]
    del ccallcapadouble[5]
    del ccallcapadouble[5]
    del ccallcapadouble[5]
    del ccallcapadouble[5]
    del ccallcapadouble[5]

    # dunkelblau: 1, 3, 4, 5, 7
    del ccallcapaC5[0]
    del ccallcapaC5[1]
    del ccallcapaC5[1]
    del ccallcapaC5[1]
    del ccallcapaC5[2]
    '''

    allCapcities = [ccallcapa, ccallcapapulse,ccallcapadouble,  cc823allcapa,   ccallcapaC5 , ccallcapaC58CC20,ccallcapapulseC585s]
    tempallCapcities = copy.deepcopy(allCapcities)
    #allCapcities = copy.deepcopy(tempallCapcities)

    for i in range(len(allCapcities)):
        #print(len(allCapcities[i]))
        for o in range(len(allCapcities[i])):
            reference = copy.deepcopy(allCapcities[i][o][1])
            for l in range(len(allCapcities[i][o])):
                temp = allCapcities[i][o][l] / reference
                allCapcities[i][o][l] = temp

    fig, ax = plt.subplot_mosaic(
    [["top left", "top centre", "top right"],
        ["middle left", "middle centre", "middle right"],
     ["bottom row", "bottom row", "bottom row"]], sharey=True
    )
    plt.style.use(['science','nature','no-latex'])
    fig.tight_layout(h_pad=2.0)
    
    x = ["top left", "top centre", "top right",
        "middle left", "middle centre", "middle right",
        "bottom row"]
    labels = ['N(CC C/20) = ' + str(len(ccallcapa)),'N(Pulsed C/5.7) = ' + str(len(ccallcapapulse)), 'N(Pulsed C/5.7double) = ' + str(len(ccallcapadouble)), 'N(CC C/8.2) = ' + str(len(cc823allcapa)),
                        'N(Pulsed C/5) = ' + str(len(ccallcapaC5)), 'N(Mixed C/5.7 CC20) = ' + str(len(ccallcapaC58CC20)),'N(Pulsed C/5.7 5s) = ' + str(len(ccallcapapulseC585s))
                        ]
    colors = ['r', 'cornflowerblue','cadetblue', 'steelblue', 'forestgreen', 'darkgreen', 'darkolivegreen' ]
    for o in range(len(x)):
        for i in range(len(allCapcities[o])):
            cycles = np.linspace(1,len(allCapcities[o][i])-1,len(allCapcities[o][i])-1 )
            cccvref=ax[x[o]].scatter(cycles, allCapcities[o][i][1:], s = 0.05, label = labels[o], color = colors[o])
        line= ax[x[o]].plot(cycles, [0.8]*len(cycles), color= 'k', linewidth = 0.75)
        ax[x[o]].set_xlim([0,300])
        ax[x[o]].set_ylim([0.5,1.2])
        ax[x[o]].set_xticks([0,100,200,300])
        ax[x[o]].set_yticks([0.8,1,1.2])
        ax[x[o]].set_yticklabels([0.8,1,1.5],fontsize=4)
        ax[x[o]].set_xticklabels([0,100,200,300],fontsize=4)
        ax[x[o]].set_title(labels[o], fontdict={'color': colors[o]}, fontsize = 5)
    
    #ax["middle left"].set_ylabel('Discharge Capacity [mAh]', fontsize=6)
    ax["middle left"].set_ylabel('relative discharge capacity [%]', fontsize=6)
    ax["bottom row"].set_xlabel('cycle number', fontsize=6)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('DischargecapacityVSCyclenumberALL.png', dpi=800)
    os.chdir(origpath)
    plt.show()

def getCycleCapaplot5cellsABSOLUTE():

    # to include all data change path to C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20Call\PulsedC58\LongtermC58

    def getCCCapaCycles():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'ZSWInform'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa

    def getpulseCapaCycles():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        
        hdf5path = pathC57ref1
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)

        hdf5path = pathC57ref2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa
    
    def getCC83CapaCycles():
        hdf5path = pathC83CC
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa
            
    def getpulseCapaCyclesC575s():
        hdf5path = path5spulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa

    def getpulseCapaCyclesC58CC20():
        hdf5path = pathC58CC20
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa

    def getdoublepulseCapaCycles():
        hdf5path = pathdoublepulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        
        hdf5path = pathdoublepulse2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa
    
    def getpulseCapaCyclesC5():
        hdf5path = pathC5pulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test_'

        CCdccapa = [[],[],[],[]]
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                CCdccapa[0].append(capacities[0])
                CCdccapa[1].append(capacities[1])
                CCdccapa[2].append(capacities[2])
                CCdccapa[3].append(capacities[3])
                allCapa.append(capacities)
        return CCdccapa, legendsCC, allCapa


    cccvdccapa, legendscccv, ccallcapa = getCCCapaCycles()
    pulsedccapa, legendspulse, ccallcapapulse = getpulseCapaCycles()
    cccvC83dccapa, legendscccvC83, cc823allcapa = getCC83CapaCycles()
    pulsedC585scapa, legendspulseC585s, ccallcapapulseC585s = getpulseCapaCyclesC575s()
    pulsedC58CC20ccapa, legendspulseC58CC20, ccallcapaC58CC20 = getpulseCapaCyclesC58CC20()
    doublepulsedcapa, doublelegendspulse, ccallcapadouble = getdoublepulseCapaCycles()
    pulsedC5ccapa, legendspulseC5, ccallcapaC5  = getpulseCapaCyclesC5()

    allCapcities = [ccallcapa, ccallcapadouble, ccallcapapulse, cc823allcapa,   ccallcapaC5 , ccallcapaC58CC20,ccallcapapulseC585s]
    tempallCapcities = copy.deepcopy(allCapcities)

    fig, ax = plt.subplot_mosaic(
    [["top left", "top centre", "top right"],
        ["middle left", "middle centre", "middle right"],
     ["bottom row", "bottom row", "bottom row"]], sharey=True
    )
    plt.style.use(['science','nature','no-latex'])
    fig.tight_layout(h_pad=2.0)
    
    x = ["top left", "top centre", "top right",
        "middle left", "middle centre", "middle right",
        "bottom row"]

    labels = ['N(CC C/20) = ' + str(len(ccallcapa)),'N(Pulsed C/5.7) = ' + str(len(ccallcapapulse)), 'N(Pulsed C/5.7double) = ' + str(len(ccallcapadouble)), 'N(CC C/8.2) = ' + str(len(cc823allcapa)),
                        'N(Pulsed C/5) = ' + str(len(ccallcapaC5)), 'N(Mixed C/5.7 CC20) = ' + str(len(ccallcapaC58CC20)),
                        'N(Pulsed C/5.7 5s) = ' + str(len(ccallcapapulseC585s))
                        ]
    labels =  ['Reference CC C/20', 'Pulsed C/5.7' ,'Pulsed double length', 'CC C/8.2', 'Pulsed C/5',  'Mixed C/5.7+C/20', 'Pulsed 5s' ]
    colors = ['r','cornflowerblue' , 'cadetblue' ,'steelblue','forestgreen', 'darkgreen', 'darkolivegreen']
    for o in range(len(x)):
        for i in range(len(allCapcities[o])):
            cycles = np.linspace(1,len(allCapcities[o][i])-1,len(allCapcities[o][i])-1 )
            cccvref=ax[x[o]].scatter(cycles, allCapcities[o][i][1:], s = 0.05, label = labels[o], color = colors[o])
        ax[x[o]].set_xlim([0,300])
        ax[x[o]].set_xticks([0,100,200,300])
        ax[x[o]].set_xticklabels([0,100,200,300],fontsize=4)
        ax[x[o]].set_title(labels[o], fontdict={'color': colors[o]}, fontsize = 5)
    
    #ax["middle left"].set_ylabel('Discharge Capacity [mAh]', fontsize=6)
    ax["middle left"].set_ylabel('discharge capacity [mAh]', fontsize=6)
    ax["bottom row"].set_xlabel('cycle number', fontsize=6)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('DischargecapacityVSCyclenumberABSOLUTEALL5cells.png', dpi=800)
    os.chdir(origpath)
    plt.show()


def WLCTshowPlot():
    pathdata = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\WLTC Drive Test"
    orig = os.getcwd()
    os.chdir(pathdata)
    wlctfile = "WLCTPowerClass3.txt"
    with open(wlctfile, 'r') as file:
        currentdata = file.readlines()

    wlctfile = "WLTC_Drive_CycleClass3Speed.txt"
    with open(wlctfile, 'r') as file:
        speeddata = file.readlines() 

    time = []
    currentvalue = []
    for line in currentdata:
        parts = line.split()
        if len(parts) == 2:
            time.append(float(parts[0]))
            currentvalue.append(float(parts[1]))

    timespeed = []
    speedvalue = []
    for line in speeddata:
        parts = line.split()
        if len(parts) == 2:
            timespeed.append(float(parts[0]))
            speedvalue.append(float(parts[1]))

    fig,ax = plt.subplots(1, 1, figsize=(4, 2), sharex=False)
    plt.style.use(['science','nature','no-latex'])
    ax.axhspan(0.00001, 0.0065, xmin=0.1, xmax=0.99 , color='#008000', alpha=0.2) # red
    ax.axhspan(-0.00001, -0.0085, xmin=0.1, xmax=0.99, color='#d62728', alpha=0.2) # green 
    ax.axhline(xmin =0.1,xmax= 0.99, linewidth=2, color='#808080') # grey
    ax.axvline(x = 589, linewidth=0.75, color='#000000') # black
    ax.axvline(x = 1022, linewidth=0.75, color='#000000') # black
    ax.axvline(x = 1477, linewidth=0.75, color='#000000') # black
    file = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\WLTC Drive Test\AutoPP.png"
    logo = image.imread(file)
    img_extent = (-150,-5,-0.001, 0.0013)
    ax.imshow(logo, extent=img_extent, aspect='auto', zorder=-1)
    #ax.title.set_text("WLCT drive cycle")
    #ax.title.set_size(12)
    ax.set_xlabel("Time [s]", fontsize=8)   
    ax.set_ylabel("Current [mAh]", fontsize=8)   
    ax.set_xlim([-200,max(time)+25])
    ax.set_ylim([-0.0085,0.0065])
    ax.minorticks_off()
    ax.plot(time,currentvalue, linewidth=0.8, color = 'k', alpha= 0.5)

    ax.text(0.13, 0.95, "Charge\nRecuperation", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'green')
    
    ax.text(0.13, 0.1, "Discharge\nDrive&Accelerate", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'red')
    
    ax.text(0.22, 1.05, "low", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'k')
    ax.text(0.445, 1.05, "medium", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'k')
    ax.text(0.67, 1.05, "high", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'k')
    ax.text(0.838, 1.05, "extra high", transform=ax.transAxes, fontsize=14,
        verticalalignment='top',  alpha=0.75, color = 'k')
    
    ax2 = ax.twinx()
    ax2.minorticks_off() 
    ax2.set_ylim([-212.5,162.5])
    # Set ticks and tick labels for the second y-axis
    new_ticks = [-200,-150,-100,-50, 0,50,100, 150 ]
    ax2.set_yticks(new_ticks)
    ax2.yaxis.label.set_color('red') 
    ax2.tick_params(axis='y', colors='red')
    yticks = ax2.yaxis.get_major_ticks()

    yticks[0].set_visible(False) # -200
    yticks[1].set_visible(False) # -150
    yticks[2].set_visible(False) # -50
    yticks[3].set_visible(False) # -100
    ax2.plot(timespeed,speedvalue, linewidth=1, color = 'r')
    ax2.set_ylabel("Speed [km/h]", fontsize=8, color= 'r')  

    os.chdir(r'C:\Users\Leon Fischer\Desktop\PhD\My papers\Paper 1')
    plt.savefig('WLCTCycle.png', dpi=800)
    os.chdir(orig)
    plt.show()
    os.chdir(orig)


def WLCTcyclePlot():
    pathdata = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\WLTC Drive Test"
    orig = os.getcwd()
    os.chdir(pathdata)
    namecriteria = '_CCCV.hdf5'

    voltlists = []
    timelists = []
    currentlists = []
    dcCapacheckup = []
    
    content = os.listdir()
    content.sort(key=natural_keys)
    for i in range(len(content)):
        tempcapa = []
        if namecriteria in content[i]:
            print(content[i])
            data = loadHDF5(content[i])
            rawdata = data['raw']    
            testTime=rawdata['Testtime(s)']              
            stepIndex = rawdata['Step_Index']
            I = rawdata['Current(A)']
            V = rawdata['Voltage(V)']
            cyclelimitlist = []
            currentlist = []
            for l in range(len(I)-1):
                currentlist.append(I[l])
                if (stepIndex[l]== 12 and stepIndex[l-1]!= 12 or stepIndex[l]==13 and stepIndex[l-1]==12):
                    cyclelimitlist.append(l)
            for i in range(1, len(cyclelimitlist), 2):
                cycletesttime = []
                if i == 0:
                    minlimit = 0
                else:
                    minlimit = cyclelimitlist[i-1]
                correction = testTime[minlimit]
                for l in range(minlimit,cyclelimitlist[i]):
                    cycletesttime.append(testTime[l]-correction)
                currentlist = []
                voltlist = []
                cycleend = cyclelimitlist[i]
                cyclestart = cyclelimitlist[i-1]
                if i == 0:
                    cycleend = cyclelimitlist[0]
                    cyclestart = 0
                for l in range(cyclestart,cycleend):
                    voltlist.append(V[l])
                    currentlist.append(I[l])
                voltlists.append(voltlist)
                currentlists.append(currentlist)
                timelists.append(cycletesttime)
                try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del cycletesttime[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del cycletesttime[-1]
                except:
                        pass
                c = abs(np.trapz(currentlist, cycletesttime) * (1000/3600))
                tempcapa.append(c)
            dcCapacheckup.append(tempcapa)


    fig,ax = plt.subplots(1, 1, figsize=(4, 2), sharex=True)
    plt.style.use(['science','nature','no-latex'])    
    for u in range(len(dcCapacheckup)):
        if len(dcCapacheckup[u]) == 6:
            checkupcycles = [50, 100, 150, 200, 250, 300]
        else:
            checkupcycles = [50, 100, 150, 200, 250]
        ax.plot(checkupcycles, dcCapacheckup[u], color = 'r')
        ax.set_xlabel('Check-Cycle', fontsize=8)
        ax.set_ylabel('Discharge capacity [mAh]', fontsize=8)
        #ax.set_xlim([3.6, 4.2])
    plt.show()




    orig = os.getcwd()
    os.chdir(pathdata)
    namecriteriapulsed = 'Pulsed.hdf5'

    voltlistspulse = []
    timelistspulse = []
    currentlistspulse = []
    dcCapacheckuppulse = []
    
    content = os.listdir()
    content.sort(key=natural_keys)
    for i in range(len(content)):
        tempcapa = []
        if namecriteriapulsed in content[i]:
            print(content[i])
            data = loadHDF5(content[i])
            rawdata = data['raw']    
            testTime = rawdata['Testtime(s)']              
            stepIndex = rawdata['Step_Index']
            I = rawdata['Current(A)']
            V = rawdata['Voltage(V)']
            cyclelimitlistpulse = []
            for l in range(len(I)-1):
                if (stepIndex[l]== 13 and stepIndex[l-1]!= 13 or stepIndex[l]==14 and stepIndex[l-1]==13):
                    cyclelimitlistpulse.append(l)
            for i in range(1, len(cyclelimitlistpulse), 2):
                cycletesttimepulse = []
                if i == 0:
                    minlimit = 0
                else:
                    minlimit = cyclelimitlistpulse[i-1]
                correction = testTime[minlimit]
                for l in range(minlimit,cyclelimitlistpulse[i]):
                    cycletesttimepulse.append(testTime[l]-correction)
                currentlist = []
                voltlist = []
                cycleend = cyclelimitlistpulse[i]
                cyclestart = cyclelimitlistpulse[i-1]
                if i == 0:
                    cycleend = cyclelimitlistpulse[0]
                    cyclestart = 0
                for l in range(cyclestart,cycleend):
                    voltlist.append(V[l])
                    currentlist.append(I[l])
                voltlistspulse.append(voltlist)
                currentlistspulse.append(currentlist)
                timelistspulse.append(cycletesttimepulse)
                try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del cycletesttimepulse[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del cycletesttimepulse[-1]
                except:
                        pass
                c = abs(np.trapz(currentlist, cycletesttimepulse) * (1000/3600))
                tempcapa.append(c)
            dcCapacheckuppulse.append(tempcapa)
    
    fig,ax = plt.subplots(1, 1, figsize=(4, 2), sharex=True)
    plt.style.use(['science','nature','no-latex'])   
    for u in range(len(dcCapacheckup)):
        if len(dcCapacheckup[u]) == 6:
            checkupcycles = [50, 100, 150, 200, 250, 300]
        else:
            checkupcycles = [50, 100, 150, 200, 250]
        ax.plot(checkupcycles, dcCapacheckup[u], color = 'r')
        ax.set_xlabel('Check-Cycle', fontsize=8)
        ax.set_ylabel('Discharge capacity [mAh]', fontsize=8)
    
    for u in range(len(dcCapacheckuppulse)):
        if len(dcCapacheckuppulse[u]) == 6:
            checkupcycles = [50, 100, 150, 200, 250, 300]
        else:
            checkupcycles = [50, 100, 150, 200]
        ax.plot(checkupcycles, dcCapacheckuppulse[u], color = 'b')
        ax.set_xlabel('Check-Cycle', fontsize=8)
        ax.set_ylabel('Discharge capacity [mAh]', fontsize=8)
    legend = ['CCCV','_','_','_','Pulsed']
    plt.legend(legend)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('WLCTCheckupCyclesCapa.png', dpi=800)
    os.chdir(origpath)
    plt.show()


######################################           ZSW          ########################################################

pathpulse = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\PulsedC58"
pathC57ref1 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 2 C58 CC20rest batch 1'
pathC57ref2= r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 1 C58 5s pulse batch 2'
pathCCref = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\CCCV_C20"
pathCCrefC82 = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\C823 CCCV"
#pathCCref = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C20 CC - eff C20 pulsed ref'

pathdataZSW = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Pulsed\Pulsed_Formation"
pathdataRefZSW = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Reference\Formierung'
pathdataZSWEOL = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Pulsed\End of Line Test"


plotfolder = r'C:\Users\Leon\Desktop\PhD\KIT\My papers\Formation Strategies'


def VoltagevstransferCapacityFormation():

    areacoin = 0.7*0.7*np.pi
    areaphev1 = 11737  # aus 2,13 mAh/cm^2 = 21.3 Ah/m^2 und Capa= 25 Ah da doppelt coated Fläche x2

    def getpulseform():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapacharge = []
        pulseVoltagecharge=[]
        pulsecapadischarge = []
        pulseVoltagedischarge=[]
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000)/3600)/areacoin)
                pulsecapacharge.append(tempcapacharge)
                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelist,voltlist)
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]>= -0.000164:
                #        print("Len discharge CV: "+str(timelistcv[o]/3600))
                #        break
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))   
                tempcapadischarge.append(tempcapadischarge[-1])
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000 )/(3600*areacoin)))
                pulsecapadischarge.append(tempcapadischarge)
            #c = np.trapz(currentlist, timelist) /(3600)
            # charge 0.0036540186458333306 Ah
            # dc -0.0018406808777777744
            # cv -0.0012846443027777776
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))
                pulsecapacharge.append(tempcapacharge)

                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]>= -0.000164:
                #        print("Len discharge CV: "+str(timelistcv[o]/3600))
                #        break
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))      
                tempcapadischarge.append(tempcapadischarge[-1])  
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000 )/(3600*areacoin)))
                pulsecapadischarge.append(tempcapadischarge)
            #c = np.trapz(currentlist, timelist) /(3600)

        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))
                pulsecapacharge.append(tempcapacharge)

                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]>= -0.000164:
                #        print("Len discharge CV: "+str(timelistcv[o]/3600))
                #        break
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))     
                tempcapadischarge.append(tempcapadischarge[-1])   
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000)/(3600*areacoin)))
                pulsecapadischarge.append(tempcapadischarge)
            #c = np.trapz(currentlist, timelist) /(3600)

        return pulsecapacharge, pulseVoltagecharge, pulsecapadischarge, pulseVoltagedischarge

    def getpulseformZSW():
        hdf5path = pathdataZSW
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        pulsecapacharge = []
        pulseVoltagecharge=[]
        pulsecapadischarge = []
        pulseVoltagedischarge=[]
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_2"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areaphev1)))
                pulsecapacharge.append(tempcapacharge)

                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelistcv, data['split'][cvdc]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]==1.249 or currentlistcv[o]==1.250 or currentlistcv[o]==1.246 :
                #        counter = o
                #    if currentlistcv[o]== -1.245:
                #        counter = o
                #print("Len discharge CV: "+str(timelistcv[counter]/3600))
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areaphev1)))   
                tempcapadischarge.append(tempcapadischarge[-1])
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000 )/(3600*areaphev1)))
                pulsecapadischarge.append(tempcapadischarge)

        return pulsecapacharge, pulseVoltagecharge, pulsecapadischarge, pulseVoltagedischarge

    pulsecapacharge, pulseVoltagecharge, pulsecapadischarge, pulseVoltagedischarge = getpulseform()
    pulsecapachargeZSW, pulseVoltagechargeZSW, pulsecapadischargeZSW, pulseVoltagedischargeZSW = getpulseformZSW()

    for i in range(len(pulseVoltagecharge)):
        print(len(pulseVoltagecharge[i])/3600)
    
    for i in range(len(pulseVoltagedischargeZSW)):
        print(len(pulseVoltagedischargeZSW[i])/3600)

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1,1)
    plt.style.use(['science','nature','no-latex'])

    legendcolors = ["r", "b"]
    legendtext = ['N(CR2032) = '+str(len(pulsecapacharge)),  'N(PHEV1) = '+ str(len(pulsecapachargeZSW))]

    for i in range(len(pulsecapacharge)):
        ax.plot(pulsecapacharge[i],pulseVoltagecharge[i], color= legendcolors[0] )
        coin = ax.plot(pulsecapadischarge[i],pulseVoltagedischarge[i], color= legendcolors[0])

    for i in range(len(pulsecapachargeZSW)):
        ax.plot(pulsecapachargeZSW[i],pulseVoltagechargeZSW[i], color= legendcolors[1])
        phev1 = ax.plot(pulsecapadischargeZSW[i],pulseVoltagedischargeZSW[i], color= legendcolors[1])
    

    cr2032_handle = mlines.Line2D([], [], color=legendcolors[0], label=legendtext[0])
    phev1_handle = mlines.Line2D([], [], color=legendcolors[1], label=legendtext[1])

    legend = ax.legend(handles=[cr2032_handle, phev1_handle], bbox_to_anchor=(0.54, 0, 0.6, 0.22), loc="upper left",
                     borderaxespad=0, fancybox=True, shadow=True, frameon=True, ncol=1)
    
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    ax.set_ylim([2.85, 4.25])

    ax.set_xlabel('Area specific capacity [mAh / cm\u00B2]', fontsize=8)
    ax.set_xlabel('Area specific capacity [mAh / cm\u00B2]', fontsize=8)
    ax.set_ylabel('Voltage [V]', fontsize=8)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('FormPHEV1vsCoin5cells.png', dpi=1000)
    os.chdir(origpath)
    #plt.clf()
    plt.show()


def VoltagevstransferCapacityFormation2erPlot():

    areacoin = 0.7*0.7*np.pi
    areaphev1 = 11120  # aus 2,13 mAh/cm^2 = 21.3 Ah/m^2 und Capa= 25 Ah da doppelt coated Fläche x2

    def getpulseformCoin():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapacharge = []
        pulseVoltagecharge=[]
        pulsecapadischarge = []
        pulseVoltagedischarge=[]
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000)/3600)/areacoin)
                pulsecapacharge.append(tempcapacharge)
                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelist,voltlist)
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]>= -0.000164:
                #        print("Len discharge CV: "+str(timelistcv[o]/3600))
                #        break
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))   
                tempcapadischarge.append(tempcapadischarge[-1])
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000 )/(3600*areacoin)))
                pulsecapadischarge.append(tempcapadischarge)
            #c = np.trapz(currentlist, timelist) /(3600)
            # charge 0.0036540186458333306 Ah
            # dc -0.0018406808777777744
            # cv -0.0012846443027777776
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))
                pulsecapacharge.append(tempcapacharge)

                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]>= -0.000164:
                #        print("Len discharge CV: "+str(timelistcv[o]/3600))
                #        break
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))      
                tempcapadischarge.append(tempcapadischarge[-1])  
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000 )/(3600*areacoin)))
                pulsecapadischarge.append(tempcapadischarge)
            #c = np.trapz(currentlist, timelist) /(3600)

        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))
                pulsecapacharge.append(tempcapacharge)

                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]>= -0.000164:
                #        print("Len discharge CV: "+str(timelistcv[o]/3600))
                #        break
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))     
                tempcapadischarge.append(tempcapadischarge[-1])   
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000)/(3600*areacoin)))
                pulsecapadischarge.append(tempcapadischarge)
            #c = np.trapz(currentlist, timelist) /(3600)

        return pulsecapacharge, pulseVoltagecharge, pulsecapadischarge, pulseVoltagedischarge

    def getrefformCoin():
        hdf5path = pathCCref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'EOL_'

        pulsecapacharge = []
        pulseVoltagecharge=[]
        pulsecapadischarge = []
        pulseVoltagedischarge=[]
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000)/3600)/areacoin)
                pulsecapacharge.append(tempcapacharge)
                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelist,voltlist)
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]>= -0.000164:
                #        print("Len discharge CV: "+str(timelistcv[o]/3600))
                #        break
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areacoin)))   
                tempcapadischarge.append(tempcapadischarge[-1])
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000 )/(3600*areacoin)))
                pulsecapadischarge.append(tempcapadischarge)     

        return pulsecapacharge, pulseVoltagecharge, pulsecapadischarge, pulseVoltagedischarge

    def getpulseformZSW():
        hdf5path = pathdataZSW
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        pulsecapacharge = []
        pulseVoltagecharge=[]
        pulsecapadischarge = []
        pulseVoltagedischarge=[]
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_2"
                try:
                    currentlist = data['split'][l][str(l)+'_Current(A)']
                    voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    timelist = data['split'][l][str(l)+'_Testtime(s)']
                except:
                    currentlist = data['split'][l]['I']
                    timelist = data['split'][l]['t']
                    voltlist = data['split'][l]['V']
                    #plt.plot(timelist, data['split'][l]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                pulseVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areaphev1)))
                pulsecapacharge.append(tempcapacharge)

                dc = "Step_3"
                cvdc = "Step_4"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                try:
                    currentlist = data['split'][dc][str(dc)+'_Current(A)'] 
                    currentlistcv =data['split'][cvdc][str(cvdc)+'_Current(A)']
                    voltlist = data['split'][dc][str(dc)+'_Voltage(V)']
                    voltlistcv = data['split'][cvdc][str(cvdc)+'_Voltage(V)']
                    timelist = data['split'][dc][str(dc)+'_Testtime(s)']
                    timelistcv =  data['split'][cvdc][str(cvdc)+'_Testtime(s)']
                except:
                    currentlist = data['split'][dc]['I']
                    currentlistcv = data['split'][cvdc]['I']
                    timelist = data['split'][dc]['t'] 
                    timelistcv = data['split'][cvdc]['t']
                    voltlist = data['split'][dc]['V']
                    voltlistcv = data['split'][cvdc]['V']
                    #plt.plot(timelistcv, data['split'][cvdc]['V'])
                    #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                pulseVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]==1.249 or currentlistcv[o]==1.250 or currentlistcv[o]==1.246 :
                #        counter = o
                #    if currentlistcv[o]== -1.245:
                #        counter = o
                #print("Len discharge CV: "+str(timelistcv[counter]/3600))
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff*1000 )/(3600*areaphev1)))   
                tempcapadischarge.append(tempcapadischarge[-1])
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv*1000 )/(3600*areaphev1)))
                pulsecapadischarge.append(tempcapadischarge)

        return pulsecapacharge, pulseVoltagecharge, pulsecapadischarge, pulseVoltagedischarge

    def getRefZSW():
        hdf5path = pathdataRefZSW
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        refcapacharge = []
        refVoltagecharge=[]
        refcapadischarge = []
        refVoltagedischarge=[]
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcapacharge = [0]
                l = "Step_1"
                currentlist = data['split'][l]['I']
                timelist = data['split'][l]['t']
                voltlist = data['split'][l]['V']
                #plt.plot(timelist, data['split'][l]['V'])
                #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                except:
                    pass
                refVoltagecharge.append(voltlist)
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapacharge.append(tempcapacharge[-1] + ((float(currentlist[l]) * timediff *1000)/(3600*areaphev1)))
                refcapacharge.append(tempcapacharge)

                dc = "Step_2"
                cvdc = "Step_3"
                tempcapadischarge = []
                tempcapadischarge.append(tempcapacharge[-1])
                currentlist = data['split'][dc]['I']
                currentlistcv = data['split'][cvdc]['I']
                timelist = data['split'][dc]['t'] 
                timelistcv = data['split'][cvdc]['t']
                voltlist = data['split'][dc]['V']
                voltlistcv = data['split'][cvdc]['V']
                #plt.plot(timelistcv, data['split'][cvdc]['V'])
                #plt.show()
                try:
                    if currentlist[0] != currentlist[1]:
                        del currentlist[0]
                        del timelist[0]
                        del voltlist[0]    
                    if currentlist[-1] != currentlist[-2]:
                        del currentlist[-1]
                        del timelist[-1]
                        del voltlist[-1]
                    if currentlistcv[0] != currentlistcv[1]:
                        del currentlistcv[0]
                        del timelistcv[0]
                        del voltlistcv[0]    
                    if currentlistcv[-1] != currentlistcv[-2]:
                        del currentlistcv[-1]
                        del timelistcv[-1]
                        del voltlistcv[-1]
                except:
                    pass
                refVoltagedischarge.append(voltlist + voltlistcv)
                #print("Len discharge 1C: "+ str(timelist[-1]/3600))
                #for o in range(len(currentlistcv)):
                #    if currentlistcv[o]==1.249 or currentlistcv[o]==1.250 or currentlistcv[o]==1.246 :
                #        counter = o
                #    if currentlistcv[o]== -1.245:
                #        counter = o
                #print("Len discharge CV: "+str(timelistcv[counter]/3600))
                #print("cv limit not reached")
                for l in range(1,len(currentlist)):
                    timediff = timelist[l]- timelist[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlist[l] * timediff * 1000 )/(3600*areaphev1)))   
                tempcapadischarge.append(tempcapadischarge[-1])
                for l in range(1,len(currentlistcv)):
                    timediffcv = timelistcv[l]- timelistcv[l-1]
                    tempcapadischarge.append(tempcapadischarge[-1] + ((currentlistcv[l] * timediffcv * 1000 )/(3600*areaphev1)))
                refcapadischarge.append(tempcapadischarge)

        return refcapacharge, refVoltagecharge, refcapadischarge, refVoltagedischarge


    pulsecapacharge, pulseVoltagecharge, pulsecapadischarge, pulseVoltagedischarge = getpulseformCoin()

    refcapacharge, refVoltagecharge, refcapadischarge, refVoltagedischarge = getrefformCoin()

    pulsecapachargeZSW, pulseVoltagechargeZSW, pulsecapadischargeZSW, pulseVoltagedischargeZSW = getpulseformZSW()

    refcapachargeZSW, refVoltagechargeZSW, refcapadischargeZSW, refVoltagedischargeZSW = getRefZSW()

    for i in range(len(pulseVoltagecharge)):
        print(len(pulseVoltagecharge[i])/3600)
    
    for i in range(len(pulseVoltagedischargeZSW)):
        print(len(pulseVoltagedischargeZSW[i])/3600)
    
    for i in range(len(refVoltagedischargeZSW)):
        print(len(refVoltagedischargeZSW[i])/3600)


    fig,ax = plt.subplots(1,2 , figsize=(5, 3), sharex=True, sharey=True)
    fig.tight_layout(w_pad=2.0)
    plt.style.use(['science','nature','no-latex'])

    legendcolors = ["r", "violet", "cornflowerblue", "purple"]
    legendtext = ['CC C/20', 'CC C/20', 'Pulsed C/5.7',  'Pulsed C/5.7']

    ax[0].set_title('Coin cell')
    ax[1].set_title('PHEV1 cell')

    for i in range(len(refcapacharge)):
        print(max(refcapacharge[i]))
        ax[0].plot(refcapacharge[i],refVoltagecharge[i], color= legendcolors[0])
        coinref = ax[0].plot(refcapadischarge[i],refVoltagedischarge[i], color= legendcolors[0])
    
    for i in range(len(pulsecapacharge)):
        ax[0].plot(pulsecapacharge[i],pulseVoltagecharge[i], color= legendcolors[2] )
        coin = ax[0].plot(pulsecapadischarge[i],pulseVoltagedischarge[i], color= legendcolors[2])

    for i in range(len(refcapachargeZSW)):
        ax[1].plot(refcapachargeZSW[i],refVoltagechargeZSW[i], color= legendcolors[1])
        phev1ref = ax[1].plot(refcapadischargeZSW[i],refVoltagedischargeZSW[i], color= legendcolors[1])

    for i in range(len(pulsecapachargeZSW)):
        ax[1].plot(pulsecapachargeZSW[i],pulseVoltagechargeZSW[i], color= legendcolors[3])
        phev1 = ax[1].plot(pulsecapadischargeZSW[i],pulseVoltagedischargeZSW[i], color= legendcolors[3])
    
    cr2032ref_handle = mlines.Line2D([], [], color=legendcolors[0], label=legendtext[0])
    phev1ref_handle = mlines.Line2D([], [], color=legendcolors[1], label=legendtext[1])

    coin_handle = mlines.Line2D([], [], color=legendcolors[2], label=legendtext[2])
    phev1_handle = mlines.Line2D([], [], color=legendcolors[3], label=legendtext[3])

    legend = ax[0].legend(handles=[cr2032ref_handle, coin_handle], bbox_to_anchor=(0.57, 0, 0.4, 0.22), loc="upper left",
                     borderaxespad=0, fancybox=True, shadow=True, frameon=True, ncol=1, fontsize=6)
    
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    legend2 = ax[1].legend(handles=[ phev1ref_handle, phev1_handle], bbox_to_anchor=(0.57, 0, 0.6, 0.22), loc="upper left",
                     borderaxespad=0, fancybox=True, shadow=True, frameon=True, ncol=1, fontsize=6)
    
    legend2.get_frame().set_edgecolor('k')
    legend2.get_frame().set_linewidth(0.5)
    ax[0].text(-0.05, 1.1, 'a)', transform=ax[0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
    ax[1].text(-0.05, 1.1, 'b)', transform=ax[1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

    ax[0].set_ylim([2.85, 4.25])
    
    #fig.supxlabel(r'Area specific capacity [mAh / cm$^2$]', fontsize=8)
    fig.text(0.5, -0.02, r'Area specific capacity [mAh / cm$^2$]', ha='center', fontsize=8)


    ax[0].set_ylabel('Voltage [V]', fontsize=8)

    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('FormPHEV1vsCoin5cellstest.png', dpi=1000)
    os.chdir(origpath)
    #plt.clf()
    plt.show()


######################################           ZSW EOL

def ZSWEOLAuswertung():

    def getpulseEOL():
        hdf5path = pathpulse
        hdf5path2 = pathC57ref1
        hdf5path3 = pathC57ref2
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '_PulsedC58Temp'

        pulsecapa = []
        
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]["V"]
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[0]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,voltlist)
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                
        os.chdir(hdf5path2)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]["V"]
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[0]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                
        os.chdir(hdf5path3)
        namecriteria = 'test'
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                cycles = []
                capacities = []
                for n in data['split'].keys():
                    if re.search('Step_' +'.+',n) != None:
                        cycles.append(n)
                cycles.sort(key=natural_keys)
                for l in cycles:
                    try:
                        currentlist = data['split'][l][str(l)+'_Current(A)']
                        timelist = data['split'][l][str(l)+'_Testtime(s)']
                        voltlist = data['split'][l][str(l)+'_Voltage(V)']
                    except:
                        currentlist = data['split'][l]['I']
                        timelist = data['split'][l]['t']
                        voltlist = data['split'][l]["V"]
                        #plt.plot(timelist, data['split'][l]['V'])
                        #plt.show()
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[0]
                    except:
                        pass
                    if l == "Step_4":
                        for l in range(len(currentlist)):
                            if currentlist[l]> 0:
                                del currentlist[l:]
                                del timelist[l:]
                                del voltlist[l:]
                                break
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(c)
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                pulsecapa.append(capacities)
                
        return pulsecapa
    allCapaPulse = getpulseEOL()

    areaphev1 = 111120

    def getZSWEOL():
        hdf5path = pathdataZSWEOL
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        allCapa = []
        irvalues = []
        voltagesd = []
        timevoltsd = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                capacities=[]
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles),1):
                    print(l)
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                    #plt.plot(timelist, voltlist)
                    #plt.show()               
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) / (3600)
                    capacities.append(c)
                    if l == 12:
                        #print(c)
                        #print(voltlist[0])
                        #print(voltlist[-1])
                        #print(currentlist[2])
                        irvalues.append(abs((voltlist[0]- voltlist[-1])/currentlist[2])/areaphev1)
                    if l == 17:
                        voltagesd.append(voltlist)
                        timevoltsd.append(timelist)
                allCapa.append(capacities)

        return allCapa, irvalues, voltagesd, timevoltsd

    allCapa, irvalues, voltagesd, timevoltsd = getZSWEOL()

    #print(len(allCapa[2])) 
    #print(irvalues)
    #o = 4
    #plt.plot(timevoltsd[o], voltagesd[o])
    #plt.show()

    reverscapa = []
    for i in range(len(allCapa)):
        reverscapa.append((allCapa[i][21]-allCapa[i][17])*1000)


    irrreverscapa = []
    for i in range(len(allCapa)):
        irrreverscapa.append((allCapa[i][6] - allCapa[i][21])*1000)
    
    sdcurrent = []
    for i in range(len(allCapa)):
        sdcurrent.append(((reverscapa[i] + irrreverscapa[i])/(5*24)))

    print(reverscapa)
    print(irrreverscapa)
    print(sdcurrent)

    
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.top'] = False
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.bottom'] = False
    fig,ax = plt.subplots(2, 2)
    plt.style.use(['science','nature','no-latex'])
    fig.tight_layout(w_pad=1.0)
    fontsizelabelaxis = 4

    medianwidth = 0.5
    markersize = 2
    marker = 'D'
    boxlinewidth = 0.5

    IRZSW =ax[0,0].boxplot(irvalues, positions = [0], patch_artist=True, boxprops=dict(facecolor='cadetblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    ax[0,0].set_ylabel('Area specific impedance [Ohm / cm\u00B2]', fontsize=fontsizelabelaxis)
    ax[0,0].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,0].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,0].xaxis.set_ticklabels([])
    for median in IRZSW['medians']:
        median.set( linewidth = 1)
    
    vioZSWrev = ax[0,1].violinplot(reverscapa, positions=[0],showmeans=False, showextrema=True, showmedians=True)
    ax[0,1].set_ylabel('reversible Loss [mAh]', fontsize=fontsizelabelaxis)
    ax[0,1].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,1].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[0,1].xaxis.set_ticklabels([])
    for settings in vioZSWrev['bodies']:
        settings.set_facecolor('cadetblue')
        settings.set_alpha(1)
    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = vioZSWrev[partname]
        vp.set_edgecolor("k")
        vp.set_linewidth(0.5)
    vp.set_edgecolor('lime')
    vp.set_linewidth(1)

    vioZSWirrev = ax[1,0].violinplot(irrreverscapa, positions=[0], showmeans=False, showextrema=True, showmedians=True)
    ax[1,0].set_ylabel('irreversible Loss [mAh]', fontsize=fontsizelabelaxis)
    ax[1,0].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,0].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,0].xaxis.set_ticklabels([])
    for settings in vioZSWirrev['bodies']:
        settings.set_facecolor('cadetblue')
        settings.set_alpha(1)
    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = vioZSWirrev[partname]
        vp.set_edgecolor("k")
        vp.set_linewidth(0.5)
    vp.set_edgecolor('lime')
    vp.set_linewidth(1)
    
    SDcurrentZSW = ax[1,1].boxplot(sdcurrent, positions = [0], patch_artist=True, boxprops=dict(facecolor='cadetblue', linewidth=boxlinewidth), medianprops =dict(linewidth = medianwidth),whiskerprops=dict(linewidth = boxlinewidth), capprops = dict(linewidth = boxlinewidth), flierprops={'marker': marker, 'markersize': markersize, "markeredgewidth":boxlinewidth})
    ax[1,1].set_ylabel('self-discharge current [mA]', fontsize=fontsizelabelaxis)
    ax[1,1].yaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,1].xaxis.set_tick_params(labelsize=fontsizelabelaxis)
    ax[1,1].xaxis.set_ticklabels([])
    for median in SDcurrentZSW['medians']:
        median.set( linewidth = 1)
    
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('ZSWEOLAuswertung.png', dpi=1000)
    os.chdir(origpath)
    plt.show()


def ZSWEOL_IR_Auswertung():
    pathZSWEOLref = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Reference\EOL_Referenzformierung'
    areaphev1 = 111120


    def getZSWReferenceEOL_IR_SD():

        # auch nur 5 Zellen als vergleich daher alle anderen raus
        # alle vernachlässigten haben auch komische zwischenzyklen

        hdf5path = pathZSWEOLref
        
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        allCapa = []
        irvalues = []
        voltagesd = []
        timevoltsd = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                capacities=[]
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles),1):
                    #print(l)
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                    #print(len(cycles))
                    #plt.plot(timelist, voltlist)
                    #plt.show()               
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) / (3600)
                    capacities.append(c)
                    if len(cycles)==27 and l == 11:
                        #print(len(cycles))
                        #plt.plot(timelist, voltlist)
                        #plt.show() 
                        irvalues.append(abs((voltlist[0]- min(voltlist))/currentlist[2])/areaphev1)
                    if len(cycles)==28 and l == 12:
                        #print(len(cycles))
                        #plt.plot(timelist, voltlist)
                        #plt.show() 
                        irvalues.append(abs((voltlist[0]- min(voltlist))/currentlist[2])/areaphev1)
                    if len(cycles)==27 and l == 16 and len(voltlist)<400000:
                        voltagesd.append(voltlist)
                    if len(cycles)==28 and l == 16:
                        voltagesd.append(voltlist)
                        timevoltsd.append(timelist)
                allCapa.append(capacities)

        return allCapa, irvalues, voltagesd, timevoltsd
    
    allCapa, irvalues, voltagesd, timevoltsd  = getZSWReferenceEOL_IR_SD()

    print(irvalues)
    print(np.mean(irvalues))

    def getZSWPulsedEOL_IR_SD():
        hdf5path = pathdataZSWEOL
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        allCapa = []
        irvalues = []
        voltagesd = []
        timevoltsd = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                capacities=[]
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles),1):
                    print(l)
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                    #plt.plot(timelist, voltlist)
                    #plt.show()               
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) / (3600)
                    capacities.append(c)
                    if l == 12:
                        #print(c)
                        #print(voltlist[0])
                        #print(voltlist[-1])
                        #print(currentlist[2])
                        irvalues.append(abs((voltlist[0]- voltlist[-1])/currentlist[2])/areaphev1)
                    if l == 17:
                        voltagesd.append(voltlist)
                        timevoltsd.append(timelist)
                allCapa.append(capacities)

        return allCapa, irvalues, voltagesd, timevoltsd

    allCapaPulse, irvaluesPulse, voltagesdPulse, timevoltsdPulse = getZSWPulsedEOL_IR_SD()

    print(irvaluesPulse)
    print(np.mean(irvaluesPulse))


def ZSWEOLAuswertung():
    areaphev1 = 111120
    pathZSWEOLref = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Reference\EOL_Referenzformierung'
    
    def getZSWReferenceEOL_IR_SD():

        # auch nur 5 Zellen als vergleich daher alle anderen raus
        # alle vernachlässigten haben auch komische zwischenzyklen

        hdf5path = pathZSWEOLref
        
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        allCapa = []
        irvalues = []
        voltagesd = []
        timevoltsd = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                capacities=[]
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles),1):
                    #print(l)
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                    #print(len(cycles))
                    #plt.plot(timelist, voltlist)
                    #plt.show()               
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) / (3600)
                    capacities.append(c)
                    if len(cycles)==27 and l == 11:
                        #print(len(cycles))
                        #plt.plot(timelist, voltlist)
                        #plt.show() 
                        irvalues.append(abs((voltlist[0]- min(voltlist))/currentlist[2])/areaphev1)
                    if len(cycles)==28 and l == 12:
                        #print(len(cycles))
                        #plt.plot(timelist, voltlist)
                        #plt.show() 
                        irvalues.append(abs((voltlist[0]- min(voltlist))/currentlist[2])/areaphev1)
                    if len(cycles)==27 and l == 16 and len(voltlist)<400000:
                        voltagesd.append(voltlist)
                    if len(cycles)==28 and l == 16:
                        voltagesd.append(voltlist)
                        timevoltsd.append(timelist)
                allCapa.append(capacities)

        return allCapa, irvalues, voltagesd, timevoltsd
    
    allCapa, irvalues, voltagesd, timevoltsd  = getZSWReferenceEOL_IR_SD()

    def getZSWEOLPulsed():
        hdf5path = pathdataZSWEOL
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'

        allCapa = []
        irvalues = []
        voltagesd = []
        timevoltsd = []

        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                capacities=[]
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles),1):
                    print(l)
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                    #plt.plot(timelist, voltlist)
                    #plt.show()               
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) / (3600)
                    capacities.append(c)
                    if l == 12:
                        #print(c)
                        #print(voltlist[0])
                        #print(voltlist[-1])
                        #print(currentlist[2])
                        irvalues.append(abs((voltlist[0]- voltlist[-1])/currentlist[2])/areaphev1)
                    if l == 17:
                        voltagesd.append(voltlist)
                        timevoltsd.append(timelist)
                allCapa.append(capacities)

        return allCapa, irvalues, voltagesd, timevoltsd


    allCapapulse, irvaluespulse, voltagesdpulse, timevoltsdpulse = getZSWEOLPulsed()

    # already exclued the first 24 h
    def change_time_sqrt(timelist, voltlist):
        fittimepulse = []
        fitvoltpulse = []
        for p in range(len(timelist)):
            temptime = []
            tempvolt = []
            for time in range(len(timelist[p])):
                if timelist[p][time]/(60*60) > 24 :
                    temptime.append(np.sqrt(timelist[p][time]/(60*60)))
                    tempvolt.append(voltlist[p][time])
                #timelist[p][time] = timelist[p][time]/(60*60)
            fittimepulse.append(temptime)
            fitvoltpulse.append(tempvolt)
        return fittimepulse, fitvoltpulse

    fittimepulseZSW, fitvoltpulseZSW = change_time_sqrt(timevoltsdpulse, voltagesdpulse)
    fittimeZSW, fitvoltZSW = change_time_sqrt(timevoltsd, voltagesd)
    

    def linear_sqrt_function(x, a, b):
        return a * x + b
    
    def linear_fit(timelist, voltlist):
        voltdecayref = []
        rsqrt = []
        for r in range(len(timelist)):
            params, covariance = curve_fit(linear_sqrt_function, timelist[r], voltlist[r])
            a, b = params
            residuals = []
            totals = []
            for u in range(len(voltlist)):
                residuals.append((voltlist[r][u]- linear_sqrt_function(timelist[r][u], a, b))**2)
                totals.append((voltlist[r][u]-np.mean(voltlist[r]))**2)
            ss_res = np.sum(residuals)
            ss_tot = np.sum(totals)
            r_sq = np.round(1 - (ss_res / ss_tot), 5)
            voltdecayref.append(a *1000)
            rsqrt.append(r_sq)
        return voltdecayref, rsqrt

    voltdecayZSWpulsed, rsqrtZSWpulsed = linear_fit(fittimepulseZSW, fitvoltpulseZSW)

    voltdecayZSW, rsqrtZSW= linear_fit(fittimeZSW, fitvoltZSW)

    print(voltdecayZSWpulsed)
    print(np.mean(voltdecayZSWpulsed))
    print(rsqrtZSWpulsed)
    
    print(voltdecayZSW)
    print(np.mean(voltdecayZSW))
    print(rsqrtZSW)
    


######################################           ZSW Longterm

plotfolder = r'C:\Users\Leon\Desktop\PhD\KIT\My papers\Formation Strategies'
pathpulse = r"C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\PulsedC58\LongtermC58"
pathC57ref1 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 1 C58 5s pulse batch 2\300 Cycles'
pathC57ref2= r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C58 Pulse ref double pulse batch 2 C58 CC20rest batch 1\300 Cycles'
pathref = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\CCCV_C20\CCCV_InForm_300Cycles'
#pathref = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\EOL C20 CC - eff C20 pulsed ref\300Cycles'
pathdataZSWLongterm = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Pulsed\Zyklisierung'
pathdataZSWLongtermreference = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\ZSW Reference\Langzeit'
pathC82 = r'C:\Users\Leon\Desktop\PhD\Data\INFORM\EOL 20C\C823 CCCV\Longterm'

def longtermZSW():

    def getZSWLongtermCycles():
        hdf5path = pathdataZSWLongterm
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'
        checkupcycles = []
        checkupcyclenumber = []
        regularcycles = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcheckupcycles = []
                tempregularcycles = []
                tempcheckupcyclenumber = []
                counter = 1
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles)):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                                    
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    try:
                        if currentlist[0]<0 and currentlist[1]<0 and currentlist[2]<0:
                            if -1.255< currentlist[0] <-1.245:
                                tempcheckupcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                tempcheckupcyclenumber.append(counter)
                                counter = len(tempcheckupcyclenumber) * 50
                            else:
                                # need to check
                                tempcapa = np.trapz(currentlist, timelist) / (3600)
                                if tempcapa < -20:
                                    #print(tempcapa)
                                    if len(tempregularcycles)==0 or len(tempregularcycles)==1:
                                        tempregularcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                    elif len(tempregularcycles)>0 and abs(tempregularcycles[-1]-(abs(np.trapz(currentlist, timelist) / 3600)))<0.3 :
                                        tempregularcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                else:
                                    pass
                    except:
                        pass
                checkupcycles.append(tempcheckupcycles)
                checkupcyclenumber.append(tempcheckupcyclenumber)
                regularcycles.append(tempregularcycles)
                #print(checkupcycles)
                #print(checkupcyclenumber)
                #print(len(regularcycles[1]))
        return checkupcycles, checkupcyclenumber, regularcycles

    def getLongtermZSWReference():
        os.chdir(pathdataZSWLongtermreference)
        data = pd.read_excel("Cycling_Reference_ohne_checkup.xlsx")
        caparef = []
        cells = ["Cell 3014", "Cell 3015", "Cell 3017", "Cell 3018", "Cell 3019"]
        for cell in cells:
            tempcaparef = []
            capa_cell = list(data[cell])
            refvalue = copy.deepcopy(capa_cell[0]/1000)
            for i in range(0,300):
                tempcaparef.append(np.round((capa_cell[i]/1000)/refvalue, 3))
            caparef.append(tempcaparef)
        return caparef

    def getpulseCapaCycles():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'
        chargetimes = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        
        hdf5path = pathC57ref1
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)

        hdf5path = pathC57ref2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        return  legendsCC, allCapa, chargetimes

    def getrefCapaCycles():
        hdf5path = pathref
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'ZSWInform_CCCVFormC20'
        chargetimes = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,voltlist)
                    #plt.show()
                    #print(c)
                allCapa.append(capacities)
        return  legendsCC, allCapa, chargetimes


    checkupcycles, checkupcyclenumber, regularcycles = getZSWLongtermCycles()

    legendspulse, ccallcapapulse, chargetimes = getpulseCapaCycles()

    zswrefcapa = getLongtermZSWReference()

    legendsref, ccallcapadcref, chargetimesref = getrefCapaCycles()

    print(regularcycles[0][200])

    pulsecapabackup = copy.deepcopy(ccallcapapulse)

    # del cell 8 cause crazy behavior
    #del ccallcapapulse[1]
    #del ccallcapapulse[7]

    # del first dc after storage
    for i in range(len(regularcycles)):
        del regularcycles[i][0]

    cycleref = 0
    roundnumber = 3
    for i in range(len(checkupcycles)):
            reference = copy.deepcopy(regularcycles[i][cycleref])
            for l in range(len(checkupcycles[i])):
                temp = np.round(checkupcycles[i][l] / reference,roundnumber)
                checkupcycles[i][l] = temp

    for i in range(len(regularcycles)):
            reference = copy.deepcopy(regularcycles[i][cycleref])
            for l in range(len(regularcycles[i])):
                temp = np.round(regularcycles[i][l] / reference,roundnumber)
                regularcycles[i][l] = temp
    
    for i in range(len(ccallcapadcref)):
            reference = copy.deepcopy(ccallcapadcref[i][cycleref])
            for l in range(len(ccallcapadcref[i])):
                temp = np.round(ccallcapadcref[i][l]/ reference,roundnumber)
                ccallcapadcref[i][l] = temp
    
    for i in range(len(ccallcapapulse)):
            reference = copy.deepcopy(ccallcapapulse[i][cycleref])
            for l in range(len(ccallcapapulse[i])):
                temp = np.round(ccallcapapulse[i][l]/ reference,roundnumber)
                ccallcapapulse[i][l] = temp

    for i in range(len(regularcycles)):
        print(len(regularcycles[i]))
        #print(regularcycles[i][101])

    for i in range(len(ccallcapapulse)):
        print(ccallcapapulse[i][-1])

    # delete outlier in test 57
    del ccallcapapulse[1][1]

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1,1)
    plt.style.use(['science','nature','no-latex'])

    legendcolors = ["r", "darkviolet", "cornflowerblue", "purple"]
    legendtext = ['Coin cell ref', 'PHEV1 ref', 'Coin Pulsed',  'PHEV1 Pulsed']

    for i in range(len(ccallcapadcref)):
        cycles = np.linspace(1,len(ccallcapadcref[i]),len(ccallcapadcref[i]) )
        ax.plot(cycles, ccallcapadcref[i], color = legendcolors[0])
    
    for i in range(len(zswrefcapa)):
        cycles = np.linspace(1,len(zswrefcapa[i]),len(zswrefcapa[i]) )
        ax.plot(cycles, zswrefcapa[i], color = legendcolors[1])
    
    for i in range(len(ccallcapapulse)):
        cycles = np.linspace(1,len(ccallcapapulse[i]),len(ccallcapapulse[i]) )
        ax.plot(cycles, ccallcapapulse[i], color= legendcolors[2])

    for i in range(len(regularcycles)):
        cycles = np.linspace(1,len(regularcycles[i]),len(regularcycles[i]) )
        ax.plot(cycles, regularcycles[i], color = legendcolors[3])
       
    #for i in range(len(checkupcycles)):
    #    scat = ax.scatter(checkupcyclenumber[i], checkupcycles[i],color = legendcolors[2], marker= "*" ,s = 10)
    
    ax.set_xlim([0,300])
    ax.set_ylim([0.7,1.01])
    ax.set_xticks([0,100,200,300])
    #ax.set_yticks([0.8,1,1.5])
    #ax.set_yticklabels([0.8,1,1.5],fontsize=4)
    ax.set_xticklabels([0,100,200,300])

    cr2023_ref_handle = mlines.Line2D([], [], color=legendcolors[0], label=legendtext[0])
    phev1_ref_handle = mlines.Line2D([], [], color=legendcolors[1], label=legendtext[1])
    cr2032_handle = mlines.Line2D([], [], color=legendcolors[2], label=legendtext[2])
    phev1_handle = mlines.Line2D([], [], color=legendcolors[3], label=legendtext[3])
    
    #checkup_handle = mlines.Line2D([], [], marker="*", color='w', label=legendtext[2],  markerfacecolor=legendcolors[2],markersize=10)

    legend = ax.legend(handles=[cr2023_ref_handle, phev1_ref_handle, cr2032_handle, phev1_handle], bbox_to_anchor=(0.04, 0, 0.6, 0.29), loc="upper left",
                     borderaxespad=0, fancybox=True, shadow=True, frameon=True, ncol=1)
    
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    ax.set_xlabel('cycle number', fontsize=8)
    ax.set_ylabel('relative discharge capacity', fontsize=8)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('FormPHEV1vsCoinLongterm5cellsnew.png', dpi=1000)
    os.chdir(origpath)
    plt.show()
    #plt.clf()


def longtermZSWErrorbar():

    def getZSWLongtermCycles():
        hdf5path = pathdataZSWLongterm
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'
        checkupcycles = []
        checkupcyclenumber = []
        regularcycles = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcheckupcycles = []
                tempregularcycles = []
                tempcheckupcyclenumber = []
                counter = 1
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles)):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                                    
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    try:
                        if currentlist[0]<0 and currentlist[1]<0 and currentlist[2]<0:
                            if -1.255< currentlist[0] <-1.245:
                                tempcheckupcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                tempcheckupcyclenumber.append(counter)
                                counter = len(tempcheckupcyclenumber) * 50
                            else:
                                # need to check
                                tempcapa = np.trapz(currentlist, timelist) / (3600)
                                if tempcapa < -20:
                                    #print(tempcapa)
                                    if len(tempregularcycles)==0 or len(tempregularcycles)==1:
                                        tempregularcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                    elif len(tempregularcycles)>0 and abs(tempregularcycles[-1]-(abs(np.trapz(currentlist, timelist) / 3600)))<0.3 :
                                        tempregularcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                else:
                                    pass
                    except:
                        pass
                checkupcycles.append(tempcheckupcycles)
                checkupcyclenumber.append(tempcheckupcyclenumber)
                regularcycles.append(tempregularcycles)
                #print(checkupcycles)
                #print(checkupcyclenumber)
                #print(len(regularcycles[1]))
        return checkupcycles, checkupcyclenumber, regularcycles

    def getLongtermZSWReference():
        os.chdir(pathdataZSWLongtermreference)
        data = pd.read_excel("Cycling_Reference_ohne_checkup.xlsx")
        caparef = []
        cells = ["Cell 3014", "Cell 3015", "Cell 3017", "Cell 3018", "Cell 3019"]
        for cell in cells:
            tempcaparef = []
            capa_cell = list(data[cell])
            refvalue = copy.deepcopy(capa_cell[0]/1000)
            for i in range(0,300):
                tempcaparef.append(np.round((capa_cell[i]/1000)/refvalue, 3)*100)
            caparef.append(tempcaparef)
        return caparef

    def getpulseCapaCycles():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'
        chargetimes = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        
        hdf5path = pathC57ref1
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)

        hdf5path = pathC57ref2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        return  legendsCC, allCapa, chargetimes

    def getrefCapaCycles():
        hdf5path = pathC82
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'
        chargetimes = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,voltlist)
                    #plt.show()
                    #print(c)
                allCapa.append(capacities)
        return  legendsCC, allCapa, chargetimes


    checkupcycles, checkupcyclenumber, regularcycles = getZSWLongtermCycles()

    legendspulse, ccallcapapulse, chargetimes = getpulseCapaCycles()

    zswrefcapa = getLongtermZSWReference()

    legendsref, ccallcapadcref, chargetimesref = getrefCapaCycles()

    print(regularcycles[0][200])

    pulsecapabackup = copy.deepcopy(ccallcapapulse)

    # del cell 8 cause crazy behavior
    #del ccallcapapulse[1]
    #del ccallcapapulse[7]

    # del first dc after storage
    for i in range(len(regularcycles)):
        del regularcycles[i][0]

    cycleref = 0
    roundnumber = 3

    for i in range(len(regularcycles)):
            reference = copy.deepcopy(regularcycles[i][cycleref])
            for l in range(len(regularcycles[i])):
                temp = np.round(regularcycles[i][l] / reference,roundnumber)*100
                regularcycles[i][l] = temp
    
    for i in range(len(ccallcapadcref)):
            reference = copy.deepcopy(ccallcapadcref[i][cycleref])
            for l in range(len(ccallcapadcref[i])):
                temp = np.round(ccallcapadcref[i][l]/ reference,roundnumber)*100
                ccallcapadcref[i][l] = temp
    
    for i in range(len(ccallcapapulse)):
            reference = copy.deepcopy(ccallcapapulse[i][cycleref])
            for l in range(len(ccallcapapulse[i])):
                temp = np.round(ccallcapapulse[i][l]/ reference,roundnumber)*100
                ccallcapapulse[i][l] = temp

    # delete outlier in test 57
    del ccallcapapulse[1][1]

    def remove_outliers(data):
        filtered_data = []
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1 * IQR
        upper_bound = Q3 + 1 * IQR
        # Filter out values outside the IQR bounds
        filtered_row = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_data.append(filtered_row)
        return filtered_data

    cycleappear = 2

    pulsedZSWlong_avg =[]
    pulsedZSWlong_std = []
    for i in range(0,len(regularcycles[0]), cycleappear):
        temp =np.array([regularcycles[0][i],
                regularcycles[1][i],
                regularcycles[2][i],
                regularcycles[3][i],
                regularcycles[4][i]])
        tempfilterd = remove_outliers(temp)
        pulsedZSWlong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        pulsedZSWlong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))
    
    CCZSWreflong_avg =[]
    CCZSWreflong_std = []
    for i in range(0, len(zswrefcapa[0]), cycleappear):
        temp =np.array([zswrefcapa[0][i],
                zswrefcapa[1][i],
                zswrefcapa[2][i],
                zswrefcapa[3][i],
                zswrefcapa[4][i]])
        tempfilterd = remove_outliers(temp)
        CCZSWreflong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        CCZSWreflong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))
    
    pulsedCoinlong_avg =[]
    pulsedCoinlong_std = []
    for i in range(0, len(ccallcapapulse[0]), cycleappear):
        temp =np.array([ccallcapapulse[0][i],
                ccallcapapulse[1][i],
                ccallcapapulse[2][i],
                ccallcapapulse[3][i],
                ccallcapapulse[4][i]])
        tempfilterd = remove_outliers(temp)
        pulsedCoinlong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        pulsedCoinlong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))

    CCCoinlong_avg =[]
    CCCoinlong_std = []
    for i in range(0,len(ccallcapadcref[0]), cycleappear):
        temp = np.array([ccallcapadcref[0][i],
                ccallcapadcref[1][i],
                ccallcapadcref[2][i],
                ccallcapadcref[3][i],
                ccallcapadcref[4][i]])
        tempfilterd = remove_outliers(temp)
        CCCoinlong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        CCCoinlong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))

    for i in range(len(regularcycles)):
        print(len(regularcycles[i]))
        #print(regularcycles[i][101])

    for i in range(len(ccallcapapulse)):
        print(ccallcapapulse[i][-1])


    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1,1)
    plt.style.use(['science','nature','no-latex'])

    legendcolors = ["steelblue", "purple", "cornflowerblue", "violet"]
    legendtext = ['Coin cell CC C/8.2', 'PHEV1 CC C/20', 'Coin Pulsed C/5.7',  'PHEV1 Pulsed C/5.7']    

    alphaval = 0.5
    capsizeval = 1.5
    sval = 0.5

    cycles = np.linspace(1, int(len(CCCoinlong_avg)*2), len(CCCoinlong_avg), dtype= int)
    ax.errorbar(cycles, CCCoinlong_avg, yerr = CCCoinlong_std, color = legendcolors[0], alpha=alphaval, capsize = capsizeval )
    ax.scatter(cycles, CCCoinlong_avg, s= sval, color = legendcolors[0])

    cycles = np.linspace(1, int(len(pulsedZSWlong_avg)*2), len(pulsedZSWlong_avg), dtype= int)
    ax.errorbar(cycles, pulsedZSWlong_avg,yerr = pulsedZSWlong_std ,color = legendcolors[3], alpha=alphaval, capsize = capsizeval)
    ax.scatter(cycles, pulsedZSWlong_avg, s= sval, color = legendcolors[3])

    cycles = np.linspace(1, int(len(CCZSWreflong_avg)*2), len(CCZSWreflong_avg), dtype= int)
    ax.errorbar(cycles, CCZSWreflong_avg, yerr=CCZSWreflong_std ,color = legendcolors[1], alpha=alphaval, capsize = capsizeval)
    ax.scatter(cycles, CCZSWreflong_avg, s= sval, color = legendcolors[1])

    cycles = np.linspace(1, int(len(pulsedCoinlong_avg)*2), len(pulsedCoinlong_avg), dtype= int)
    ax.errorbar(cycles, pulsedCoinlong_avg, yerr= pulsedCoinlong_std ,color = legendcolors[2], alpha=alphaval, capsize = capsizeval)
    ax.scatter(cycles, pulsedCoinlong_avg, s= sval, color = legendcolors[2])


    sval_zoom = 0.7
    alphaval_zoom = 0.7
    ax_inset2 = inset_axes(ax, width="30%", height="30%", bbox_to_anchor=(0.005, 0.16, 1.0, 0.6), bbox_transform=ax.transAxes)

    cycle_zoom = 75
    ax_inset2.errorbar(cycles[cycle_zoom:101], pulsedZSWlong_avg[cycle_zoom:101],yerr = pulsedZSWlong_std[cycle_zoom:101] ,color = legendcolors[3], alpha=alphaval_zoom, capsize = capsizeval)
    ax_inset2.scatter(cycles[cycle_zoom:101], pulsedZSWlong_avg[cycle_zoom:101], s= sval_zoom, color = legendcolors[3])

    ax_inset2.errorbar(cycles[cycle_zoom:101], CCZSWreflong_avg[cycle_zoom:101], yerr=CCZSWreflong_std[cycle_zoom:101] ,color = legendcolors[1], alpha=alphaval_zoom, capsize = capsizeval)
    ax_inset2.scatter(cycles[cycle_zoom:101], CCZSWreflong_avg[cycle_zoom:101], s= sval_zoom, color = legendcolors[1])
    ax.indicate_inset_zoom(ax_inset2, edgecolor="black")
    #for i in range(len(checkupcycles)):
    #    scat = ax.scatter(checkupcyclenumber[i], checkupcycles[i],color = legendcolors[2], marker= "*" ,s = 10)
    
    ax.set_xlim([0,300])
    #ax.set_ylim([0.75,1.01])
    ax.set_xticks([0,100,200,300])
    #ax.set_yticks([0.8, 0.85, 0.9, 0.95, 1.0])
    #ax.set_yticklabels([0.8,1,1.5],fontsize=4)
    ax.set_xticklabels([0,100,200,300])

    cr2023_ref_handle = mlines.Line2D([], [], color=legendcolors[0], label=legendtext[0])
    phev1_ref_handle = mlines.Line2D([], [], color=legendcolors[1], label=legendtext[1])
    cr2032_handle = mlines.Line2D([], [], color=legendcolors[2], label=legendtext[2])
    phev1_handle = mlines.Line2D([], [], color=legendcolors[3], label=legendtext[3])
    
    #checkup_handle = mlines.Line2D([], [], marker="*", color='w', label=legendtext[2],  markerfacecolor=legendcolors[2],markersize=10)

    legend = ax.legend(handles=[cr2023_ref_handle,  cr2032_handle,phev1_ref_handle, phev1_handle], bbox_to_anchor=(0.03, -0.03, 0.6, 0.29), loc="upper left",
                     borderaxespad=0, fancybox=True, shadow=True, frameon=True, ncol=1, fontsize = 4)
    
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlabel('cycle number', fontsize=8)
    ax.set_ylabel('relative discharge capacity [%]', fontsize=8)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('FormPHEV1vsCoinLongterm5cellsErrorBars.png', dpi=1000)
    os.chdir(origpath)
    plt.show()
    #plt.clf()

def longtermZSWErrorbarPHEV1():

    def getZSWLongtermCycles():
        hdf5path = pathdataZSWLongterm
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'
        checkupcycles = []
        checkupcyclenumber = []
        regularcycles = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcheckupcycles = []
                tempregularcycles = []
                tempcheckupcyclenumber = []
                counter = 1
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles)):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                                    
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    try:
                        if currentlist[0]<0 and currentlist[1]<0 and currentlist[2]<0:
                            if -1.255< currentlist[0] <-1.245:
                                tempcheckupcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                tempcheckupcyclenumber.append(counter)
                                counter = len(tempcheckupcyclenumber) * 50
                            else:
                                # need to check
                                tempcapa = np.trapz(currentlist, timelist) / (3600)
                                if tempcapa < -20:
                                    #print(tempcapa)
                                    if len(tempregularcycles)==0 or len(tempregularcycles)==1:
                                        tempregularcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                    elif len(tempregularcycles)>0 and abs(tempregularcycles[-1]-(abs(np.trapz(currentlist, timelist) / 3600)))<0.3 :
                                        tempregularcycles.append(np.around(abs(np.trapz(currentlist, timelist) / 3600),2))
                                else:
                                    pass
                    except:
                        pass
                checkupcycles.append(tempcheckupcycles)
                checkupcyclenumber.append(tempcheckupcyclenumber)
                regularcycles.append(tempregularcycles)
                #print(checkupcycles)
                #print(checkupcyclenumber)
                #print(len(regularcycles[1]))
        return checkupcycles, checkupcyclenumber, regularcycles

    def getLongtermZSWReference():
        os.chdir(pathdataZSWLongtermreference)
        data = pd.read_excel("Cycling_Reference_ohne_checkup.xlsx")
        caparef = []
        cells = ["Cell 3014", "Cell 3015", "Cell 3017", "Cell 3018", "Cell 3019"]
        for cell in cells:
            tempcaparef = []
            capa_cell = list(data[cell])
            refvalue = copy.deepcopy(capa_cell[0]/1000)
            for i in range(0,300):
                tempcaparef.append(np.round((capa_cell[i]/1000)/refvalue, 3))
            caparef.append(tempcaparef)
        return caparef

    def getpulseCapaCycles():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'
        chargetimes = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        
        hdf5path = pathC57ref1
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)

        hdf5path = pathC57ref2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1])
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2])
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        return  legendsCC, allCapa, chargetimes

    def getrefCapaCycles():
        hdf5path = pathC82
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'
        chargetimes = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                print(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.around(np.trapz(currentlist, timelist),2 )* (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,voltlist)
                    #plt.show()
                    #print(c)
                allCapa.append(capacities)
        return  legendsCC, allCapa, chargetimes


    checkupcycles, checkupcyclenumber, regularcycles = getZSWLongtermCycles()

    legendspulse, ccallcapapulse, chargetimes = getpulseCapaCycles()

    zswrefcapa = getLongtermZSWReference()

    legendsref, ccallcapadcref, chargetimesref = getrefCapaCycles()

    print(regularcycles[0][200])

    pulsecapabackup = copy.deepcopy(ccallcapapulse)

    # del cell 8 cause crazy behavior
    #del ccallcapapulse[1]
    #del ccallcapapulse[7]

    # del first dc after storage
    for i in range(len(regularcycles)):
        del regularcycles[i][0]

    cycleref = 0
    roundnumber = 3

    for i in range(len(regularcycles)):
            reference = copy.deepcopy(regularcycles[i][cycleref])
            for l in range(len(regularcycles[i])):
                temp = np.round(regularcycles[i][l] / reference,roundnumber)
                regularcycles[i][l] = temp
    
    for i in range(len(ccallcapadcref)):
            reference = copy.deepcopy(ccallcapadcref[i][cycleref])
            for l in range(len(ccallcapadcref[i])):
                temp = np.round(ccallcapadcref[i][l]/ reference,roundnumber)
                ccallcapadcref[i][l] = temp
    
    for i in range(len(ccallcapapulse)):
            reference = copy.deepcopy(ccallcapapulse[i][cycleref])
            for l in range(len(ccallcapapulse[i])):
                temp = np.round(ccallcapapulse[i][l]/ reference,roundnumber)
                ccallcapapulse[i][l] = temp

    # delete outlier in test 57
    del ccallcapapulse[1][1]

    def remove_outliers(data):
        filtered_data = []
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1 * IQR
        upper_bound = Q3 + 1 * IQR
        # Filter out values outside the IQR bounds
        filtered_row = data[(data >= lower_bound) & (data <= upper_bound)]
        filtered_data.append(filtered_row)
        return filtered_data

    cycleappear = 2

    pulsedZSWlong_avg =[]
    pulsedZSWlong_std = []
    for i in range(0,len(regularcycles[0]), cycleappear):
        temp =np.array([regularcycles[0][i],
                regularcycles[1][i],
                regularcycles[2][i],
                regularcycles[3][i],
                regularcycles[4][i]])
        tempfilterd = remove_outliers(temp)
        pulsedZSWlong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        pulsedZSWlong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))
    
    CCZSWreflong_avg =[]
    CCZSWreflong_std = []
    for i in range(0, len(zswrefcapa[0]), cycleappear):
        temp =np.array([zswrefcapa[0][i],
                zswrefcapa[1][i],
                zswrefcapa[2][i],
                zswrefcapa[3][i],
                zswrefcapa[4][i]])
        tempfilterd = remove_outliers(temp)
        CCZSWreflong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        CCZSWreflong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))
    
    pulsedCoinlong_avg =[]
    pulsedCoinlong_std = []
    for i in range(0, len(ccallcapapulse[0]), cycleappear):
        temp =np.array([ccallcapapulse[0][i],
                ccallcapapulse[1][i],
                ccallcapapulse[2][i],
                ccallcapapulse[3][i],
                ccallcapapulse[4][i]])
        tempfilterd = remove_outliers(temp)
        pulsedCoinlong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        pulsedCoinlong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))

    CCCoinlong_avg =[]
    CCCoinlong_std = []
    for i in range(0,len(ccallcapadcref[0]), cycleappear):
        temp = np.array([ccallcapadcref[0][i],
                ccallcapadcref[1][i],
                ccallcapadcref[2][i],
                ccallcapadcref[3][i],
                ccallcapadcref[4][i]])
        tempfilterd = remove_outliers(temp)
        CCCoinlong_avg.append(np.around(np.mean(tempfilterd),roundnumber))
        CCCoinlong_std.append(np.around((np.std(tempfilterd)/2), roundnumber))

    for i in range(len(regularcycles)):
        print(len(regularcycles[i]))
        #print(regularcycles[i][101])

    for i in range(len(ccallcapapulse)):
        print(ccallcapapulse[i][-1])


    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1,1)
    plt.style.use(['science','nature','no-latex'])

    legendcolors = ["r", "purple", "cornflowerblue", "violet"]
    legendtext = ['PHEV1 reference CC C/20','PHEV1 Pulsed C/5.7']  
    alphaval = 0.7
    capsizeval = 1.5
    sval = 1

    cycles = np.linspace(1, int(len(pulsedZSWlong_avg)*2), len(pulsedZSWlong_avg), dtype= int)
    ax.errorbar(cycles, pulsedZSWlong_avg,yerr = pulsedZSWlong_std ,color = legendcolors[3], alpha=alphaval, capsize = capsizeval)
    ax.scatter(cycles, pulsedZSWlong_avg, s= sval, color = legendcolors[3])

    cycles = np.linspace(1, int(len(CCZSWreflong_avg)*2), len(CCZSWreflong_avg), dtype= int)
    ax.errorbar(cycles, CCZSWreflong_avg, yerr=CCZSWreflong_std ,color = legendcolors[1], alpha=alphaval, capsize = capsizeval)
    ax.scatter(cycles, CCZSWreflong_avg, s= sval, color = legendcolors[1])
    
    #for i in range(len(checkupcycles)):
    #    scat = ax.scatter(checkupcyclenumber[i], checkupcycles[i],color = legendcolors[2], marker= "*" ,s = 10)
    
    ax.set_xlim([0,300])
    #ax.set_ylim([0.75,1.01])
    ax.set_xticks([0,100,200,300])
    #ax.set_yticks([0.8, 0.85, 0.9, 0.95, 1.0])
    #ax.set_yticklabels([0.8,1,1.5],fontsize=4)
    ax.set_xticklabels([0,100,200,300])

    phev1_ref_handle = mlines.Line2D([], [], color=legendcolors[1], label=legendtext[0])
    phev1_handle = mlines.Line2D([], [], color=legendcolors[3], label=legendtext[1])
    
    #checkup_handle = mlines.Line2D([], [], marker="*", color='w', label=legendtext[2],  markerfacecolor=legendcolors[2],markersize=10)

    legend = ax.legend(handles=[phev1_ref_handle, phev1_handle], bbox_to_anchor=(0.04, 0.01, 0.6, 0.29), loc="upper left",
                     borderaxespad=0, fancybox=True, shadow=True, frameon=True, ncol=1, fontsize = 6)
    
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlabel('cycle number', fontsize=8)
    ax.set_ylabel('relative discharge capacity', fontsize=8)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('FormPHEV1vsCoinLongterm5cellsErrorBarsPHEV1.png', dpi=1000)
    os.chdir(origpath)
    plt.show()
    #plt.clf()


def Chargetimeanalysis():

    def getZSWLongtermCycles():
        hdf5path = pathdataZSWLongterm
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = '.hdf5'
        checkupcycles = []
        checkupcyclenumber = []
        regularcycles = []
        allCapa = []
        legendsCC = []
        chargetimes = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                print(content[i])
                data = loadHDF5(content[i])
                tempcheckupcycles = []
                tempregularcycles = []
                tempcheckupcyclenumber = []
                tempchargetimes = []
                counter = 1
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles),1):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    #if currentlist[0]>=1.2 and currentlist[0]<1.3:
                    #    print(l)
                                    
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    try:
                        if currentlist[0]<0 and currentlist[1]<0 and currentlist[2]<0:
                            if -1.255< currentlist[0] <-1.245:
                                tempcheckupcycles.append(abs(np.trapz(currentlist, timelist)) / 3600)
                                tempcheckupcyclenumber.append(counter)
                                counter = len(tempcheckupcyclenumber) * 50
                            else:
                                # need to check
                                tempcapa = np.trapz(currentlist, timelist) / (3600)
                                if tempcapa < -20:
                                    #print(tempcapa)
                                    if len(tempregularcycles)==0 or len(tempregularcycles)==1:
                                        tempregularcycles.append(abs(np.trapz(currentlist, timelist)) / 3600)
                                    elif len(tempregularcycles)>0 and abs(tempregularcycles[-1]-(abs(np.trapz(currentlist, timelist)) / 3600))<0.3 :
                                        tempregularcycles.append(abs(np.trapz(currentlist, timelist)) / 3600)
                                else:
                                    pass
                    except:
                        pass
                    
                    if currentlist[0]>10 and currentlist[1]>10 and currentlist[2]>10:
                        #print(l)
                        #plt.plot(timelist, voltlist)
                        #plt.show()    
                        temptemptimes = []
                        temptemptimes.append(timelist[-1]/3600)
                        counter = 0
                        toggle = True
                        for p in range(len(currentlist)):
                            if toggle == True and currentlist[p]< 0.00328:
                                toggle = False
                                temptemptimes.append(timelist[p-2]/3600)
                            if currentlist[p]<0:
                                counter = counter+1
                        temptemptimes.append(counter)
                        tempchargetimes.append(temptemptimes)
                chargetimes.append(tempchargetimes)
                checkupcycles.append(tempcheckupcycles)
                checkupcyclenumber.append(tempcheckupcyclenumber)
                regularcycles.append(tempregularcycles)
                #print(checkupcycles)
                #print(checkupcyclenumber)
                #print(len(regularcycles[1]))
        return checkupcycles, checkupcyclenumber, regularcycles, chargetimes

    def getpulseCapaCycles():
        hdf5path = pathpulse
        #hdf5path = "C:\\Users\\Leon Fischer\\Desktop\\PhD\\Data"
        namecriteria = 'test'
        chargetimes = []
        allCapa = []
        legendsCC = []
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                    #print(content[i])
                    #plt.plot(timelist,data['split'][l]['V'])
                    #plt.show()
                    #print(c)
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1]/3600)
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2]/3600)
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        
        hdf5path = pathC57ref1
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1]/3600)
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2]/3600)
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)

        hdf5path = pathC57ref2
        namecriteria = 'test_'
        os.chdir(hdf5path)
        originpath = os.getcwd()
        content = os.listdir()
        content.sort(key=natural_keys)
        for i in range(len(content)):
            if namecriteria in content[i]:
                data = loadHDF5(content[i])
                capacities = []
                tempchargetimes = []
                cycles = list(data['split'].keys())
                cycles.sort(key=natural_keys)
                for l in range(1,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    c = np.trapz(currentlist, timelist) * (1000/3600)
                    capacities.append(abs(c))
                for l in range(0,len(cycles), 2):
                    try:
                        currentlist = data['split'][str(cycles[l])][str(cycles[l])+'_Current(A)']
                        timelist = data['split'][str(cycles[l])][str(cycles[l])+'_Testtime(s)']
                        voltlist = data['split'][str(cycles[l])][str(cycles[l])+'Voltage(V)']
                    except:
                        currentlist = data['split'][str(cycles[l])]['I']
                        timelist = data['split'][str(cycles[l])]['t']
                        voltlist = data['split'][str(cycles[l])]['V']
                    try:
                        if currentlist[0] != currentlist[1]:
                            del currentlist[0]
                            del timelist[0]
                            del voltlist[0]
                        if currentlist[-1] != currentlist[-2]:
                            del currentlist[-1]
                            del timelist[-1]
                            del voltlist[-1]
                    except:
                        pass
                    #plt.plot(timelist,currentlist)
                    #plt.show()
                    temptemptimes = []
                    temptemptimes.append(timelist[-1]/3600)
                    counter = 0
                    toggle = True
                    for p in range(len(currentlist)):
                        if toggle == True and currentlist[p]< 0.00328:
                            toggle = False
                            temptemptimes.append(timelist[p-2]/3600)
                        if currentlist[p]<0:
                            counter = counter+1
                    temptemptimes.append(counter)
                    tempchargetimes.append(temptemptimes)
                allCapa.append(capacities)
                chargetimes.append(tempchargetimes)
        return  legendsCC, allCapa, chargetimes
    

    checkupcycles, checkupcyclenumber, regularcycles, chargetimesZSW = getZSWLongtermCycles()

    legendspulse, ccallcapapulse, chargetimes = getpulseCapaCycles()

    

    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
    fig,ax = plt.subplots(1,1)
    plt.style.use(['science','nature','no-latex'])

    legendcolors = ["r", "b", "cornflowerblue"]
    legendtext = ['N(CR2032) = '+str(len(ccallcapapulse)),  'N(PHEV1) = '+ str(len(regularcycles))]

    # check chargetimes zelle 10 weil über 300 cycles

    for i in range(len(chargetimes)):
        for o in range(len(chargetimes[i])):
            timeratio = chargetimes[i][o][1] / chargetimes[i][o][0]
            if timeratio >= 0.5:
                ax.scatter([o+1],chargetimes[i][o][0], color = 'r', s=0.2 )
            else:
                ax.scatter([o+1],chargetimes[i][o][0], color = 'darkred', s=0.2 )
    

    for i in range(len(chargetimesZSW)):
        for o in range(len(chargetimesZSW[i])):
            timeratio = chargetimesZSW[i][o][1] / chargetimesZSW[i][o][0]
            if timeratio >= 0.5:
                ax.scatter([o+1],chargetimesZSW[i][o][0], color = 'b', s=0.001 )
            else:
                ax.scatter([o+1],chargetimesZSW[i][o][0], color = 'cornflowerblue', s=0.001 )

    ax.set_xlim([0,300])
    #ax.set_ylim([0.52,1.12])
    ax.set_xticks([0,100,200,300])
    #ax.set_yticks([0.8,1,1.5])
    #ax.set_yticklabels([0.8,1,1.5],fontsize=4)
    ax.set_xticklabels([0,100,200,300])

    from matplotlib.patches import Patch
    pa1 = Patch(facecolor='r') 
    pb1 = Patch(facecolor='b')

    pa2 = Patch(facecolor='darkred') 
    pb2 = Patch(facecolor='cornflowerblue') 

    legend = ax.legend(handles=[pa1, pb1, pa2, pb2],
                        labels=['', '','N(CR2032) = '+str(len(ccallcapapulse)),  'N(PHEV1) = '+ str(len(regularcycles))] ,
                        handletextpad=0.5, handlelength=1.0, columnspacing=-0.5,
                        bbox_to_anchor=(0.5, 0, 0.6, 0.8), loc="upper left",borderaxespad=0, fancybox=True, shadow=True, frameon=True, ncol=2)
    
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_linewidth(0.5)

    ax.set_xlabel('cycle number [1]', fontsize=8)
    ax.set_ylabel('total time [h]', fontsize=8)
    origpath = os.getcwd()
    os.chdir(plotfolder)
    plt.savefig('PHEV1vsCoinLongtermChargetime5cells.png', dpi=1000)
    os.chdir(origpath)
    plt.clf()















"""
#ZSW    
    1.3588888888888888
Len discharge 1C: 0.890772222222224
Len discharge CV: 0.46790833333333165

    1.3811111111111112
Len discharge 1C: 0.8893555555555557
Len discharge CV: 0.4683111111111106

    1.395
Len discharge 1C: 0.8932888888888899
Len discharge CV: 0.46830833333333227

    1.3708333333333333
Len discharge 1C: 0.8929944444444441
Len discharge CV: 0.468313888888889

    1.3677777777777778


len1D = [
    0.5611891666666654,0.49451249999999974,0.7278800000000006, 0.7253247222222207,0.7378233333333345,0.7684066666666675,0.8295388888888899,0.7892358333333322,0.16590444444444427,0.46268861111111115,0.6475030555555552,
    0.6280566666666669,0.4279580555555559, 0.7058986111111113, 0.61417611111111, 0.5794336111111119, 0.42519361111111115, 0.8809591666666671,  0.5683158333333328, 0.422399722222222, 0.6863972222222219,
    0.8350380555555542, 0.7558400000000014, 0.7697233333333335, 0.8627127777777767, 0.8293761111111113, 0.9182847222222215, 0.8126991666666679, 0.5334613888888887
]
    
lenCV= [ 2.133590833333334, 1.3335402777777785 , 0.9724002777777767,1.2283716666666684,0.9573830555555549, 0.996316944444445,0.8962927777777784,1.1241758333333343,
        0.7142072222222229,  1.3797319444444445,  0.9490349999999994, 1.9495147222222213, 1.3186344444444442,1.7397055555555543,1.2589433333333322,1.6702180555555546,
        2.0175341666666666,  0.9976924999999998,  1.682747500000001,  1.8925202777777779,  1.6910186111111114,  0.903126666666665,  1.1643988888888896,  1.4463819444444461,
        0.8280469444444457, 0.91416111111111,  0.7599563888888891,  0.8919430555555563,  1.5711933333333319
        ]
"""