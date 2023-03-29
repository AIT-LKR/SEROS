import numpy as np
from os import listdir
from os.path import isfile, join
import os
import matplotlib.pyplot as plt
from scipy import signal

# dataToPlot = ['time', 'smoothedDensityDiff']

class dataSet:
    data = []
    def __init__(self, parameter):
        self.parameter = parameter

def str2IntOrFloat(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

currentDirectory = os.path.dirname(os.path.realpath(__file__))

print('')
print('Data directory:')
print('---------------')
print(currentDirectory)
print('')

parameterFile = os.path.dirname(os.path.realpath(__file__)) + '/output.dat'

allData = []
content = []
with open(parameterFile) as file:
    parameterList = file.readline().split(' ')

# print(parameterList)

time, density, energy = np.genfromtxt(parameterFile, delimiter=' ', skip_header=1, usecols=(0,8,10), unpack = True)

#time = time[1128:]
#density = density[1128:]
#energy = energy[1128:]

b, a = signal.butter(12,0.04)

densityFiltered = signal.filtfilt(b, a, density)

energyFiltered = signal.filtfilt(b, a, energy)

print('')
print('Plot density over time!')
print('')

for it in range(1128,len(time),10):
    plt.plot(time[1128:], densityFiltered[1128:],c='r')
    plt.axis([100000,370000,0.00004,0.00013])
    plt.gcf().subplots_adjust(left=0.15)
    plt.xlabel('iteration', fontsize=14)
    plt.ylabel('$\Delta \\rho$', fontsize=14)
    plt.scatter(time[it],densityFiltered[it],c='r')
    plt.grid(True)
    plt.savefig(os.path.join(currentDirectory,'densityOverTime/' + str(it) + '.png'))
    plt.close()

print('')
print('Plot energy over time!')
print('')

for it in range(1128,len(time),10):
    plt.plot(time[1128:], energy[1128:],c='r') 
    plt.axis([100000,370000,1.25,1.65])
    plt.xlabel('iteration', fontsize=14)
    plt.ylabel('$\Delta E_{kin}$', fontsize=14)
    plt.scatter(time[it], energy[it],c='r')
    plt.grid(True)
    plt.savefig(os.path.join(currentDirectory,'energyOverTime/' + str(it) + '.png'))
    plt.close()


#print(allData[0].data)
#print('')
#print(allData[8].data)
#plt.plot(allData[0].data,allData[8].data)
#plt.show()

# for data in allData:
#    plt.figure(data.time)
#    plt.hist(data.distributions,bins=numberOfBins,range=(minimum,maximum))
#    plt.axis([minimum,maximum,0,500])
#    plt.xlabel('$||\sigma_{i,j}||$', fontsize=14)
#    plt.ylabel('number of near-wall cells', fontsize=14)
#    plt.savefig(os.path.join(histogramDirectory,str(data.time).zfill(6) + '.png'))
#    plt.close()
# plt.show()
