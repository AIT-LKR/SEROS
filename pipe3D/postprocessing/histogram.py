import numpy as np
from os import listdir
from os.path import isfile, join
import os
import matplotlib.pyplot as plt

class dataSet:
    histogram = []

    def __init__(self, time, distributions):
        self.time = time
        self.distributions = distributions

def mysplit(string):
    head = string.rstrip('0123456789')
    tail = string[len(head):]
    return tail

maximum = 0
minimum = 1
numberOfBins = 20
minTimeStep = 10000

currentDirectory = os.path.dirname(os.path.realpath(__file__)) + '/distribution'
histogramDirectory = os.path.dirname(os.path.realpath(__file__)) + '/histogram' 

print('')
print('Data directory:')
print('---------------')
print(currentDirectory)

distributionFiles = [f for f in listdir(currentDirectory) if (isfile(join(currentDirectory, f)) and f.endswith('.dat'))]

print('')
print('Read ', len(distributionFiles), ' distribution files!')


allData = []

for distributionFile in distributionFiles:
    currentFile = os.path.join(currentDirectory, distributionFile)
    readData = np.genfromtxt(currentFile)
    fileName = os.path.splitext(os.path.basename(currentFile))
    time = mysplit(fileName[0]) 
    if int(time) >= minTimeStep:
        data = dataSet(int(time),readData)    
        currentMax = max(data.distributions)
        currentMin = min(data.distributions)
        if maximum < currentMax:
            maximum = currentMax 
        if minimum > currentMin:
            minimum = currentMin
        allData.append(data)
#        print('Data for timestep ', time, ' added!')

allData.sort(key=lambda x: x.time)

print('')
print('Number of data sets:', len(allData), ' added!')
print('')
print('Viscous stresses:')
print('-----------------')
print('Maximum is:', maximum)
print('Minimum is:', minimum)
print('')
#print('Calculate histograms!')
#print('')
#for data in allData:
#     histogram = np.histogram(data.distributions, numberOfBins, range=(minimum,maximum))
#    data.histogram = histogram
#    print(histogram[0])

# print(allData[100].histogram[0], allData[100].histogram[1])

print('Plot histograms!')
print('')

for data in allData:
    plt.figure(data.time)
    plt.hist(data.distributions,bins=numberOfBins,range=(minimum,maximum))
    plt.axis([minimum,maximum,0,500])
    plt.xlabel('$||\sigma_{i,j}||$', fontsize=14)
    plt.ylabel('number of near-wall cells', fontsize=14)
    plt.savefig(os.path.join(histogramDirectory,str(data.time).zfill(6) + '.png'))
    plt.close()
# plt.show()
print('Histograms are plotted!')
print('')
