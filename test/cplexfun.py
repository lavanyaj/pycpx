import logging
from pycpx import CPlexModel
import numpy as np

m = CPlexModel(verbosity=3)

sizeItems = [1,4,3,2,1,3, 4, 1, 4, 5, 2, 5, 6, 2, 8, 4, 3, 4, 1, 1, 1, 2,\
                 7, 6, 5, 5, 3, 2, 2, 1, 6, 5, 7, 5, 4, 3, 3, 2, 2, 1, 2, 3]
sizeBin = max(sizeItems)


numItems = len(sizeItems)
numBins = numItems

assign = m.new((numItems,numBins), vtype=int, lb=0, ub=1,\
                                  name='assign')

binUsed = m.new((numBins), vtype=int, lb=0, ub=1, name='binUsed')

startAssign = np.zeros((numItems,numBins))
startAssign[20,8] = 1
startAssign[32,8] = 1
startAssign[30,9] = 1
startAssign[38,9] = 1
startAssign[11,11] = 1
startAssign[26,11] = 1
startAssign[0,14] = 1
startAssign[4,14] = 1
startAssign[7,14] = 1
startAssign[29,14] = 1
startAssign[39,14] = 1
startAssign[41,14] = 1
startAssign[1,17] = 1
startAssign[37,17] = 1
startAssign[40,17] = 1
startAssign[15,19] = 1
startAssign[34,19] = 1
startAssign[3,20] = 1
startAssign[10,20] = 1
startAssign[17,20] = 1
startAssign[21,22] = 1
startAssign[23,22] = 1
startAssign[5,25] = 1
startAssign[33,25] = 1
startAssign[9,27] = 1
startAssign[28,27] = 1
startAssign[6,29] = 1
startAssign[13,29] = 1
startAssign[27,29] = 1
startAssign[12,33] = 1
startAssign[18,33] = 1
startAssign[19,33] = 1
startAssign[31,34] = 1
startAssign[36,34] = 1
startAssign[2,35] = 1
startAssign[24,35] = 1
startAssign[25,36] = 1
startAssign[35,36] = 1
startAssign[22,37] = 1
startAssign[14,38] = 1
startAssign[8,39] = 1
startAssign[16,39] = 1


for b in range(numBins):
    totalSize = sum([assign[i,b] * sizeItems[i] for i in range(numItems)])
    m.constrain(totalSize <= sizeBin)
    upperBound = sum(sizeItems)
    m.constrain(totalSize <= binUsed[b] * upperBound)
    m.constrain(binUsed[b] <= totalSize)
    pass

# each item assigned to exactly one bin
for i in range(numItems):
    numBinsForItem = sum([assign[i,b] for b in range(numBins)])
    m.constrain(numBinsForItem == 1)
    pass

numBinsUsed = sum(binUsed)

try:
    m.minimize(numBinsUsed, starting_dict={assign:startAssign}, emphasis=1,\
                   time_limit=100)
    pass
except Exception, e:
    logging.exception(e)
    pass

print "Number of constraints: %.1f" % m.getNRows()
print "Number of variables: %.1f" % m.getNCols()
print "Number of quadratic constraints: %.1f" % m.getNQCs()

m.minimize(numBinsUsed, starting_dict={assign:startAssign}, emphasis=1,\
               tree_limit=2, work_dir="mapper/output")


print m[binUsed]
print int(m[numBinsUsed])

for b in range(numBins):
    print "\nIn bin " + str(b)
    totalSize = sum([m[assign][i,b] * sizeItems[i] for i in range(numItems)])
    print "Total size: " + str(totalSize)
    for i in range(numItems):
        if m[assign][i,b] == 1:
            print " Item " + str(i) + " with size " + str(sizeItems[i]),
            pass
        pass
    pass

for b in range(numBins):
    for i in range(numItems):
        if m[assign][i,b] == 1:
            print "assign["+str(i)+","+str(b)+"] = 1"
            pass
        pass
    pass
