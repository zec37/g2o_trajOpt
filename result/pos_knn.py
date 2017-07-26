import numpy as np
import matplotlib.pyplot as plt
from sklearn import neighbors

filename = "./optimized.txt"
optimized = open(filename)

X = []
label = []
n = 2 

for line in optimized:
    data = line.split()
    point = (float(data[1]), float(data[2]))
    X.append(point)
    label.append(data[3])

dataLen = len(X)
test = X[int(0.9 * dataLen):]
validate = label[int(0.9 * dataLen):]
clf = neighbors.KNeighborsClassifier(n, weights='uniform')
clf.fit(X, label)
result = clf.predict(test)

error_rate = (1 - np.sum(result == validate) / float(len(test))) 

print "\nTotal error rate is: %f" % error_rate

for i in range(0, len(test)):
    if(result[i] == validate[i]):
        plt.scatter(test[i][0], test[i][1], color = 'blue')
    else:
        plt.scatter(test[i][0], test[i][1], color = 'red')
plt.show()
