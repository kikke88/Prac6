import sys
with open(sys.argv[1]) as f:
    array = []
    for line in f:
        array.append([float(x) for x in line.split()])
arr = []
for y in array:
    arr.append(y[0])
array = []
arr = sorted(arr)

sum = 0.0
import numpy as np
x = np.array(arr)
for y in x:
    sum += y
print(sum / x.size)


'''
print(arr)
import matplotlib.pyplot as plt
import numpy as np
x = np.array(arr)
bins = np.linspace(x[0], x[-1], x.size / 2)
plt.hist(x, bins, alpha = 0.7, histtype='bar', ec='black')
plt.xticks(fontsize = 7)
#plt.ylim(0, 5)
plt.grid(color = 'k', linestyle = '-', linewidth = 0.2)
plt.title('Распределение 1 - F, ' + sys.argv[1][2:4] + ' кубитов, ' + str(x.size) + ' запусков, eps = ' + sys.argv[2])
plt.show()
'''
