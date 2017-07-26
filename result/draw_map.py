import matplotlib.pyplot as plt

file_aligned = './aligned.txt'
file_optimized = './optimized.txt'
aligned = open(file_aligned)
optimized = open(file_optimized)

x_a = []
y_a = []
x_o = []
y_o = []
error = []

for line in aligned:
    data = line.split()
    x_a.append(float(data[1]))
    y_a.append(float(data[2]))
for line in optimized: 
    data = line.split()
    x_o.append(float(data[1]))
    y_o.append(float(data[2]))
    error.append(float(data[5]))

plt.subplot(311)
plt.scatter(y_a, x_a)
plt.subplot(312)
plt.scatter(y_o, x_o)
plt.subplot(313)
plt.errorbar(y_o, x_o, yerr=error, fmt = 'o')
plt.show()
