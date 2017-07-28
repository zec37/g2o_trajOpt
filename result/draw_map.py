import matplotlib.pyplot as plt

file_aligned = './aligned.txt'
file_optimized = './optimized.txt'
aligned = open(file_aligned)
optimized = open(file_optimized)

x_a = []
y_a = []
x_o = []
y_o = []

for line in aligned:
    data = line.split()
    x_a.append(float(data[1]))
    y_a.append(float(data[2]))
for line in optimized: 
    data = line.split()
    x_o.append(float(data[1]))
    y_o.append(float(data[2]))


plt.subplot(211)
plt.scatter(y_a, x_a)
plt.subplot(212)
plt.scatter(y_o, x_o)
plt.show()
