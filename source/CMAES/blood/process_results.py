import sys
import numpy as np
import matplotlib.pyplot as plt

param = sys.argv[1]

print("running param " + param)
filename = "p" + param + "_results.txt"

print("opening file: " + filename)
f=open(filename, 'r')

lines=f.readlines()[3:] # skip first
result=[]

for x in lines:
    result.append(x.rstrip().split('\t')[1:])

f.close()

get = lambda idx : list(map(float, list(zip(*result))[idx]))

x1 = get(0)
x2 = get(1)
x3 = get(2)
x4 = get(3)
x5 = get(4)
x6 = get(5)
f  = get(6)

main = get(int(param)-1)

plt.subplot(711)
plt.plot(main, x1)
plt.ylim(-0.1,1.1)
plt.subplot(712)
plt.plot(main, x2)
plt.ylim(-0.1,1.1)
plt.subplot(713)
plt.plot(main, x3)
plt.ylim(-0.1,1.1)
plt.subplot(714)
plt.plot(main, x4)
plt.ylim(-0.1,1.1)
plt.subplot(715)
plt.plot(main, x5)
plt.ylim(-0.1,1.1)
plt.subplot(716)
plt.plot(main, x6)
plt.ylim(-0.02,0.22)
plt.subplot(717)
plt.plot(main, np.log(f))
#plt.ylim(-30,-25)

plt.show()
