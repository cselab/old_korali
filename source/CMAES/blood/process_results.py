import matplotlib.pyplot as plt

f=open("p4_results.txt", 'r')

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

plt.subplot(611)
plt.plot(x1, x2)
plt.ylim(0,1)
plt.subplot(612)
plt.plot(x1, x3)
plt.ylim(0,1)
plt.subplot(613)
plt.plot(x1, x4)
plt.ylim(0,1)
plt.subplot(614)
plt.plot(x1, x5)
plt.ylim(0,1)
plt.subplot(615)
plt.plot(x1, x6)
plt.ylim(0,0.02)
plt.subplot(616)
plt.plot(x1, f)
plt.ylim(-30,-25)

plt.show()
