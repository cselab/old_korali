import matplotlib.pyplot as plt

f=open("p1_results.txt", 'r')

lines=f.readlines()[3:] # skip first
result=[]

for x in lines:
    result.append(x.rstrip().split('\t')[1:])

f.close()

get = lambda idx : list(map(float, list(zip(*result))[idx]))

print(result)

x1 = get(0)
x2 = get(1)
x3 = get(2)
x4 = get(3)
x5 = get(4)
x6 = get(5)
f  = get(6)

plt.plot(x1,x2,'--b')
plt.show()
