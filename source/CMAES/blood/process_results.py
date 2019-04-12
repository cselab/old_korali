import matplotlib.pyplot as plt

f=open("p1_results.txt", 'r')

lines=f.readlines()[3:] # skip first
result=[]

for x in lines:
    result.append(x.rstrip().split('\t')[1:])

f.close()

x1 = list(zip(*result))[0]
x2 = list(zip(*result))[1]
x3 = list(zip(*result))[2]
x4 = list(zip(*result))[3]
x5 = list(zip(*result))[4]
x6 = list(zip(*result))[5]

print(result)
print(x1)

plt.plot(x2,x3,x4,x5,x6)
plt.axis(x1)
