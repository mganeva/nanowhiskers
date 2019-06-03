import numpy as np
import matplotlib.pyplot as plt

rot0 = np.loadtxt('intensity1_000.txt', dtype=int)
rot5 = np.loadtxt('intensity1_005.txt', dtype=int)
rot10 = np.loadtxt('intensity1_010.txt', dtype=int)
rot15 = np.loadtxt('intensity1_015.txt', dtype=int)
rot20 = np.loadtxt('intensity1_020.txt', dtype=int)
rot25 = np.loadtxt('intensity1_025.txt', dtype=int)
rot30 = np.loadtxt('intensity1_030.txt', dtype=int)
rot35 = np.loadtxt('intensity1_035.txt', dtype=int)
rot40 = np.loadtxt('intensity1_040.txt', dtype=int)
rot45 = np.loadtxt('intensity1_045.txt', dtype=int)
rot50 = np.loadtxt('intensity1_050.txt', dtype=int)
rot55 = np.loadtxt('intensity1_055.txt', dtype=int)
rot60 = np.loadtxt('intensity1_060.txt', dtype=int)

output_rot0 = [sum([i[b] for i in rot0]) for b in range(len(rot0[0]))]
output_rot5 = [sum([i[b] for i in rot5]) for b in range(len(rot5[0]))]
output_rot10 = [sum([i[b] for i in rot10]) for b in range(len(rot10[0]))]
output_rot15 = [sum([i[b] for i in rot15]) for b in range(len(rot15[0]))]
output_rot20 = [sum([i[b] for i in rot20]) for b in range(len(rot20[0]))]
output_rot25 = [sum([i[b] for i in rot25]) for b in range(len(rot25[0]))]
output_rot30 = [sum([i[b] for i in rot30]) for b in range(len(rot30[0]))]
output_rot35 = [sum([i[b] for i in rot35]) for b in range(len(rot35[0]))]
output_rot40 = [sum([i[b] for i in rot40]) for b in range(len(rot40[0]))]
output_rot45 = [sum([i[b] for i in rot45]) for b in range(len(rot45[0]))]
output_rot50 = [sum([i[b] for i in rot50]) for b in range(len(rot50[0]))]
output_rot55 = [sum([i[b] for i in rot55]) for b in range(len(rot55[0]))]
output_rot60 = [sum([i[b] for i in rot60]) for b in range(len(rot60[0]))]

x = [i for i in range(0,256)]

fig = plt.figure()
graph1 = plt.plot(x, output_rot60, color=(0, 0, 0.6), label = "rot60")
#graph2 = plt.plot(x, output_rot5, color=(0, 0, 0.9), label = "rot5")
#graph3 = plt.plot(x, output_rot10, color=(0, 0.6, 0), label = "rot10")
#graph4 = plt.plot(x, output_rot15, color=(0, 0.9, 0), label = "rot15")
#graph5 = plt.plot(x, output_rot20, color=(0.6, 0, 0), label = "rot20")
#graph6 = plt.plot(x, output_rot25, color=(1, 0, 0), label = "rot25")
plt.legend ()
plt.savefig('Sum_rot60.png')
#print('Plot: ', len(graph1), graph1)
#print('Plot: ', len(graph2), graph2)
#print('Plot: ', len(graph3), graph3)
#print('Plot: ', len(graph4), graph4)
#print('Plot: ', len(graph5), graph5)
#print('Plot: ', len(graph6), graph6)

plt.show()
