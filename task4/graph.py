import numpy as np
import matplotlib.pyplot as plt


#ax_X = np.array([1, 2, 4, 8, 16, 32, 64])

#p = 1, T = 100
ax_X = np.arange(0, 3000)

ax_Y0 = np.loadtxt('2stat0.txt')
ax_Y1 = np.loadtxt('2stat1.txt')
# ax_Y2 = np.loadtxt('2stat2.txt')
# ax_Y3 = np.loadtxt('2stat3.txt')
# ax_Y4 = np.loadtxt('8stat4.txt')
# ax_Y5 = np.loadtxt('8stat5.txt')
# ax_Y6 = np.loadtxt('8stat6.txt')
# ax_Y7 = np.loadtxt('8stat7.txt')
# ax_Y8 = np.loadtxt('16stat8.txt')
# ax_Y9 = np.loadtxt('16stat9.txt')
# ax_Y10 = np.loadtxt('16stat10.txt')
# ax_Y11 = np.loadtxt('16stat11.txt')
# ax_Y12 = np.loadtxt('16stat12.txt')
# ax_Y13 = np.loadtxt('16stat13.txt')
# ax_Y14 = np.loadtxt('16stat14.txt')
# ax_Y15 = np.loadtxt('16stat15.txt')

fig = plt.figure()

plt.plot(ax_X, ax_Y0)
plt.plot(ax_X, ax_Y1)
# plt.plot(ax_X, ax_Y2)
# plt.plot(ax_X, ax_Y3)
# plt.plot(ax_X, ax_Y4)
# plt.plot(ax_X, ax_Y5)
# plt.plot(ax_X, ax_Y6)
# plt.plot(ax_X, ax_Y7)
# plt.plot(ax_X, ax_Y8)
# plt.plot(ax_X, ax_Y9)
# plt.plot(ax_X, ax_Y10)
# plt.plot(ax_X, ax_Y11)
# plt.plot(ax_X, ax_Y12)
# plt.plot(ax_X, ax_Y13)
# plt.plot(ax_X, ax_Y14)
# plt.plot(ax_X, ax_Y15)
plt.xscale('log')

plt.xlabel('Номер итерации')
plt.ylabel('Загрузка процесса')
plt.title('  w->a, a->ab [0.01], b->a [0.01]');

plt.legend()

plt.show()

