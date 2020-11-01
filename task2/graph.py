import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('best2.txt')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(data, label =  "Наименьшая ошибка")
data2 = np.loadtxt('av2.txt')
plt.plot(data2, label =  "Средняя ошибка")

plt.title('')
plt.xlabel('Номер итерации')
plt.ylabel('Величина ошибки')
plt.legend()
plt.show()