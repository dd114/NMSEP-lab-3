import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

if __name__ == '__main__':

    f = open('trueSolution.txt', 'r') # требования к файлу: первая колонка - время, вторая - координата, третья - значения функции на них
    tempT = list()
    tempX = list()
    tempU = list()

    for line in f:
        tempValue = [float(i) for i in line.strip().split()]
        tempT += [tempValue[0]]
        tempX += [tempValue[1]]
        tempU += [tempValue[2]]
    f.close()

    tempT = np.array(tempT)
    tempX = np.array(tempX)
    tempU = np.array(tempU)

    print(tempT)
    print(tempX)
    print(tempU)

    # t = np.linspace(tempT.min(), tempT.max(), len(np.unique(tempT)))
    # x = np.linspace(tempX.min(), tempX.max(), len(np.unique(tempX)))

    xgrid, tgrid = np.meshgrid(np.unique(tempX), np.unique(tempT))

    ax = plt.axes(projection='3d')
    ax.plot_wireframe(xgrid, tgrid, tempU.reshape((len(np.unique(tempT)), len(np.unique(tempX)))))
    # ax.plot(tgrid, xgrid, tempU.reshape((len(np.unique(tempT)), len(np.unique(tempX)))))
    plt.xlabel("ось X")
    plt.ylabel("ось T")
    plt.show()