import numpy as np


# Ядерные функции. По умолчанию используется первая.
def nuclearFunction(z, x=1, r=2, s=3):

    p = 0

    if (x == 1):
        p = (1-(z**r))**s

    elif (x == 2):
        p = np.exp(-s*(z**r))

    elif (x == 3):
        p = 1/(s*(z**r)+1)

    elif (x == 4):
        p = 1/(z**s)

    return(p)


def f(x, y):

    def z(x): return -(1/((x-1)**2 + 0.2)) - \
        (1/(2*(x-2)**2 + 0.15)) - (1/(3*(x-3)**2 + 0.3))

    return z(x)+z(y)


def ruban_algorithm(x, deltaX, f, lower, upper, n=500, e=0.001, M=1000, y=1, q=2, numberOfNuclearFunc=1, r=2, s=100):
    if (len(x) != len(deltaX) or len(x) != len(lower) or len(x) != len(upper)):
        print("Ошибка, не соблюдение размерности у входящих параметров (x, deltaX, upper, lower)!")
        return
    if (n < 0):
        print("Ошибка, количество пробных точек должно быть больше 0!")
        return
    if (e < 0):
        print("Ошибка, константа точности не должна быть меньше 0!")
        return
    if (M < 0):
        print("Ошибка, максимальное количество итераций не должно быть меньше 0!")
        return
    if (y < 0):
        print("Ошибка, (y) не должен быть меньше 0!")
        return

    k = 1  # Будем использовать для проверки первого критерия останова - превышения числа итераций
    # В матрицу будут записываться пробные точки.
    testX = np.zeros((n, len(x)))

    # Где n - количество пробных точек,
    # length(x) - количество координатых направлений
    # А эта матрица содержит значения функции в пробных точках
    functionValues = np.zeros(n)
    # В этой матрице будут храниться значения u для генерации пробных точек
    uValues = np.zeros((n, len(x)))
    p = np.zeros(n)  # Ядра
    pNorm = np.zeros(n)  # нормированные ядра
    allX = np.zeros((M, len(x)))  # Сюда записываем все X
    allX[0] = x
    # Сюда записываем все результаты. Если число циклов превысит M, то мы выберем из всех результатов наименьший
    allResults = [f(x[0], x[1])]
    while(k < M):
        # Шаг первый - инициализация пробных точек и расчет значений функции в них:
        for i in range(n):
            for j in range(len(x)):
                # Генерируем u от -1 до 1
                uValues[i, j] = np.random.rand()*2 - 1
            testX[i] = x + deltaX * uValues[i, ]  # Генерируем пробную точку
            # Находит значение функции в этой точке(f - какая-то функция)
            functionValues[i] = f(testX[i, 0], testX[i, 1])

        gmin = np.zeros(n)
        # Шаг второй - рассчитываем ядра и нормированные ядра:
        for i in range(n):
            gmin[i] = (functionValues[i] - min(functionValues)) / \
                (max(functionValues) - min(functionValues))
        for i in range(n):  # Вычисляем ядра
            p[i] = nuclearFunction(x=numberOfNuclearFunc, z=gmin[i], r=r, s=s)
        for i in range(n):  # Вычисляем нормированные ядра
            pNorm[i] = p[i]/sum(p)
        # Шаг третий - вычисляем новый центр прямоугольника и его длину:
        for i in range(len(x)):
            x[i] = x[i] + deltaX[i] * \
                sum(map(lambda x: uValues[x, i]*pNorm[x], range(n)))
        print("---------------")
        print("Результат ", k, "итерации: x=", x, "\n")
        if (x[0] < lower[0]):
            x[0] < -lower[0]
        if (x[1] < lower[1]):
            x[1] < -lower[1]
        if (x[0] > upper[0]):
            x[0] < -upper[0]
        if (x[1] > upper[1]):
            x[1] < -upper[1]
        for i in range(len(x)):
            deltaX[i] = y*deltaX[i] * \
                ((sum(
                    map(lambda x: (abs(uValues[x, i])**(q)) * pNorm[x], range(n))))**(1/q))
        print("deltaX=", deltaX, "\n")
        k = k+1  # Увеличиваем счетчик итераций
        allX[k] = x
        allResults.append(f(x[0], x[1]))
        if (max(deltaX) <= e):
            break
        # if (sqrt(sum(deltaX**2)) < e)
        #  break

    print("Сделано шагов: ", k, "\n")
    return allX[allResults.index(min(allResults))]


def f(x, y):
    z = [7*(abs(x)**2) + 7*(abs(y)**2)]
    z.append(5*(abs(x-3))**0.8 + 5*(abs(y-3)**0.6) + 6)
    z.append(5*(abs(x-6))**1.3 + 5*(abs(y-6)**1.3) + 2)
    z.append(5*(abs(x-6))**1 + 5*(abs(y+6)**1) + 8)
    z.append(4*(abs(x+6))**1.5 + 4*(abs(y+6)**1.5) + 7)
    z.append(5*(abs(x+3))**1.8 + 5*(abs(y)**1.8) + 9)
    z.append(6*(abs(x+6))**0.6 + 6*(abs(y-6)**0.9))
    return(min(z))


print("Минимум в точке (0;0)")
print(f(0, 0))
x = ruban_algorithm(x=[10, 10], deltaX=[20, 20], lower=[-10, -10], upper=[30, 30],
                    f=f, n=500, y=1, q=2, s=100, e=0.0001, r=2, numberOfNuclearFunc=1)
print(x)
