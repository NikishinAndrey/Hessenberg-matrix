import numpy as np
import random as r
import matplotlib.pyplot as plt

length = 10  # Размерность матрицы
length3 = 10  # Количество коэффициентов отделимости чисел
length4 = 10  # Количество дефформируемых матриц


# Создаем диагональную матрицу
def do_matrix(koef):
    mat_c = [[0 for i in range(length)] for j in range(length)]
    for i in range(length):
        mat_c[i][i] = (i * koef + 1)
    mat_c = np.array(mat_c)  # Диагональная матрица
    a__ = np.array([[r.randint(1, 100) for i in range(length)] for j in range(length)])  # Создаем рандомную матрицу
    a_ = np.array(np.linalg.inv(a__))  # Находим обратную матрицу
    matrix = a__.dot(mat_c).dot(a_)  # Получаем конечную симметричную матрицу
    return matrix


# eig_value1 = np.array(result1[1])  # Собственные числа матрицы
# eig_vector1 = np.array(result1[2])  # Матрица с точным значением собственных векторов
# print('Заданые собственные числа: ')
# print(eig_value1)
# print('Заданные собственные вектора: ')
# print(eig_vector1)


def deforming_matrix(a):
    massive_a = [[[0 for i in range(length)] for j in range(length)] for k in range(length4)]
    for k in range(length4):
        for i in range(length):
            for j in range(length):
                massive_a[k][i][j] = a[i][j] * (1 + 0.01 * k)
    return massive_a

# Реализуем метод приведения матрицы к специальному виду(Хессенберга)


def method_special_matrix(a):
    for k in range(length - 1):
        w_ = [[0 for i in range(length)] for j in range(length)]
        w = [0 for j in range(length)]
        summa = 0
        for j in range(k + 1, length):
            summa += a[j][k] ** 2
        s = np.sign(-a[k + 1][k]) * (summa ** 0.5)
        m = 1 / ((2 * s * (s - a[k + 1][k])) ** 0.5)
        w[k + 1] = m * (a[k + 1][k] - s)
        for int in range(k + 2, length):
            w[int] = m * a[int][k]
        for i in range(length):
            w_[i][i] += 1
            for int2 in range(length):
                w_[i][int2] += -2 * w[i] * w[int2]
        P = np.array(w_)
        a = P.dot(a).dot(P)
    return a

# Получили матрицу Хессенберга


# Считаем след матрицы, квадрата матрицы, куба матрицы... для нахождения
# в дальнейшем собственных чисел


def trace_matrix(a):
    matrix_a = [0 for i in range(length - 1)]
    matrix_a[0] = np.trace(a)
    b = a
    for i in range(1, length - 1):
        a = a.dot(b)
        matrix_a[i] = np.trace(a)
    return matrix_a


# Получили значение следа для матрицы n-ой степени
# Определяем длину для массива,


# заполненного значениями следа для матрицы n-ой степени


# Считаем коэффициенты в характеристическом многочлене


def coefficients_of_character_polynomial(res, l):
    p = [0 for i in range(l)]
    p[0] = res[0]
    for i in range(1, l):
        p[i] = res[i] / (i + 1)
        for j in range(i):
            p[i] -= p[j] * res[i - 1 - j] / (i + 1)
    return p

# Сохраняем результат предыдущей функции


# Досчитываем коэффициенты в характеристическом многочлене

def do_polynomial(res, mat, l):
    k = [-1 for i in range(l)]
    for j in range(1, l):
        k[j] = res[j - 1]
    k.append(((-1) ** (l - 1)) * np.linalg.det(mat))
    return k


# Сохраняем окончательный результат с коэффициентами
# Получаем корни характерестического уравнения
# Получаем собственные числа с помощью решения уравнения n-ой степени
# print('Получившиеся собственные числа: ')
# print(eig_value2)


# Создаем систему из (n-1)-ого уравнений, закладывая x(n) = 1
def character_polynomial(mat, eig_v):
    a = [[[0 for i in range(length - 1)] for j in range(length - 1)] for k in range(length)]
    b = [[0 for i in range(length - 1)] for k in range(length)]
    for k in range(length):
        for i in range(length - 1):
            b[k][i] = -mat[i][length - 1]
            for j in range(length - 1):
                a[k][i][j] = mat[i][j]
            a[k][i][i] = mat[i][i] - eig_v[k]
    return a, b


# Сохраняем результат функции
# Массив матриц для разных собственных чисел
# Массив свободных векторов для разных собственных чисел


# Находим собственные вектора для каждого собственного числа


def find_eigenvectors(m_a, m_b):
    matrix_x = [[0 for i in range(length - 1)] for j in range(length)]  # Массив собственных векторов
    for k in range(length):
        matrix_x[k] = np.linalg.solve(m_a[k], m_b[k]).tolist()
        matrix_x[k].append(1)
        matrix_x[k] = matrix_x[k] / np.linalg.norm(matrix_x[k])
    return matrix_x


# Собственные вектора для соответсвующих собственных чисел
# print('Получившиеся собственные вектора: ')
# print(eig_vector2)


# Находим значения синуса между заданным и получившимся собственнымыми векторами


def sin_phi(eig_v1, eig_v2):
    sin_phi = [0 for i in range(length)]
    for i in range(length):
        phi = np.arccos((eig_v1[i].dot(eig_v2[i])))
        sin_phi[i] = np.sin(phi)
    sin_phi_ = min(sin_phi)
    return sin_phi_


# print("Значение синуса: ")
# print(relative_error)


def chart1():
    coefficients = [0 for i in range(length3)]
    sinus = [0 for i in range(length3)]
    eig_vector1 = np.eye(length)
    rel = [0 for i in range(length3)]
    for j in range(length3):
        coefficients[j] = 1 / (10 ** j)
        mat_cf = do_matrix(coefficients[j])
        mat_hes = method_special_matrix(mat_cf)
        mat_tr = trace_matrix(mat_hes)
        length2 = len(mat_tr)
        mat_ch = coefficients_of_character_polynomial(mat_tr, length2)
        res1 = do_polynomial(mat_ch, mat_hes, length)
        roots = np.roots(res1)
        eig_value2 = roots[:: -1]
        res2 = character_polynomial(mat_cf, eig_value2)
        mat_a = res2[0]
        mat_b = res2[1]
        eig_vector2 = find_eigenvectors(mat_a, mat_b)
        sinus[j] = sin_phi(eig_vector1, eig_vector2)
        rel[j] = 1 - (sinus[j] / sinus[0])
    return coefficients, sinus, rel

result3 = chart1()
coefficient_ = result3[0][::-1]
sinus_ = result3[1]
chart2 = result3[2]


plt.figure(1)
plt.grid()
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Отделимость собственных чисел')
plt.ylabel('Максимальный синус угла между собственными векторами')
plt.title("График зависимости максимального синуса угла между собственными векторами от отделимости собственных чисел")
plt.plot(coefficient_, sinus_)

def chart2():
    eig_vector1 = np.eye(length)
    rel = [0 for i in range(length4)]
    for j in range(length4):
        mat_cf = do_matrix(1)
        mat_dm = np.array(deforming_matrix(mat_cf))
        mat_hes = method_special_matrix(mat_dm[j])
        mat_tr = trace_matrix(mat_hes)
        length2 = len(mat_tr)
        mat_ch = coefficients_of_character_polynomial(mat_tr, length2)
        res1 = do_polynomial(mat_ch, mat_hes, length)
        roots = np.roots(res1)
        eig_value2 = roots[:: -1]
        res2 = character_polynomial(mat_cf, eig_value2)
        mat_a = res2[0]
        mat_b = res2[1]
        eig_vector2 = find_eigenvectors(mat_a, mat_b)
        sinus = sin_phi(eig_vector1, eig_vector2)
        if j != 0:
            rel[j] = abs(1 - (sinus / sinus_0))
        else:
            sinus_0 = sinus
    return rel



error = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

plt.figure(2)
plt.grid()
plt.xlabel('Изменение исходной матрицы в процентах')
plt.ylabel('Относительная погрешность относительно исходной матрицы')
plt.title("График зависимости относительной погрешности от изменения исходной матрицы")
plt.plot(error, chart2(), label='Отделимость собсвтенных чисел = 1')
plt.legend()
plt.show()