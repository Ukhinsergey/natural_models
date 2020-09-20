# Задание 1. Методы Монте-Карло

Требуется написать параллельную программу, реализующую метод Монте-Карло для модели случайных блужданий.
Программа должна принимать на вход (через аргументы командной строки) пять пара-метров:

* a, b — границы интервала,
* x — начальная позиция,
* p — вероятность перехода частицы вправо,
* N — число частиц.
Результатом работы программы должны быть два файла. В файле output.txt должны
быть записаны полученные данные (вероятность достижения состояния b и среднее время

жизни одной частицы t). В файле stat.txt должно быть сохранено время T работы ос-
новного цикла программы с указанием всех параметров запуска (см. выше), в том числе
параметр P — число параллельных процессов.

После того, как программа протестирована, необходимо выполнить следующие серии рас-
четов, результаты которых следует оформить в виде графиков:

1. T(N) (время), S(N) (ускорение) и E(N) (эффективность) при фиксированном значении P.
2. T(P), S(P) и E(P) при фиксированном значении N.
3. T(P), S(P) и E(P) при условии N = 103P.
