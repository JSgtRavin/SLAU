clc
clear
clearAllMemoizedCaches
size = input("Введите необходимую размерность матрицы коэффициентов: ");
fprintf("Размерность матрицы коэффициентов: %d\n", size);
vectors = input("Введите необходимуе количество векторов свободных членов: ");
fprintf("Количество векторов свободных членов: %d\n", vectors);

A = rand(size, size)*100
b = rand(vectors, size)*100
fprintf("Выберите метод решения. [GS] для метода Гаусса или [LU] для LU-разложения \n");
choise = 0;
while choise == 0
  method = input("Введите код метода: ", 's');
  if method == "GS"
    choise = 1;
elseif method == "LU"
    choise = 1;
else
    disp("Неправильно выбран метод")
end
end
if method == "GS"
    disp("Вы выбрали метода Гаусса")
    tic
    result = _gauss_ (A, b);
    elapsed_time = toc;

elseif method == "LU"
    disp("Вы выбрали метода LU-разложения")
    tic
    result = _lu_ (A, b)
    elapsed_time = toc;
    disp("elapsed_time")
end
