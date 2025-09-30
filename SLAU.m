clc
clear
clearAllMemoizedCaches
size_ = input("Введите необходимую размерность матрицы коэффициентов: ");
fprintf("Размерность матрицы коэффициентов: %d\n", size_);
vectors = input("Введите необходимое количество векторов свободных членов: ");
fprintf("Количество векторов свободных членов: %d\n", vectors);

A = rand(size_, size_)*100
b = rand(vectors, size_)*100
fprintf("Выберите метод решения. [GS] для метода Гаусса или [LU] для LU-разложения \n");
choice = 0;
while choice == 0
  method = input("Введите код метода: ", 's');
  if method == "GS"
    choice = 1;
elseif method == "LU"
    choice = 1;
else
    disp("Неправильно выбран метод")
end
end
if method == "GS"
    disp("Вы выбрали метода Гаусса")
    tic
    result = _gauss_ (A, b);
    disp(result);
    elapsed_time = toc;

elseif method == "LU"
    disp("Вы выбрали метода LU-разложения")
    tic
    result = _lu_ (A, b);
    disp(result);
    elapsed_time = toc;

endif
fprintf("elapsed_time: %d\n", elapsed_time);

% Проверка A * X = b
disp("Проверка A * X = b:");
for i = 1:vectors
    computed_b = A * result(i,:)';
    original_b = b(i,:)';
    disp(original_b);
    disp(" ");
    disp(computed_b);
    error = norm(computed_b - original_b);
    fprintf("Вектор свободных членов %d: норма ошибки = %e\n", i, error);
    if error < 1e-6
        disp("Проверка пройдена: A * X ≈ b (в пределах погрешности).");
    else
        disp("Проверка не пройдена: A * X ≠ b.");
    endif
endfor

