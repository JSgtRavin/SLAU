clc
clear
clearAllMemoizedCaches

matrix_size = input("Введите необходимую размерность матрицы коэффициентов: ");
fprintf("Размерность матрицы коэффициентов: %d\n", matrix_size);
vector_size = input("Введите необходимое количество векторов свободных членов: ");
fprintf("Количество векторов свободных членов: %d\n", vector_size);

A = rand(matrix_size, matrix_size) * 100;
b = rand(vector_size, matrix_size) * 100;

fprintf("Выберите метод решения. [GS] для метода Гаусса или [LU] для LU-разложения \n");

choice = 0;
while choice == 0
    method = input("Введите код метода: ", 's');
    if method == "GS"
        choice = 1;
    elseif method == "LU"
        choice = 1;
    else
        disp("Неправильно выбран метод");
    endif
end

if method == "GS"
    disp("Вы выбрали метод Гаусса.\nВычисления...");
    tic;
    result = _gauss_ (A, b);
    elapsed_time = toc;
elseif method == "LU"
    disp("Вы выбрали метод LU-разложения.\nВычисления...");
    tic;
    result = _lu_ (A, b);
    elapsed_time = toc;
endif
fprintf("Прошло времени: %d\n", elapsed_time);

disp("Проверка A * X = b...");
for i = 1:vector_size
    computed_b = A * result(i,:)';
    original_b = b(i,:)';
    error = norm(computed_b - original_b);
    if error > 1e-4
        fprintf("Вектор свободных членов %d: норма ошибки = %e\n", i, error);
        disp("Проверка не пройдена: A * X ≠ b.");
    % else
    %     disp("Проверка пройдена: A * X ≈ b (в пределах погрешности).");
    endif
endfor
