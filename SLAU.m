clc
clear
clearAllMemoizedCaches

A_sizes = [100 1000 2000 3000 4000 5000 10000];
b_sizes = [10 100 1000];

% Запрос метода решения один раз
fprintf("Выберите метод решения. [GS] для метода Гаусса, [LU] для LU-разложения, [SD] для метода Зейделя или [SDR] для улучшенного метода Зейделя \n");
choice = 0;
while choice == 0
  method = input("Введите код метода: ", 's');
  if strcmp(method, "GS")
    choice = 1;
  elseif strcmp(method, "LU")
    choice = 1;
  elseif strcmp(method, "SD")
    choice = 1;
  elseif strcmp(method, "SDR")
    choice = 1;
  else
    disp("Неправильно выбран метод")
  end
end

for i = 1:length(A_sizes)
    sizeA = A_sizes(i);
    for j = 1:length(b_sizes)
        sizeB = b_sizes(j);

        fprintf('Матрица A: %d x %d, матрица b: %d x %d\n', sizeA, sizeA, sizeB, sizeA);

        % Генерируем данные
        A = rand(sizeA, sizeA)*100;
        b = rand(sizeB, sizeA)*100;

        tic;
        if strcmp(method, "GS")
            disp("Решаем методом Гаусса")
            result = _gauss_(A, b);
        elseif strcmp(method, "LU")
            disp("Решаем методом LU-разложения")
            result = _lu_(A, b);
        elseif strcmp(method, "SD")
            disp("Решаем методом Зейделя")
            result = _seidel_(A, b);
        elseif strcmp(method, "SDR")
            disp("Решаем улучшенным методом Зейделя")
            result = _seidel_robust_(A, b);
        end
        elapsed_time = toc;

        % Проверка результата
        if strcmp(method, "SD") || strcmp(method, "SDR")
            % Для методов Зейделя result уже имеет правильную размерность
            residual = A * result' - b';
        else
            % Для методов Гаусса и LU нужно транспонировать результат
            residual = A * result' - b';
        end
        error_norm = norm(residual);
        relative_error = error_norm / norm(b');

        fprintf('Время решения: %.4f секунд\n', elapsed_time);
        fprintf('Норма невязки: %e, Относительная ошибка: %e\n\n', error_norm, relative_error);
    end
end
