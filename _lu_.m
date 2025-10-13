% Функция решает СЛАУ через LU-разложение.
% factors – матрица коэффициентов.
% free – матрица, каждая строка которой является вектором свободных членов.
% Возвращает матрицу, каждая строка которой является вектором решения
% или -1 при ошибке.
function retval = _lu_(factors, free)
    retval = -1;

    % Проверка входных данных
    if (ndims(factors) ~= 2)
        disp("Массив коэффициентов не является матрицей.");
        return;
    end

    if (ndims(free) ~= 1 && ndims(free) ~= 2)
        disp("Свободные члены поданы неверно.");
        return;
    end

    if (max(size(factors)) ~= min(size(factors)))
        disp("Матрица коэффициентов не является квадратной.");
        return;
    end

    if (length(free(1, :)) ~= length(factors))
        disp("Неверная длина вектора свободных членов.");
        return;
    end

    n = length(factors);
    
    % Проверка на вырожденность
    if abs(det(factors)) < 1e-12
        disp("Матрица близка к вырожденной. Решение может быть неточным.");
    end

    % Создаем копию матрицы для LU-разложения
    LU = factors;

    % LU-разложение методом Гаусса
    for i = 1:n
        % Проверка на ноль на диагонали
        if abs(LU(i, i)) < 1e-12
            disp("Обнаружен нулевой или близкий к нулю диагональный элемент.");
            return;
        end
        
        for j = i + 1:n
            LU(j, i) = LU(j, i) / LU(i, i);
            
            for k = i + 1:n
                LU(j, k) = LU(j, k) - LU(j, i) * LU(i, k);
            end
        end
    end

    % Решение системы для каждого вектора свободных членов
    if length(free(:, 1)) > 1
        % Множественные правые части
        for b = transpose(free)
            y = zeros(1, n);
            x = zeros(1, n);

            % Прямой ход (Ly = b)
            for i = 1:n
                minuend = 0;
                
                for j = 1:i - 1
                    minuend = minuend + LU(i, j) * y(j);
                end
                
                y(i) = b(i) - minuend;
            end

            % Обратный ход (Ux = y)
            for i = n:-1:1
                minuend = 0;
                
                for j = i + 1:n
                    minuend = minuend + LU(i, j) * x(j);
                end
                
                x(i) = (y(i) - minuend) / LU(i, i);
            end

            if retval == -1
                retval = x;
            else
                retval = vertcat(retval, x);
            end
        end
    else
        % Одна правая часть
        y = zeros(1, n);
        x = zeros(1, n);

        % Прямой ход (Ly = b)
        for i = 1:n
            minuend = 0;
            
            for j = 1:i - 1
                minuend = minuend + LU(i, j) * y(j);
            end
            
            y(i) = free(i) - minuend;
        end

        % Обратный ход (Ux = y)
        for i = n:-1:1
            minuend = 0;
            
            for j = i + 1:n
                minuend = minuend + LU(i, j) * x(j);
            end
            
            x(i) = (y(i) - minuend) / LU(i, i);
        end

        retval = x;
    end
end

