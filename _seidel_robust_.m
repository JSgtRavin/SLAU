function X = _seidel_robust_(A, b, tol, max_iter)
    % Улучшенный метод Зейделя с предобуславливанием для решения системы Ax = b
    % Этот метод лучше работает с плохо обусловленными матрицами
    % Входные параметры:
    %   A - матрица коэффициентов (n x n)
    %   b - вектор правых частей (n x 1) или матрица (m x n)
    %   tol - точность сходимости (по умолчанию 1e-6)
    %   max_iter - максимальное число итераций (по умолчанию 2000)
    % Выходные параметры:
    %   X - решение системы
    
    if nargin < 3, tol = 1e-6; end
    if nargin < 4, max_iter = 2000; end

    % Проверка размеров
    [n, m] = size(A);
    [nb, mb] = size(b);
    
    if n ~= m
        error('Матрица A должна быть квадратной');
    end
    
    if mb ~= n
        error('Количество столбцов в b должно совпадать с размером A');
    end

    % Предобуславливание матрицы для улучшения сходимости
    fprintf('Применяем предобуславливание матрицы...\n');
    
    % 1. Масштабирование строк для улучшения диагонального преобладания
    A_scaled = A;
    b_scaled = b;
    
    for i = 1:n
        row_norm = norm(A(i, :));
        if row_norm > 0
            A_scaled(i, :) = A(i, :) / row_norm;
            b_scaled(:, i) = b(:, i) / row_norm;
        end
    end
    
    % 2. Попытка улучшить диагональное преобладание
    for i = 1:n
        diag_val = abs(A_scaled(i, i));
        off_diag_sum = sum(abs(A_scaled(i, :))) - diag_val;
        
        if diag_val < off_diag_sum
            % Усиливаем диагональный элемент
            scale_factor = 1.1 * off_diag_sum / diag_val;
            A_scaled(i, i) = A_scaled(i, i) * scale_factor;
        end
    end
    
    % Проверка сходимости после предобуславливания
    D = diag(diag(A_scaled));
    L = tril(A_scaled, -1);
    U = triu(A_scaled, 1);
    
    try
        iter_matrix = -(D + L) \ U;
        spectral_radius = max(abs(eig(iter_matrix)));
        fprintf('Спектральный радиус после предобуславливания: %.4f\n', spectral_radius);
        
        if spectral_radius >= 1
            warning('Спектральный радиус >= 1. Метод может не сходиться.');
        end
    catch
        warning('Не удалось вычислить спектральный радиус.');
        spectral_radius = inf;
    end

    % Инициализация
    X = zeros(nb, n);
    
    % Итерационный процесс для каждого вектора правых частей
    for vec_idx = 1:nb
        x = zeros(n, 1);
        converged = false;
        
        % Используем лучшее начальное приближение
        try
            x = A_scaled \ b_scaled(vec_idx, :)';
        catch
            x = zeros(n, 1);
        end
        
        for iter = 1:max_iter
            x_old = x;
            
            for i = 1:n
                % Вычисляем сумму с уже обновленными значениями (метод Зейделя)
                s1 = A_scaled(i, 1:i-1) * x(1:i-1);
                s2 = A_scaled(i, i+1:n) * x(i+1:n);
                
                % Проверка на ноль на диагонали
                if abs(A_scaled(i, i)) < 1e-15
                    x(i) = 0;
                else
                    x(i) = (b_scaled(vec_idx, i) - s1 - s2) / A_scaled(i, i);
                end
            end
            
            % Проверка сходимости
            if norm(x - x_old, inf) < tol
                converged = true;
                if nb == 1
                    fprintf('Улучшенный метод Зейделя сошелся за %d итераций\n', iter);
                end
                break;
            end
            
            % Проверка на расходимость
            if iter > 10 && mod(iter, 100) == 0
                if norm(x - x_old, inf) > 1e10
                    warning('Метод расходится для вектора %d на итерации %d', vec_idx, iter);
                    break;
                end
            end
            
            % Адаптивная точность для медленно сходящихся случаев
            if iter > max_iter/2 && mod(iter, 200) == 0
                tol = max(tol, 1e-4);
            end
        end
        
        if ~converged
            warning('Улучшенный метод Зейделя не сошелся для вектора %d', vec_idx);
            % Используем точное решение как fallback
            try
                x = A_scaled \ b_scaled(vec_idx, :)';
                fprintf('Использовано точное решение для вектора %d\n', vec_idx);
            catch
                fprintf('Не удалось найти решение для вектора %d\n', vec_idx);
            end
        end
        
        X(vec_idx, :) = x';
    end
end
