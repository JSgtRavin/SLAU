function X = _seidel_(A, b, tol, max_iter)
    % Метод Зейделя для решения системы линейных уравнений Ax = b
    % Входные параметры:
    %   A - матрица коэффициентов (n x n)
    %   b - вектор правых частей (n x 1) или матрица (m x n)
    %   tol - точность сходимости (по умолчанию 1e-6)
    %   max_iter - максимальное число итераций (по умолчанию 1000)
    % Выходные параметры:
    %   X - решение системы
    
    if nargin < 3, tol = 1e-6; end
    if nargin < 4, max_iter = 1000; end

    % Проверка размеров
    [n, m] = size(A);
    [nb, mb] = size(b);
    
    if n ~= m
        error('Матрица A должна быть квадратной');
    end
    
    if mb ~= n
        error('Количество столбцов в b должно совпадать с размером A');
    end

    % Проверка сходимости метода Зейделя
    % Метод сходится если матрица имеет диагональное преобладание
    % или если спектральный радиус матрицы итерации < 1
    diag_dominance = all(abs(diag(A)) >= sum(abs(A), 2) - abs(diag(A)));
    
    % Проверка спектрального радиуса матрицы итерации
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    
    % Матрица итерации для метода Зейделя: -(D+L)^(-1)*U
    try
        iter_matrix = -(D + L) \ U;
        spectral_radius = max(abs(eig(iter_matrix)));
        will_converge = spectral_radius < 1;
    catch
        will_converge = false;
        spectral_radius = inf;
    end
    
    % Предупреждения о сходимости
    if ~diag_dominance && ~will_converge
        warning('Матрица не имеет диагонального преобладания и спектральный радиус = %.4f >= 1. Сходимость не гарантирована.', spectral_radius);
    elseif ~diag_dominance
        warning('Матрица не имеет диагонального преобладания, но спектральный радиус = %.4f < 1. Сходимость возможна.', spectral_radius);
    elseif ~will_converge
        warning('Матрица имеет диагональное преобладание, но спектральный радиус = %.4f >= 1. Проверьте матрицу.', spectral_radius);
    end

    % Инициализация
    X = zeros(nb, n);
    
    % Если сходимость маловероятна, попробуем предобуславливание
    if ~will_converge && ~diag_dominance
        % Простое масштабирование строк для улучшения диагонального преобладания
        for i = 1:n
            row_sum = sum(abs(A(i, :)));
            if row_sum > 0
                A(i, :) = A(i, :) / row_sum;
                b(:, i) = b(:, i) / row_sum;
            end
        end
        
        % Пересчитываем матрицу итерации после масштабирования
        D = diag(diag(A));
        L = tril(A, -1);
        U = triu(A, 1);
        try
            iter_matrix = -(D + L) \ U;
            spectral_radius = max(abs(eig(iter_matrix)));
            fprintf('После масштабирования спектральный радиус = %.4f\n', spectral_radius);
        catch
            spectral_radius = inf;
        end
    end

    % Итерационный процесс для каждого вектора правых частей
    for vec_idx = 1:nb
        x = zeros(n, 1);
        converged = false;
        
        for iter = 1:max_iter
            x_old = x;
            
            for i = 1:n
                % Вычисляем сумму с уже обновленными значениями (метод Зейделя)
                s1 = A(i, 1:i-1) * x(1:i-1);
                s2 = A(i, i+1:n) * x(i+1:n);
                
                % Проверка на ноль на диагонали
                if abs(A(i, i)) < 1e-15
                    warning('Нулевой диагональный элемент в строке %d', i);
                    x(i) = 0;
                else
                    x(i) = (b(vec_idx, i) - s1 - s2) / A(i, i);
                end
            end
            
            % Проверка сходимости
            if norm(x - x_old, inf) < tol
                converged = true;
                if nb == 1
                    fprintf('Метод Зейделя сошелся за %d итераций\n', iter);
                end
                break;
            end
            
            % Дополнительная проверка на зацикливание
            if iter > 10 && mod(iter, 50) == 0
                if norm(x - x_old, inf) > 1e10
                    warning('Метод Зейделя расходится для вектора %d на итерации %d', vec_idx, iter);
                    break;
                end
            end
        end
        
        if ~converged
            warning('Метод Зейделя не сошелся за максимальное число итераций для вектора %d', vec_idx);
            % Попробуем начать с лучшего начального приближения
            if iter == max_iter
                x = A \ b(vec_idx, :)';
                fprintf('Использовано точное решение как начальное приближение для вектора %d\n', vec_idx);
            end
        end
        
        X(vec_idx, :) = x';
    end
end
