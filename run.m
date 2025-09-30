clc
clear
clearAllMemoizedCaches

functions = {@_gauss_, @_lu_};
f_names = {"Гаусса", "LU-разложения"};
hilbert_m_sizes = [10, 100, 1000];
m_sizes = [100, 1000, 2000, 3000, 4000, 5000, 10000];
v_sizes = [10, 100, 1000];

% Смотрим матрицы Гильберта
for n = 1:2
    for m_size = hilbert_m_sizes
        for v_size = v_sizes
            A = hilb(m_size);
            B = rand(v_size, m_size);

            fprintf("Метод %s...\n", f_names{n});
            fprintf("Вычисление матрицы Гильберта размерностью %d", m_size);
            fprintf(" для вектора свободных членов длиной %d...\n", v_size);

            if (ndims(A) != 2)
                disp("Массив коэффициентов не является матрицей.");
                continue;
            endif

            if (ndims(B) != 1 && ndims(B) != 2)
                disp("Свободные члены поданы неверно.");
                continue;
            endif

            if (max(size(A)) != min(size(A)))
                disp("Матрица коэффициентов не является квадратной.");
                continue;
            endif

            if (length(B(1, :)) != length(A))
                disp("Неверная длина вектора свободных членов.");
                continue;
            endif

            tic;
            result = functions{n}(A, B);
            fprintf("Продолжительность: %d\n", toc);

            disp("Проверка A * X = b...");
            error_count = 0;
            for i = 1:v_size
                computed_b = A * result(i,:)';
                original_b = B(i,:)';
                if norm(computed_b - original_b) > 1e-4
                    error_count++;
                endif
            endfor
            if error_count
                fprintf("Проверка не пройдена для %d векторов свободных членов.\n\n", error_count);
            else
                fprintf("Проверка пройдена: A * X ≈ b (в пределах погрешности).\n\n")
            endif
        endfor
    endfor
endfor

% Нормальные матрицы
for n = 1:2
    for m_size = m_sizes
        for v_size = v_sizes
            A = rand(m_size, m_size);
            B = rand(v_size, m_size);

            fprintf("Метод %s...\n", f_names{n});
            fprintf("Вычисление матрицы размерностью %d", m_size);
            fprintf(" для вектора свободных членов длиной %d...\n", v_size);

            if (ndims(A) != 2)
                disp("Массив коэффициентов не является матрицей.");
                continue;
            endif

            if (ndims(B) != 1 && ndims(B) != 2)
                disp("Свободные члены поданы неверно.");
                continue;
            endif

            if (max(size(A)) != min(size(A)))
                disp("Матрица коэффициентов не является квадратной.");
                continue;
            endif

            if (length(B(1, :)) != length(A))
                disp("Неверная длина вектора свободных членов.");
                continue;
            endif

            tic;
            result = functions{n}(A, B);
            fprintf("Продолжительность: %d\n", toc);

            disp("Проверка A * X = b...");
            error_count = 0;
            for i = 1:v_size
                computed_b = A * result(i,:)';
                original_b = B(i,:)';
                if norm(computed_b - original_b) > 1e-4
                    error_count++;
                endif
            endfor
            if error_count
                fprintf("Проверка не пройдена для %d векторов свободных членов.\n\n", error_count);
            else
                fprintf("Проверка пройдена: A * X ≈ b (в пределах погрешности).\n\n")
            endif
        endfor
    endfor
endfor
