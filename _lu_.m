% Функция решает СЛАУ через LU-разложение.
% factors – матрица коэффициентов.
% free – матрица, каждая строка которой является вектором свободных членов.
% Возвращает матрицу, каждая строка которой является вектором решения
% или -1 при ошибке.
function retval = _lu_ (factors, free)
    retval = -1;

    if (ndims(factors) != 2)
        disp("Массив коэффициентов не является матрицей.");
        return;
    endif

    if (ndims(free) != 1 && ndims(free) != 2)
        disp("Свободные члены поданы неверно.");
        return;
    endif

    if (max(size(factors)) != min(size(factors)))
        disp("Матрица коэффициентов не является квадратной.");
        return;
    endif

    if (length(free(1, :)) != length(factors))
        disp("Неверная длина вектора свободных членов.");
        return;
    endif

    n = length(factors);

    % TODO: проверка на вырожденность: угловые миноры должны быть ненулевыми.

    % Заполнение первого столбца.
    %for i = 2:n
    %    factors(i, 1) = factors(i, 1) / factors(1, 1);
    %endfor

    % Заполнение остальных позиций матрицы.
    for i = 1:n

        for j = i + 1:n
            factors(j, i) = factors(j, i) / factors(i, i);

            for k = i + 1:n
                factors(j, k) = factors(j, k) - factors(j, i) * factors(i, k);
            endfor

        endfor

    endfor

    factors

    % Перебор векторов свободных членов.
    if length(free(:, 1)) > 1

        for b = transpose(free)
            y = zeros(1, n);

            for i = 1:n;
                minuend = 0;

                for j = 1:i - 1
                    minuend += factors(i, j) * y(j);
                endfor

                y(i) = b(i) - minuend;
            endfor

            for i = n:-1:1
                minuend = 0;

                for j = i + 1:n
                    minuend += factors(i, j) * x(j);
                endfor

                x(i) = (y(i) - minuend) / factors(i, i);
            endfor

            if retval == -1
                retval = x;
            else
                retval = vertcat(retval, x);
            endif

        endfor

    else
        y = zeros(1, n);
        x = zeros(1, n);

        for i = 1:n;
            minuend = 0;

            for j = 1:i - 1
                minuend += factors(i, j) * y(j);
            endfor

            y(i) = free(i) - minuend;
        endfor

        for i = n:-1:1
            minuend = 0;

            for j = i + 1:n
                minuend += factors(i, j) * x(j);
            endfor

            x(i) = (y(i) - minuend) / factors(i, i);
        endfor

        if retval == -1
            retval = x;
        else
            retval = vertcat(retval, x);
        endif

    endif

end

