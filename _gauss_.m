function X = _gauss_(factors, free)
  % if det(factors) == 0
  %     fprintf("Матрица вырожденная (det(A) = 0).\n");
  % endif
  B = free';
  A = factors;
  N = size(A)(1);
  for k = 1:N - 1
    for i = N:-1:k + 1
      t = A(i, k) / A(i - 1, k);
      A(i,:) -= A(i - 1,:) * t;
      B(i,:) -= B(i - 1,:) * t;
    endfor
  endfor
  X = NaN(size(B));
  X(N,:) = B(N,:) / A(N,N);
  for k = N - 1:-1:1
    summ = B(k,:);
    for i = k + 1:N
      summ -= X(i,:) * A(k, i);
    endfor
    X(k,:) = summ / A(k, k);
  endfor
  X = X';
end

