function X = _gauss_(factors, free)
  if det(factors) == 0
      error('Матрица вырожденная (det(A) = 0), нужна обратимая.');
  end
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
    X(k,:) = (B(k,:) - sum(X(k + 1: N,:) .* A(k,:)(k + 1:N)')) / A(k,k);
  endfor
  X = X';
end
