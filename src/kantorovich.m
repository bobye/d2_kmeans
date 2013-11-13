function [fval, x] = kantorovich(X, wX, Y, wY, x0)

  global A;
  global optim_options;
  

  
  n = size(X,2);
  m = size(Y,2);

  wX = wX/sum(wX);
  wY = wY/sum(wY);

  if (length(wX) ~= n || length(wY) ~= m ) 
     error('format not correct');
  end

  D = pdist2(X', Y', 'sqeuclidean');
  f = reshape(D, n*m, 1);
  Aeq = A{n,m}(1:end-1,:);
  beq = [wX'; wY(1:end-1)'];
  %Aeq = A{n,m};
  %beq = [wX'; wY'];

  if nargin == 4
      x0 = [];
  else
      x0 = reshape(x0,[n*m,1]);
  end
  
  [x, fval, exitflag] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), [], x0, optim_options );
  if exitflag < 0
      error('linprog no search direction [%d, %f]', exitflag, fval);
  end
  x = reshape(x, n, m);
end
