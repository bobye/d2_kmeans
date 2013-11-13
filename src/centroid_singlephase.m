function [c] = centroid_singlephase(stride, supp, w) 
% Single phase centroid

% Re-prepare
  global A B;
  global stdoutput qpoptim_options;
  
  dim = size(supp,1);
  n = length(stride);
  m = length(w);

% Compute initial guess
  avg_stride = ceil(mean(stride));
  c_means = supp * w' / n;
  zm = supp - repmat(c_means, [1, m]);
  c_covs = zm * diag(w) * zm' / n;

  c.supp = mvnrnd(c_means', c_covs, avg_stride)';
  %c.w = rand(1,avg_stride); c.w = c.w/sum(c.w);
  c.w = 1/avg_stride * ones(1,avg_stride);

  %load cstart.mat
  save(['cstart' num2str(n) '.mat'], 'c', 'avg_stride');
  
  X = zeros(avg_stride, m);
  D = zeros(n,1);

  iter=0;
  function  obj = d2energy(warm)
  pos=1;
  for it=1:n  
    if warm
    [D(it), X(:, pos:pos+stride(it) -1)] = ... 
    kantorovich(c.supp, c.w, supp(:,pos:pos+stride(it)-1), w(pos:pos+stride(it)-1), ...
    X(:, pos:pos+stride(it) -1));
    else
    [D(it), X(:, pos:pos+stride(it) -1)] = ... 
    kantorovich(c.supp, c.w, supp(:,pos:pos+stride(it)-1), w(pos:pos+stride(it)-1));        
    end
    pos = pos + stride(it);
  end
  obj = sum(D);
  fprintf(stdoutput, '\n\t\t %d\t %e', iter, obj );      
  end

  d2energy(false);


% optimization

  nIter = 20; 
  suppIter = 1;
  admmIter = 10;
  cterm = Inf;
  statusIter = zeros(nIter,1);
  for iter=1:nIter
      
    for xsupp=1:suppIter
    % update c.supp
    c.supp = supp * X' ./ repmat(n*c.w, [dim, 1]);

    % setup initial guess for X in ADMM
    d2energy(true);
    end
    
    % update c.w as well as X, using ADMM

    % empirical parameters
    pho = 50*mean(D);
    rho = 1.;

    % precompute linear paramters
    C = pdist2(c.supp', supp', 'sqeuclidean');
    %C=zeros(avg_stride,m);

    % lagrange multiplier
    lambda =  zeros(avg_stride, n);
    mu = 0;
    
    for admm=1:admmIter
      toc;tic;
      % step 1, update X
      pos=1;
      for i=1:n
	  vecsize = [avg_stride * stride(i), 1];
	  strips = pos:pos+stride(i)-1;
	  
	  x0 = reshape(X(:, strips), vecsize);
	  H = pho * B{avg_stride, stride(i)}; 
	  q = reshape(pho * repmat(lambda(:,i) - c.w', [1, stride(i)]) + C(:,strips), vecsize);
	  Aeq = A{avg_stride,stride(i)}(avg_stride+1:end, :);
	  beq = w(strips)';
	  [xtmp] = ... 
	  quadprog(H, q, [], [], Aeq, beq, zeros(vecsize), [], x0, qpoptim_options);
	  X(:,strips) = reshape(xtmp,[avg_stride, stride(i)]);
      
      
	  pos = pos + stride(i);
      end
      

      % step 2, update c.w
      w2 = c.w;
      %fprintf(stdoutput, '\n%f', w2);
      
      %H = n*eye(avg_stride);
      %q = - (sum(X, 2) + sum(lambda, 2));
      %[c.w] = quadprog(H, q, [], [], ones(1,avg_stride), 1, zeros(avg_stride,1), [], c.w', optim_options)';
      
      H = n * eye(avg_stride) + rho * ones(avg_stride);
      q = - (sum(X, 2) + sum(lambda, 2) + rho*(1 - mu));
      [c.w] = quadprog(H, q, [], [], [], [], zeros(avg_stride,1), [], c.w', qpoptim_options)';

      % step 3, update lambda and mu
      lambda2 = lambda; mu2 = mu;
      
      pos=1;
      for i=1:n
	  lambda(:, i) = lambda(:, i) + sum(X(:, pos:pos+stride(i)-1),2) - c.w';
      pos = pos + stride(i);
      end
      mu = mu + sum(c.w) - 1;
      
      dualres = norm(w2 - c.w);
      primres1 = norm(lambda2 - lambda, 'fro')/sqrt(n*avg_stride);
      primres2 = norm(mu2 - mu);
      %fprintf(stdoutput, '%e\t%e\t%e', primres1, primres2, dualres);
      
%       if primres1 > 10*dualres
%           pho = 2 * pho;
%           lambda = lambda/2;
%           mu = mu/2;
%           fprintf(stdoutput,' *2');
%       elseif 10*primres1 < dualres
%           pho = pho / 2;
%           lambda = lambda*2;
%           mu = mu*2;
%           fprintf(stdoutput,' /2');
%       end
      
      
    end

    % sum2one(c.w)
    c.w = c.w/sum(c.w);
    % output status
    tic;statusIter(iter) = d2energy(false);toc;
    % pause;
  end

  global statusIterRec;
  statusIterRec = statusIter;
  
  h = figure;
  plot(statusIter);
  print(h, '-dpdf', 'centroid_singlephase.pdf');
  
  fprintf(stdoutput, ' %f', c.w);
  fprintf(stdoutput, '\n');
end

