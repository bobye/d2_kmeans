function [c] = centroid_sphADMM(stride, supp, w)
% Single phase centroid using ADMM

% Re-prepare
  global A B;
  global stdoutput;

  dim = size(supp,1);
  n = length(stride);
  m = length(w);
  lenOfSegment = 10;
  numOfSegments = ceil(n/lenOfSegment);
  partitions = ceil((1:n)/lenOfSegment);

% load start guess
  avg_stride = ceil(mean(stride));
  load(['cstart' num2str(n) '.mat']);

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
    % update c.supp
    for xsupp=1:suppIter
    c.supp = supp * X' ./ repmat(n*c.w, [dim,1]);
    d2energy(true);
    end

    % update c.w using ADMM
    pho = 10*lenOfSegment*mean(D);
    C = pdist2(c.supp', supp', 'sqeuclidean');
    wseg = zeros(numOfSegments, avg_stride);
    % lagrange multiplier
    lambda =  zeros(numOfSegments, avg_stride);

    
    poss = 1;  
    for nseg = 1:numOfSegments
        idx = find(partitions == nseg);
        npar = length(idx);
        mpar = sum(stride(idx));
        Cpar = C(:,poss:poss+mpar-1);
        f = reshape(Cpar, avg_stride * mpar, 1);
        ff{nseg} = [f; zeros(avg_stride,1)];
        H{nseg} = zeros(avg_stride * (mpar+1));

        
        Aeq{nseg} = zeros(avg_stride*npar + mpar + 1, ...
                    avg_stride*(mpar +1));
        beq{nseg} = zeros(avg_stride*npar + mpar + 1, 1);

        posi=poss;pos=1;posm=1;
        for i=1:npar
            stripi = posi:posi+stride(idx(i))-1;
            strips = pos:pos+stride(idx(i))+avg_stride-1;
            stripsm= posm:posm+stride(idx(i))*avg_stride-1;
            Aeq{nseg}(strips,stripsm) = A{avg_stride, stride(idx(i))};
            beq{nseg}(strips,1) = [zeros(avg_stride,1); w(stripi)'];
            Aeq{nseg}(pos:pos+avg_stride-1,avg_stride*mpar+1:end) = -eye(avg_stride);
            
            posi= posi + stride(idx(i)) ;
            pos = pos + stride(idx(i))+avg_stride;
            posm = posm + stride(idx(i))*avg_stride;
        end
        Aeq{nseg}(pos, posm:end) = ones(1,avg_stride);
        beq{nseg}(end) = 1;

        H{nseg}(posm:end, posm:end) = pho * eye(avg_stride);

        poss = poss + mpar;    
    end
    
    for admm=1:admmIter
        toc;tic;

        % step 1
        x=cell(numOfSegments,1);
        for nseg=1:numOfSegments

            idx = find(partitions == nseg);
            npar = length(idx);
            mpar = sum(stride(idx));
            ff{nseg}(end-avg_stride+1:end) = pho * (lambda(nseg,:) - c.w);

%            if nseg > 1
		% try interior-point-convex first
                optim_options = optimset('Display','off', 'Diagnostics','off', ...
                                         'Algorithm','interior-point-convex');
                [x{nseg}, fval, exitflag] = ...
                    quadprog(H{nseg}, ff{nseg}, [], [], Aeq{nseg}, beq{nseg}, ...
                             zeros(avg_stride*(mpar+1),1), [], x{nseg}, ...
                             optim_options);
%            elseif nseg == 1

	    if exitflag < 0 % if failed, use active-set
                optim_options = optimset('Display','off', 'Diagnostics','off', ...
                                         'Algorithm','active-set');
                [x{nseg}, fval, exitflag] = ...
                    quadprog(H{nseg}, ff{nseg}, [], [], Aeq{nseg}, beq{nseg}, ...
                             zeros(avg_stride*(mpar+1),1), [], [], ...
                             optim_options);                
            end
            
            if exitflag < 0
                error('quadprog no search direction [%d %f]',exitflag, fval);
            end
            wseq(nseg,:) = x{nseg}(end-avg_stride+1:end)';
        end

        % step 2, update c.w
        w2 = c.w;       
        [c.w] = mean(wseq + lambda);

        lambda2 = lambda;
        % step 3, update lambda
        lambda = lambda + wseq - repmat(c.w, numOfSegments, 1);
        
        dualres = norm(w2 - c.w);
        primres = norm(lambda2 - lambda, 'fro')/sqrt(numOfSegments ...
                                                     * avg_stride);
        %fprintf(stdoutput, '%e\t%e', primres, dualres);
    end

    c.w = c.w/sum(c.w);
    statusIter(iter) = d2energy(false);
  end

  global statusIterRec;
  statusIterRec(:,3) = statusIter;

  h=figure;
  plot(statusIter);
  print(h, '-dpdf', 'centroid_sphADMM.pdf');

  fprintf(stdoutput, ' %f', c.w);
  fprintf(stdoutput, '\n');
end
