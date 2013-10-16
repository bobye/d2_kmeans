function [clusters] = d2clusters( db, k)
  %%  
  % INPUT
  % k: number of clusters
  % db.stride: size of supports
  % db.w: prob of supports
  % db.supp: supports

  % OUTPUT
  % clusters: k convergent d2 clusters
  %

  global stdoutput;
  
  if nargin == 1
    k = 5;
  end

  clusters = [];
  n = length(db{1}.stride); % size of total samples
  labels = randi(k,1,n);

  nROUND = 1;

  for i=1:nROUND
    fprintf(stdoutput, 'Round %d ... ', i);
    % relabel based on distance

    % compute the centroids
    for j=1:k
      fprintf(stdoutput, '\n\t cluster %d - ', j);
      clusters{j} = centroid(j, labels, db);
    end

  end
end



function [c] = centroid( lb, labels, db)
% INPUT
% lb: to compute the multiphase centroid of points with label = lb

% OUTPUT
% c: multiphase centroid of selected points
  global stdoutput;   
  nphase = length(db);

  for i=1:nphase
    fprintf(stdoutput, '\n\t\tphase %d: ', i);
    warmlabels{i} = getwarm(labels, db{i}.stride);
    w{i} = db{i}.w(lb == warmlabels{i});
    supp{i} = db{i}.supp(:,lb==warmlabels{i});    
    stride{i} = db{i}.stride(lb == labels);
    c{i} = centroid_singlephase(stride{i}, supp{i}, w{i});
    c{i} = centroid_sphLP(stride{i}, supp{i}, w{i});
    pause;return;
  end
  
end


function [warmlabels] = getwarm(lbs, stride)
% INPUT
% lbs: [1,2,1,3,2]
% stride: [4,3,2,3,4]

% OUTPUT
% warmlabels: [1 1 1 1 2 2 2 1 1 3 3 3 2 2 2 2]
%

  len = sum(stride);
  warmlabels = zeros(1,len);
  pos=1;
  for j=1:length(stride)
    warmlabels(pos:pos+stride(j)-1) = lbs(j);
    pos = pos + stride(j);
  end

end



