clear;
%%
global stdoutput optim_options; 
stdoutput = 1;
optim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off');

%% Load data

fprintf(stdoutput, 'Loading data ... ');
tic;
s_modalities = 2;
d_modalities = [3, 3];


fp = fopen('../mountaindat.txt');

for i=1:s_modalities
    db{i}.stride = [];
    db{i}.w = [];
    db{i}.supp = [];
end

count = 0;
while ~feof(fp)
  count = count +1;
  for i=1:s_modalities      
    fscanf(fp, '%d', 1);
    [d check] = fscanf(fp, '%d', 1); 
    if check == 0 break; end

    db{i}.stride(end+1) = d;
    we = fscanf(fp, '%f', [1, d]);
    db{i}.w(1,(end+1):(end+d)) = we/sum(we);
    db{i}.supp(:, (end+1):(end+d)) = fscanf(fp, '%f', [d_modalities(i), d]);      
  end
  if (count == 200) break; end
end

fclose(fp);
toc;

%%
global max_stride;
max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare;
clusters = d2clusters(db, 1);
