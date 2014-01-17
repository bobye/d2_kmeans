%clear;
if ~exist('numOfSamples')
    numOfSamples = 50;
end
%%
global stdoutput ctime optim_options lpoptim_options qpoptim_options bufferc; 
stdoutput = 1;
ctime=zeros(2,1);
optim_options   = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off');
lpoptim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off', 'Simplex', 'on');
qpoptim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off', 'Algorithm','active-set');
%% Load data

fprintf(stdoutput, 'Loading data ... ');
tic;
s_modalities = 1;%2;
d_modalities = 3;%[3, 3];

%fp = fopen('../mountaindat.txt');
fp=fopen('../1500_3_5_10.txt');

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
  if (count == numOfSamples) break; end
end

fclose(fp);
toc;

%%
global max_stride;
global statusIterRec;

max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare;
clusters = d2clusters(db, 1);

n = size(statusIterRec,1);
h = figure;
%plot((1:n)', statusIterRec(:,1),'-or', ...
%     (1:n)', statusIterRec(:,2),'-dg', ...
%     (1:n)', statusIterRec(:,3),'-+b');
 
plot((1:n)', statusIterRec(:,1),'-or', ...
     (1:n)', statusIterRec(:,2),'-dg');
 
err = sqrt(kantorovich(bufferc{1}.supp, bufferc{1}.w, bufferc{2}.supp, bufferc{2}.w)) ... 
    /norm(bufferc{2}.supp,'fro');

print(h, '-dpdf', ['centroid_sphALL' num2str(numOfSamples) '.pdf']);

ctime
err

