function [] = usps_blankout(ratio)

load usps_all.mat % load data

% train/test stratified split at random
idx= cell2mat(cellfun(@(x) {randperm(1100)'}, num2cell(1:10)));
train.idx = idx(1:800,:);
test.idx = idx(801:end,:);
train.data=cell2mat(cellfun(@(i) {data(:,train.idx(:,i),i)}, num2cell(1:10)));
train.label=repmat(0:9, [800,1])(:);
test.data=cell2mat(cellfun(@(i) {data(:,test.idx(:,i),i)}, num2cell(1:10)));
test.label=repmat(0:9, [300,1])(:); 


% delete feature
train.data(rand(size(train.data)) > ratio) = 0;
test.data(rand(size(test.data)) > ratio) = 0;

% save to d2s file
print_data(train, 'train', ratio);
print_data(test, 'test', ratio);
end 

function print_data(dat, tag, ratio)
filename = ['usps_blankout', num2str(ratio*100), '_', tag];
fp = fopen([filename, '.d2s'], 'w+');
for i=1:length(dat.label)
  idx=find(dat.data(:,i));
  fprintf(fp, '2\n%d\n', length(idx));
  fprintf(fp, '%d ', dat.data(idx,i));
  fprintf(fp, '\n');
  fprintf(fp, '%d ', idx);
  fprintf(fp, '\n');
end
fclose(fp);
dlmwrite([filename, '.mat'], dat.data);
dlmwrite([filename, '.label'], dat.label);
end
