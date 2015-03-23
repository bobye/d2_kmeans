f=fopen('mnist60k.d2s_243560_c.d2', 'r');

figure;
for i=1:30
  subplot(5,6,i) % first subplot
  tmp =fscanf(f, '%d\n', 1);
  tmp =fscanf(f, '%d\n', 1);
  w = fscanf(f, '%f', 80);
  a = fscanf(f, '%f', 80*2); a=reshape(a, [2, 80]);
  scatter(a(2,:), 27-a(1,:), 400*w);
  axis equal
  set(gca, 'visible', 'off');
end

fclose(f);
