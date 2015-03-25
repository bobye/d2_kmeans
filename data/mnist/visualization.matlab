f=fopen('../usps/usps.d2s_575690_c.d2', 'r');

figure;
points=40
pts_size=300
for i=1:30
  subplot(5,6,i) % first subplot
  tmp =fscanf(f, '%d\n', 1);
  tmp =fscanf(f, '%d\n', 1);
  w = fscanf(f, '%f', points);
  a = fscanf(f, '%f', points*2); a=reshape(a, [2, points]);
  scatter(a(1,:), -a(2,:), pts_size*w);
  axis equal
  set(gca, 'visible', 'off');
end

fclose(f);
