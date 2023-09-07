x = -22500:5000:22500;
y = -22500:5000:22500;
z = 0;

fid = fopen(['sites.dat'],'w');
for j = 1:length(y)
   for i = 1:length(x)
       fprintf(fid,'%G %G %G\n',x(i),y(j),z);
   end
end
status = fclose(fid);