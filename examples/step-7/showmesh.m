function showmesh(filename)

nodefile = fopen([filename '.n'],'r');
nnodes = fscanf(nodefile,'%d');
nnodes = nnodes(1);
str = fscanf(nodefile,'%20c',1);
nodedata = fscanf(nodefile,'%f',[4,nnodes]);
fclose(nodefile);

sidefile = fopen([filename '.s'],'r');
nsides = fscanf(sidefile,'%d',1);
sidedata = fscanf(sidefile,'%d',[6,nsides]);
fclose(sidefile);

X=[nodedata(2,sidedata(2,:)+1); nodedata(2,sidedata(3,:)+1)];
Y=[nodedata(3,sidedata(2,:)+1); nodedata(3,sidedata(3,:)+1)];
plot(X,Y,'k');
