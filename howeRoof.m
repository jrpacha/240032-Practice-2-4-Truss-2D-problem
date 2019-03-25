clearvars
close all

%Exercise 2:
%Compute the displacements of the point E in the Howe roof structure 
%(see figure) when using the same bars as in the initial bridge example.

%Note: in this setting node E corresponds to node 3
numNodE=3;

%Geometry

nodes=[0.0,0.0;
       2.5e3,0.0;
       5.0e3,0.0;
       7.5e3,0.0;
       10.0e3,0.0;
       2.5e3,2.0e3;
       5.0e3,4.0e3;
       7.5e3,2.0e3];
elem=[1,2;
      2,3;
      3,4;
      4,5;
      1,6;
      2,6;
      3,6;
      6,7;
      3,7;
      7,8;
      3,8;
      4,8;
      8,5];
      
numNod=size(nodes,1);
numElem=size(elem,1);
dim=size(nodes,2);
  
plotLinkNodElem(nodes, elem);

%Real constants: Materials and sections area
Area=3250.0;     %section Area (in mm^2);
Y=200.0;         %Young modulus, in kN/mm^2 (1GP=1.0 kN/mm^2)
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Assembly
u=zeros(dim*numNod,1);
Q=zeros(dim*numNod,1);
K=zeros(dim*numNod);

for e=1:numElem
    Ke=planeLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[2*elem(e,1)-1,2*elem(e,1),2*elem(e,2)-1,2*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Loads
%node 6:
nod=6;
Q(dim*nod-1)=30.0;
Q(dim*nod)=-50.0;     %kN

%node 7:
nod=7;
Q(dim*nod-1)=0.0;
Q(dim*nod)=-40.0;     %kN

%node 8:
nod=8;
Q(dim*nod-1)=0.0;
Q(dim*nod)=-50.0;     %kN

%Boundary Conditions
fixedNods=[];
%node 1:
nod=1;
fixedNods=[fixedNods,dim*nod-1];   %(u1_x=0);
fixedNods=[fixedNods,dim*nod];     %(u1_y=0); 
u(dim*nod-1)=0.0;
u(dim*nod)=0.0;

%node 5:
nod=5;
fixedNods=[fixedNods,dim*nod];     %(u5_y=0); 
u(dim*nod)=0.0;

%Reduced system
freeNods=setdiff(1:dim*numNod,fixedNods);
Qm=Q(freeNods,1)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

%Print out displacements
fprintf('\n%6s%8s%14s\n','NOD.','UX','UY')
fprintf('%4d%14.5e%14.5e\n',[(1:numNod)',u(1:2:end),u(2:2:end)]')

%Post-process
%Show the original structure and the deformed one
figure()
esc=20; %scale factor to magnify displacements
plotLinkNodElemDespl(nodes,elem,u,esc)
%find out mac displacement and the corresponding node.
[val,idx]=max(abs(u));
dir={'X','Y'};
fprintf('\nMax.Displ., |U%s| =%12.5e, at Nod.Num:%2d\n',...
    dir{2-mod(idx,2)},val,ceil(idx/2))

fprintf('\nDisplacements of the point E (in mm):\n')
fprintf('%6s%8s%14s\n','NOD.','UX','UY')
fprintf('%4d%14.5e%14.5e\n',...
    [numNodE,u(dim*numNodE-1),u(dim*numNodE)]')       