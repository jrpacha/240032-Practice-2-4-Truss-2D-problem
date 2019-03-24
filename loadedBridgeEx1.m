clearvars
close all

%Exercise 1:
%Suppose that we double the section area of the elements on the bottom. 
%Compute the maximum displacement in this new design with same BC.

%Geometry
H=1800.0*sqrt(3.0);      %mm

nodes=[0.0,0.0;
       3600.0,0.0;
       7200.0,0.0;
       10800.0,0.0;
       1800.0,H;
       5400.0,H;
       9000.0,H];
elem=[1,2;
      2,3;
      3,4;
      1,5;
      2,5;
      5,6;
      2,6;
      3,6;
      6,7;
      3,7;
      4,7];

numNod=size(nodes,1);
numElem=size(elem,1);
dim=size(nodes,2);
  
plotLinkNodElem(nodes, elem);

%Real constants: Materials and sections area
Area=3250.0;     %section Area (in mm^2);
Y=200.0;         %Young modulus, in kN/mm^2 (1GP=1.0 kN/mm^2)
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Double the area of the elements at the bottom
A(1:3)=2.0*Area;

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
%node 1:
nod=1;
Q(dim*nod-1)=0.0;
Q(dim*nod)=-280.0;     %kN

%node 2:
nod=2;
Q(dim*nod-1)=0.0;
Q(dim*nod)=-210.0;     %kN

%node 3:
nod=3;
Q(dim*nod-1)=0.0;
Q(dim*nod)=-280.0;     %kN

%node 4:
nod=4;
Q(dim*nod-1)=0.0;
Q(dim*nod)=-360.0;     %kN

%Boundary Conditions
fixedNods=[];
%node 1:
nod=1;
fixedNods=[fixedNods,dim*nod-1];   %(u1_x=0);
fixedNods=[fixedNods,dim*nod];     %(u1_y=0); 
u(dim*nod-1)=0.0;
u(dim*nod)=0.0;

%node 4:
nod=4;
fixedNods=[fixedNods,dim*nod];     %(u4_y=0); 
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
esc=10; %scale factor to magnify displacements
plotLinkNodElemDespl(nodes,elem,u,esc)
%find out mac displacement and the corresponding node.
[val,idx]=max(abs(u));
dir={'X','Y'};
fprintf('\nMax.Displ., |U%s| =%12.5e, at Nod.Num:%2d\n',...
    dir{2-mod(idx,2)},val,ceil(idx/2))

















  
       