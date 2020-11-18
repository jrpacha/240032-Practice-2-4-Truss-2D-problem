clearvars
close all

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
ndim=size(nodes,2);
  
plotLinkNodElem(nodes, elem);

%Real constants: Materials and sections area
Area=3250.0;     %section Area (in mm^2);
Y=200.0;         %Young modulus, in kN/mm^2 (1GP=1.0 kN/mm^2)
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Assembly
u=zeros(ndim*numNod,1);
Q=zeros(ndim*numNod,1);
K=zeros(ndim*numNod);

for e=1:numElem
    Ke=planeLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[ndim*elem(e,1)-1,ndim*elem(e,1),ndim*elem(e,2)-1,ndim*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Loads
%node 1:
nod=1;
Q(ndim*nod-1)=0.0;
Q(ndim*nod)=-280.0;     %kN

%node 2:
nod=2;
Q(ndim*nod-1)=0.0;
Q(ndim*nod)=-210.0;     %kN

%node 3:
nod=3;
Q(ndim*nod-1)=0.0;
Q(ndim*nod)=-280.0;     %kN

%node 4:
nod=4;
Q(ndim*nod-1)=0.0;
Q(ndim*nod)=-360.0;     %kN

%Boundary Conditions
fixedNods=[];
%node 1:
nod=1;
fixedNods=[fixedNods,ndim*nod-1];   %(u1_x=0);
fixedNods=[fixedNods,ndim*nod];     %(u1_y=0); 
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 4:
nod=4;
fixedNods=[fixedNods,ndim*nod];     %(u4_y=0); 
u(ndim*nod)=0.0;

%Reduced system
freeNods=setdiff(1:ndim*numNod,fixedNods);
Qm=Q(freeNods,1)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

%Print out displacements
displacements = [u(1:2:end),u(2:2:end)];
fprintf('\n%6s%10s%14s\n','NOD.','UX(mm)','UY(mm)')
fprintf('%4d%14.5e%14.5e\n',[(1:numNod)', displacements]')

%Post-process
%Show the original structure and the deformed one
figure()
esc=10; %scale factor to magnify displacements
plotLinkNodElemDespl(nodes,elem,u,esc)
%find out max displacements in X and Y directions and the corresponding
%nodes
[val,idx]=max(abs(displacements)); %compute the maximum value for the
                                   %components of each column and gives
                                   %that component of the column 
fprintf('\n*Max.Displ.in X direction, max|UX| =%12.5e, at Node:%2d',...
    val(1),idx(1))
fprintf('\n*Max.Displ.in Y direction, max|UY| =%12.5e, at Node:%2d\n',...
    val(2),idx(2))

%Reaction forces
R = K*u-Q;
displacements = [R(1:2:end),R(2:2:end)];
fprintf('\n%6s%10s%14s\n','NOD.','RX(kN)','RY(kN)')
fprintf('%4d%14.5e%14.5e\n',[(1:numNod)',displacements]')