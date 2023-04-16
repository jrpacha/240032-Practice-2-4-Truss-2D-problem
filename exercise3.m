clearvars
close all

% =========================== Exercise 3 ==================================
fileSols = 'solutionEx3.txt';
fOut = fileSols;

% Data
D = 1000.0;                  %Length (mm)
F = 200.0;                   %Force (N)
Area = 200.0;                %Area (mm^2)
Y = 2000.0;                  %Young's modulus (N/mm^2)
nodDisp = 4;                 %Node disp.

fOut=fopen(fileSols,'w');
fprintf(fOut,'Length.............................. L = %e (mm)\n',      D);
fprintf(fOut,'Force............................... F = %e (N)\n',       F);
fprintf(fOut,'Area............................. Area = %e (mm^2)\n', Area);
fprintf(fOut,'Young modulus....................... Y = %e (N/mm^2)\n',  Y);
fprintf(fOut,'Num. nod. displ............... nodDisp = %d\n',     nodDisp);
fprintf(fOut,'\n');

%Goemetry
nodes=[2*D,0,0;
       -2*D,0,0;
       0,0,D;
       0,3*D,-D/2];
       
elem=[1,4;
      3,4;
      2,4];

numbering=1;
plotElementsOld(nodes,elem,numbering)
%plotElements(nodes,elem,numbering)
      
numNod=size(nodes,1);
numElem=size(elem,1);
ndim=size(nodes,2);

%Real constants
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Assembly
u=zeros(ndim*numNod,1);
Q=zeros(ndim*numNod,1);
K=zeros(ndim*numNod);

for e=1:numElem
    Ke=spatialLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[ndim*elem(e,1)-2,ndim*elem(e,1)-1,ndim*elem(e,1),...
          ndim*elem(e,2)-2,ndim*elem(e,2)-1,ndim*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Loads
%node 4:
nod=4;
Q(ndim*nod-2)=0.0;     %kN
Q(ndim*nod-1)=0.0;     %kN
Q(ndim*nod)= -F;       %kN 

%Boundary Conditions
fixedNods=[];
%node 1:
nod=1;
fixedNods=[fixedNods,ndim*nod-2];     %(u1_x=0);
fixedNods=[fixedNods,ndim*nod-1];     %(u1_y=0);
fixedNods=[fixedNods,ndim*nod];       %(u1_z=0);
u(ndim*nod-2)=0.0;
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 2:
nod=2;
fixedNods=[fixedNods,ndim*nod-2];     %(u2_x=0);
fixedNods=[fixedNods,ndim*nod-1];     %(u2_y=0);
fixedNods=[fixedNods,ndim*nod];       %(u2_z=0);
u(ndim*nod-2)=0.0;
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 3:
nod=3;
fixedNods=[fixedNods,ndim*nod-2];     %(u3_x=0);
fixedNods=[fixedNods,ndim*nod-1];     %(u3_y=0);
fixedNods=[fixedNods,ndim*nod];       %(u3_z=0);
u(ndim*nod-2)=0.0;
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%Reduced system
freeNods=setdiff(1:ndim*numNod,fixedNods);
Qm=Q(freeNods,1)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

%Reaction forces
Fr=K*u-Q;

%Plot deformed structure
esc=2
plotDeformedTruss(nodes,elem,u,esc)

%Print out displacements at node nodDisp
%fprintf(fOut,'%35s\n\n',                                      'Part (B)');
fprintf(fOut,'%38s%d (in mm)\n',             'Displacement node ',nodDisp);
fprintf(fOut,'%9s%11s%17s%17s\n',                   'NOD.','UX','UY','UZ');
fprintf(fOut,'%7d%17.5e%17.5e%17.5e\n',...
                       nodDisp,u(3*nodDisp-2),u(3*nodDisp-1),u(3*nodDisp));

%Print out all the displacements
fprintf(fOut,'\n%44s\n',                          'Displacements (in mm)');
fprintf(fOut,'%9s%11s%17s%17s\n',                   'NOD.','UX','UY','UZ');
fprintf(fOut,'%7d%17.5e%17.5e%17.5e\n',...
                          [(1:numNod)',u(1:3:end),u(2:3:end),u(3:3:end)]');

%Print out all the reaction forces                      
fprintf(fOut,'\n%45s\n',                         'Reaction forces (in N)');
fprintf(fOut,'%9s%11s%17s%17s\n',                   'NOD.','RX','RY','RZ');
fprintf(fOut,'%7d%17.5e%17.5e%17.5e\n',...
                       [(1:numNod)',Fr(1:3:end),Fr(2:3:end),Fr(3:3:end)]');
fclose(fOut);
type(fileSols);