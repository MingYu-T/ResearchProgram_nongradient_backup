function obj = objfunc(var)
% - Global Variables
    global Ncp N E0 nu nele thickness...
    gobalcoords0 edofMat fixeddof loaddof F ndof Kin Kout outdof volfrac varfail splineR 

% - check if FEM fail
    varfail = var; 

% - Arrange Variables
    sv_num = (Ncp-2)*2+1;
    Bx = zeros(N,Ncp);
    By = zeros(N,Ncp);
    Br = zeros(N,1);
    for i = 1:N
        Bx(i,2:end-1) = var( 1+(i-1)*sv_num:1+(i-1)*sv_num+(Ncp-3) );
        By(i,2:end-1) = var( (1+Ncp-2)+(i-1)*sv_num:(1+Ncp-2)+(i-1)*sv_num+(Ncp-3) );
        Br(i)=var(sv_num*i);
        if Br(i)<0.1e-3
            Br(i) = 0;
        end
    end
    
    % Spline Constrained Ends
    Bx(1,1) = 0;By(1,1) = 0;Bx(1,end) = 10e-3;By(1,end) = 5e-3;
    Bx(2,1) = 0;By(2,1) = 0;Bx(2,end) = 10e-3;By(2,end) = 5e-3;

    Bx(3,1) = 0;By(3,1) = 0;Bx(3,end) = 0;By(3,end) = 5e-3;
    Bx(4,1) = 0;By(4,1) = 0;Bx(4,end) = 0;By(4,end) = 5e-3;

    Bx(5,1) = 0;By(5,1) = 5e-3;Bx(5,end) = 10e-3;By(5,end) = 5e-3;
    Bx(6,1) = 0;By(6,1) = 5e-3;Bx(6,end) = 10e-3;By(6,end) = 5e-3;
 

% - Topology Description Function 
    Phi_s = generateTDF(Bx,By,Br);
    
% - Heaviside Function
    H = Heaviside(Phi_s);

% - volume penalty
    volnum = find(Phi_s>0);
    vol = numel(volnum)/numel(Phi_s);
    vol_p = max(0,vol-volfrac);

% - Commands generation
    generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,edofMat,fixeddof,loaddof,Kin,Kout,outdof,H);

% - Run ANSYS APDL as a Subroutine
    [sta, cmd] = dos('apdl.bat','-echo');

% - Read Nodal Displacement Solution
    U=load('NODALDISPLACEMENT.txt');
    delete('NODALDISPLACEMENT.txt');  
    obj=U(1)+vol_p;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heaviside Function
function H=Heaviside(phi)
    [nodey,nodex] = size(phi);
    num1=find(phi>0);
    H(num1)=1;
    num2=find(phi<=0);
    H(num2)=1e-3;
    H = reshape(H,nodey,nodex);
end

function  generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,edofMat,fixeddof,loaddof,Kin,Kout,outdof,H)
global gen
fid = fopen('command.txt','w');
ndot=ndof/2;
fprintf(fid,'/PREP7\n');
%define node
for i=1:ndot
    fprintf(fid,'N,%d,%G,%G,%G\n',i,gobalcoords0(2*i-1),gobalcoords0(2*i),0);
end
%define element type and the thickness
fprintf(fid,'et,1,plane182\nKEYOPT,1,1,0\nKEYOPT,1,3,3\nkeyopt,1,6,0\nR,1,%G, \nMPTEMP,1,0 \n',thickness);
edotMat=0.5*edofMat(:,[2,4,6,8]);
%calculate the element parameter
[nodey,nodex] = size(H);
nelx = nodex-1;
nely = nodey-1;
temp = zeros(nely,nelx); 
for i = 1:nely
    for j = 1:nelx
        temp(i,j) = ( H(i,j)+H(i+1,j)+H(i,j+1)+H(i+1,j+1) )/4;  % ersatz material model
    end
end
xba = flip(temp);
xba = xba(:);
if gen>10
    xt = 0.5;
else
    xt = 0.25;
end
solid = find(xba>=xt);
xba(solid) = 1;
% xba = reshape(xba,nely,nelx);
MU0=E0*xba/(2*(1+nu));
K0=E0*xba/(2*(1-2*nu));
d=2./K0;
snele=0;
for i=1:nele
    if xba(i)>=0.5
        snele=snele+1;
        %define original elements
        fprintf(fid,'TB,HYPE,%d,1,2,NEO\nTBTEMP,0\n',snele);
        fprintf(fid,'TBDATA,,%G,%G,,,,\n',MU0(i),d(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',snele);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
    end
end
%%%%%define node of the spring
coords0w=reshape(gobalcoords0,2,[]);
lx=max(coords0w(1,:))-min(coords0w(1,:));
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+1,gobalcoords0(loaddof)+2*lx,gobalcoords0(loaddof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+1);
fprintf(fid,'d,%d,uy,0\n',ndot+1);
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+2,gobalcoords0(outdof)+2*lx,gobalcoords0(outdof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+2);
fprintf(fid,'d,%d,uy,0\n',ndot+2);
%%%%%define the spring
fprintf(fid,'ET,2,LINK180\nKEYOPT,2,2,0\nR,2,1, ,0\n');
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+1,Kin*2*lx,2*nele+1,0.3,2*nele+1,(loaddof+1)/2,ndot+1);
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+2,Kout*2*lx,2*nele+2,0.3,2*nele+2,(outdof+1)/2,ndot+2);
%apply the displacement
nfix=size(fixeddof,1);
for i=1:nfix
    if mod(fixeddof(i),2)==1
        fprintf(fid,'d,%d,ux,0\n',(fixeddof(i)+1)/2);
    else       
        fprintf(fid,'d,%d,uy,0\n',fixeddof(i)/2);
    end
end
fprintf(fid,'d,%d,uz,0\n',(loaddof+1)/2);
fprintf(fid,'d,%d,uz,0\n',(outdof+1)/2);
%apply the external load
nload=size(loaddof,1);
for i=1:nload
    if mod(loaddof(i),2)==1
    fprintf(fid,'F,%d,fx,%G\n',(loaddof(i)+1)/2,full(F(loaddof(i))));
    else
    fprintf(fid,'F,%d,fy,%G\n',loaddof(i)/2,full(F(loaddof(i))));
    end
end
%solve 
fprintf(fid,'finish\n/sol\nANTYPE,0\nNLGEOM,1\nNSUBST,1,0,0\n');
fprintf(fid,'CNVTOL,U,-1, \n'); 
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',sum(abs(full(F))));
fprintf(fid,'OUTRES,ERASE\nOUTRES,ALL,ALL\n/status,solu\nsolve\n');
%post processing---the nodal displacement 
fprintf(fid,'/POST1\n');
fprintf(fid,'SET, , ,1, ,1, , \n');
fprintf(fid,'*cfopen,NODALDISPLACEMENT,txt\n');
fprintf(fid,'*get,ndtemp1,node,%d,ux\n',(outdof+1)/2);
fprintf(fid,'*vwrite,ndtemp1\n');
fprintf(fid,'(E13.5)\n');
fprintf(fid,'*cfclose\n');
fprintf(fid,'finish\n');
fclose(fid);
end