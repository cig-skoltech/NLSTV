function D = StructureTensor2D_NL(u,W,C,bc)
% Non Local 2D Structure Tensor.
%
% u: Nx x Ny x Nc array of a vector-valued image with Nc-channels, defined
% on a Nx x Ny pixel grid
%
% W: Nx x Ny x Nw array with the Nw weights of each patch centered at
% location l=(nx,ny), 1<= nx <= Nx, 1<= ny <= Ny. These weights are inverse
% proportional to the patch distance between the two patches and are
% computed by the function NL_weights.
%
% C: Nx x Ny x Nw array with the center coordinates of the Nw patches
% which are related to the weights W 
% 
% bc: boundary condition type: 'symmetric' |'circular'|'zero'.
%
% (OUTPUT):
% D: 3D array with dimensions (Nx,Ny,3), where for each pixel (i,j)
% the last dimension of D contains the 3 unique elements of the Non local
% Structure tensor, i.e, D(i,j,:)=[fx_NL(i,j)^2,fxy_NL(i,j),fy_NL(i,j)^2] 
% where 
%                 Nw           Nc
%   fx_NL(i,j)^2= S  W(i,j,k)* S  fx(k(i),k(j),c)^2
%                k=1           c=1
%
%                Nw          Nc
%   fxy_NL(i,j)= S W(i,j,k)* S  fx(k(i),k(j),c)*fy(k(i),k(j),c)
%                k=1         c=1
%
%                 Nw           Nc
%   fy_NL(i,j)^2= S  W(i,j,k)* S  fy(k(i),k(j),c)^2
%                k=1           c=1
%
%  where k(i),k(j) indicates the coordinates of the center pixel of the 
%  patch whose distance from the patch centered at (i,j) is equal to
%  W(i,j,k).

if nargin < 4
  bc='symmetric';
end

[Nx,Ny,Nc] = size(u);
% P = Nx*Ny;

%Nw=size(W,3);

% if ~(mod(Nw, 2)) % if not all [NGx,NGy] are odd numbers
%     error('The number of weights must be odd.');
% end

grad_u  = GradOp2D(u,bc);

T1=zeros(Nx,Ny);
T2=T1;
T3=T1;
for k=1:Nc
  T1=T1+grad_u(:,:,k,1).^2;
  T2=T2+grad_u(:,:,k,1).*grad_u(:,:,k,2);
  T3=T3+grad_u(:,:,k,2).^2;
end

D = zeros(Nx,Ny,3);
D(:,:,1)=sum(W.*T1(C),3);
D(:,:,2)=sum(W.*T2(C),3);
D(:,:,3)=sum(W.*T3(C),3);


function Df=GradOp2D(f,bc)

[r,c,k]=size(f);
Df=zeros(r,c,k,2);
Df(:,:,:,1)=shift(f,[-1,0,0],bc)-f; %f(i+1,j,l)-f(i,j,l)
Df(:,:,:,2)=shift(f,[0,-1,0],bc)-f; %f(i,j+1,l)-f(i,j,l)
