function [x,fun_val,residual,ISNR]=prox_STVNL_AL(y,lambda,W,C,varargin)

%Proximal map of the Non Local Structure Tensor Total Variation (NL-STV) 
%regularizer under convex constraints using Augmented Lagrangian.
%
%  argmin  0.5*||y-x||^2+ lambda*||z1||_(1,p) + i_C(z2),
%  x,z1,z2
% s.t z1=D(x)
%     z2=x
%
% D is the Non Local patch-based Jacobian and D(x)D^T(x)=J(x) where J(x) is 
% the Non Local Structure Tensor.
% i_C: indicator function of the closed convex set C.
%
% ========================== INPUT PARAMETERS (required) ==================
% Parameters    Values description
% =========================================================================
% y             Vector-valued noisy image of size Nx x Ny x Nc, 
%               where Nc: number of channels.
% lambda        Regularization penalty parameter.
% W             Nx x Ny x Nw array which is given by function NL_weights 
%               and contains the Nw weights of each patch centered at 
%               location n=(nx,ny), 1<= nx <= Nx, 1<= ny <= Ny. 
%               These weights are related to the patch distance between the 
%               patch in n=(nx,ny) and the ones centered in locations 
%               m_k=(mx_k,my_k) stored in array the C. C is provided as an 
%               output of the function NL_weights.
% C             Nx x Ny x Nw array which for each pixel location n=(nx,ny)
%               contains the centers of the Nw closest patches to the one
%               centered in n. 
% ======================== OPTIONAL INPUT PARAMETERS ======================
% Parameters    Values' description
%
% img           Original Image. (For the compution of the ISNR improvement)
% x_init        Initial solution. (Default: y)
% maxiter       Number of iterations (Default: 100)
% tol           Stopping threshold for denoising (Default:5e-5)
% verbose       If verbose is set to true then info for each iteration is
%               printed on screen. (Default: false)
% showfig       If showfig is set to true the result of the reconstruction 
%               in every iteration is shown on screen. (Default: false)
% project       A function handle for the projection onto the convex set C.
%               (Default: project=@(x)x, which means that there is no
%               constraint on x.)
% bc            Boundary conditions for the differential operators.
%               {'symmetric'|'circular'|'zero'} (Default: 'symmetric')
% p             Specifies the order of the Schatten norm with p >= 1.
%               (Default: p=1).               
% alpha         Augmented Lagrangian penalty parameter (Default: 1e-2).
% cg_iter       Congugate gradient internal iterations (Default: 10).
% cg_tol        stopping threshold for CG (Default 1e-4).
% =========================================================================
% ========================== OUTPUT PARAMETERS ============================
% x             Denoised image.
% fun_val       Evolution of the objective function over the iterations.
% residual      Evolution of the residual over the iterations. 
% ISNR          Evolution of the SNR improvement over the iterations.
% =========================================================================
%
% Author: s.lefkimmiatis@skoltech.ru
%
% =========================================================================


[maxiter,alpha,tol,verbose,img,project,bc,p,showfig,cg_iter,...
  cg_tol,x_init]=process_options(varargin,'maxiter',100,'alpha',1e-2,...
  'tol',5e-5,'verbose',false,'img',[],'project',[],'bc','symmetric',...
  'p',1,'showfig',false,'cg_iter',10,'cg_tol',1e-4,'x_init',[]);

pflag=true;
if isempty(project)
  pflag=false;
  project=@(x)x;
end

% snorm = {'spectral'|'nuclear'|'frobenius'|'lp'} %snorm=Schatten norms
% If snorm is set to 'lp', order also has to be specified ( 1<=order<=inf)

%bc: Boundary conditions for the differential operators
%    ('symmetric'|'circular'|'zero')

if p < 1
  error('prox_STVNL_AL: Invalid Schatten norm order.');
end

Nw=size(W,3);

[idx,n,I]=compute_Nsum_indices(C);

if isempty(x_init)
  x=y;
else
  x=x_init;
end

count=0;
W=sqrt(W);

fun_val=zeros(maxiter,1);
residual=fun_val;

if isempty(img)
  ISNR=0;
else
  ISNR=fun_val;
end

[nx, ny, nc]=size(y);
s1=zeros([nx ny 2 Nw*nc]);

if pflag
  s2=zeros(nx,ny,nc);
end
  
y_over_alpha=y/alpha;

if verbose
  fprintf('*********************************************************\n');
  fprintf('**         Denoising with NL-STV Regularizer           **\n');
  fprintf('*********************************************************\n');
  fprintf('#iter     relative-dif    fun_val          residual         ISNR\n')
  fprintf('====================================================================\n');
end

for i=1:maxiter
  
 %u1=D(x)+s1
 %u2=x+s2;
 
 z1=proxSpMat2xNc(JacobianOp2D_NL(x,W,C,bc)+s1,p,lambda/alpha);
 
 if pflag
   z2=project(x+s2);
   
   %v1=z1-s1
   %v2=z2-s2
   b=y_over_alpha+AdjJacobianOp2D_NL(z1-s1,W,idx,n,I,bc)+(z2-s2);
   k=(alpha+1)/alpha;
   A=@(x)JTJOp2D_NL(x,W,idx,n,I,bc)+k*x;
   %xnew=ReconCG(b, A, x, 10);
   xnew=CG_solver(A,b,x,cg_iter,cg_tol);
   tmp1=JacobianOp2D_NL(xnew,W,C,bc)-z1;
   tmp2=xnew-z2;
   s1=s1+tmp1;
   s2=s2+tmp2;
 else
   %v1=z1-s1
   b=y_over_alpha+AdjJacobianOp2D_NL(z1-s1,W,idx,n,I,bc);
   A=@(x)JTJOp2D_NL(x,W,idx,n,I,bc)+(x/alpha);
   xnew=CG_solver(A,b,x,cg_iter,cg_tol);
   tmp1=JacobianOp2D_NL(xnew,W,C,bc)-z1;
   s1=s1+tmp1;
 end
 

  re=norm(xnew(:)-x(:))/norm(xnew(:));%relative error
  if (re<tol)
    count=count+1;
  else
    count=0;
  end
 
  
  x=xnew;
  
  if verbose
    if ~isempty(img)
      if pflag
        residual(i)=sqrt(norm(tmp1(:),2)^2+norm(tmp2(:),2)^2);
      else
        residual(i)=norm(tmp1(:),2);
      end
      fun_val(i)=cost(y,x,W,C,lambda,bc,p);
      ISNR(i)=20*log10(norm(y(:)-img(:))/norm(x(:)-img(:))); %#ok<AGROW>
      % printing the information of the current iteration
      fprintf('%3d \t %10.5f \t %10.5f \t %10.5f \t %2.8f\n',i,re,fun_val(i),residual(i),ISNR(i));
    else
      if pflag
        residual(i)=sqrt(norm(tmp1(:),2)^2+norm(tmp2(:),2)^2);
      else
        residual(i)=norm(tmp1(:),2);
      end
      fun_val(i)=cost(y,x,W,C,lambda,bc,p);
      % printing the information of the current iteration
      fprintf('%3d \t %10.5f \t %10.5f \t %10.5f\n',i,re,fun_val(i),residual(i));
    end
  end
  
%   if i==1
%     fprintf('Recomputing Weights....\n');
%     fG=imfilter(x,fspecial('gaussian',[9 9],1.25),'conv','symmetric');
%     [W,C]=NL_weights(fG,[9 9],[21 21],9,0.25,bc);
%     [idx,n,I]=compute_Nsum_indices(C);
%   end
  
  if showfig
    fh=figure(1);
    figure(fh);
    if verbose && ~isempty(img)
      msg=['iteration: ' num2str(i) ' ,ISNR: ' num2str(ISNR(i))];
      set(fh,'name',msg);imshow(x,[]);
    else
      imshow(x,[]);
    end
  end
  
  if count >=5
    ISNR(i+1:end)=[];
    fun_val(i+1:end)=[];
    residual(i+1:end)=[];
    break;
  end
end



function [Q,NLSTVnorm]=cost(y,f,W,C,lambda,bc,p)

D = StructureTensor2D_NL(f,W.^2,C,bc);
D_trace=D(:,:,1)+D(:,:,3);%trace of D
D_det=sqrt((D(:,:,1)-D(:,:,3)).^2+4*D(:,:,2).^2);%sqrt of determinant of D

switch p
  case inf
    NLSTVnorm=max(sqrt((D_trace(:)+D_det(:))/2),sqrt(abs(D_trace(:)-D_det(:))/2));
    NLSTVnorm=sum(NLSTVnorm(:));
  case 2
    NLSTVnorm=sqrt(D_trace(:));
    NLSTVnorm=sum(NLSTVnorm(:));
  otherwise
    NLSTVnorm=(((D_trace(:)+D_det(:))/2).^(p/2)+(abs(D_trace(:)-D_det(:))/2).^(p/2)).^(1/p);
    NLSTVnorm=sum(NLSTVnorm(:));
end


Q=0.5*norm(y(:)-f(:),2)^2+lambda*NLSTVnorm;


function [x,delta,i]=CG_solver(A,b,x,iter,tol)

% Implements the conjugate gradient method to solve the linear system Ax=b, 
% where A is a function handle for the operation Ax.
%
% Inputs parameters are:
%    b    : particular solution b
%    x    : starting value for x
%    iter : maximum CG iterations
%    tol  : threshold for stoping CG
% Output parameters are:
%    xnew  : solution of conjudate gradient
%    delta : squared norm of the residual r=A(x)-b.
%    i     : number of iterations used


i=0;
r=b-A(x);
d=r;
delta=r(:)'*r(:);

while(i < iter && delta > tol)
  % q=A*d; 
  
  q=A(d);
  alfa=delta/(d(:)'*q(:));
  
  x=x+alfa*d;
  
  if ~mod(i,50)
    r=b-A(x);
  else
    r=r-alfa*q;
  end
  
  delta_old=delta;
  delta=r(:)'*r(:);
  beta=delta/delta_old;
  d=r+beta*d;
  i=i+1;
end
