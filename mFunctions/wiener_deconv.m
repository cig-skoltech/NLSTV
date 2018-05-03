function x_wiener=wiener_deconv(y,h,sigma,k)

if nargin < 4
  k=1.25;
end

if nargin < 3 || isempty(sigma)
  nc=size(y,3);
  if nc > 1
    sigma=0;
    for m=1:nc
      sigma=sigma+wmad(y(:,:,m),'db7');% estimation of the noise standard deviation.
    end
    sigma=sigma/nc;
  else
    sigma=wmad(y,'db7');% estimation of the noise standard deviation.
  end    
end


[nx,ny,nc]=size(y);

cflag=true;
if nc==1
  cflag=false;
end

[hx,hy]=size(h);

Y=fft2(y);
%zero padding h to match the image size.
hzp=h;
hzp(nx,ny)=0;

%Circular shift of the PSF so as to have consistent results by applying
%either one of the following 2 operations on an image x (hzp is the zero
%padded psf)
%1)imfilter(x,h,'conv','circular')
%2)ifft2(fft2(x).*fft2(hzp));

%============================ IMPORTANT NOTE ==============================
%Matlab's imfilter picks as the origin of h its central value, i.e
% chy=floor((hr+1)/2) and chx=floor((hc+1)/2). Therefore without circular
% shifting hzp, the above 2 operations are not consistent.
%==========================================================================

hzp=circshift(hzp,[-floor(hx/2),-floor(hy/2)]);

%x0 Initialization
H=fft2(hzp);
clear h hzp

if cflag
  H=repmat(H,[1 1 nc]);
end

% Initialization with Wiener deconvolution filter
x_wiener= real(ifft2((conj(H).*Y)./(abs(H).^2+k*sigma^2/var(y(:)))));
