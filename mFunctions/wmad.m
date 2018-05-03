function sigma=wmad(x,wfilter)
%median absolute deviation estimate of the noise

if nargin < 2
    wfilter='db7';
end
 

[Lo,Hi] = wfilters(wfilter,'d');

fsize=size(Hi.'*Hi);
if fsize(1)>size(x,1) || fsize(2)>size(x,2)
    error('wmad:: The size of the filter is larger than the input image');
end


x=perextend(x,fsize);
d=conv2(Hi,Hi,x);d=d(fsize(1):end-fsize(1)+1,fsize(2):end-fsize(2)+1);%Diagonal coefficients
d=d(1:2:end,1:2:end);


sigma=median(abs(d(:)))/.6745;


function xe=perextend(x,fsize)%Periodic extension of matrix x.
%fsize: filter size of the filter that will be convolved with matrix x.
xe=[x x(:,1:fsize(2)-1)];xe=[xe;xe(1:fsize(1)-1,:)];
