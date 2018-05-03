
mpath = mfilename('fullpath');
mpath = mpath(1:end-12);

addpath([mpath filesep 'mFunctions'])
addpath([mpath filesep 'mexFunctions' filesep 'source'])

% Grayscale Ground-truth image
f = double(imread([mpath filesep 'BSDS_images' filesep '23084.jpg']));
f = f/max(f(:));

% Initialize the seed for the random generator
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

% Blurred+Noisy image
PSF = load('levin_motion_blur.mat');
h = PSF.h; clear PSF
Hf = imfilter(f,h,'conv','circular');

BSNR = 20;
stdn = sqrt(var(Hf(:))/10^(BSNR/10));
noise = stdn*randn(size(f));
y = Hf+noise;


%%- Deblurring using the NL_STV of order p -%%

figure(1),imshow(y,[]);
fprintf('Input PSNR=%2.2f\n',psnr(y,f));

x_wiener=wiener_deconv(y,h,[],1.25);% Wiener deconvolution


% Options for the algorithm 
options_STVNL_AL={'maxiter',150,'tol',5e-5,'verbose',true,'project',...
  @(x)BoxProjection(x,[0 1]),'bc','symmetric','p',1,'showfig',false,...
  'cg_iter',2,'cg_tol',1e-5};

% Search for similar image patches
winsize = [11 11];
% Use a smoothed version of the noisy image instead of the noisy image
% itself. This leads to computing weights that lead to results of better
% reconstruction quality. 
yG=imfilter(x_wiener,fspecial('gaussian',[5 5],1.25),'conv','symmetric');
[D_STV_NL,C_STV_NL]=NL_wdist(yG,[7 7],winsize,'K',9,'bc','symmetric','isgrad',false);
W_STV_NL=exp(-D_STV_NL/(0.25)^2);

lambda = 0.0009;
[STV_NL,obj_STV_NL,res_STV_NL,ISNR_STV_NL]=deconv_STVNL_AL(y,h,lambda,...
  W_STV_NL,C_STV_NL,options_STVNL_AL{:},'alpha',10*lambda,'x_init',...
  x_wiener,'img',f);
% The ground-truth image is used only for the compuation of the increase in
% the SNR (ISNR) per iteration, but it can also be removed from the input
% arguments.

fprintf('\n Output PSNR=%2.2f\n',psnr(STV_NL,f));

figure(2),imshow(STV_NL,[]);


% Color Ground-truth image
f = double(imread([mpath filesep 'BSDS_images' filesep '97017.jpg']));
f = f/max(f(:));

% Initialize the seed for the random generator
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

% Blurred+Noisy image
h = fspecial('gaussian',[9 9],6);
Hf = imfilter(f,h,'conv','circular');

BSNR = 25;
stdn = sqrt(var(Hf(:))/10^(BSNR/10));
noise = stdn*randn(size(f));
y = Hf+noise;


%%- Deblurring using the NL_STV of order p -%%
figure(3),imshow(y,[]);
fprintf('\n Input PSNR=%2.2f\n',psnr(y,f));

x_wiener=wiener_deconv(y,h,[],1.25);% Wiener deconvolution

% Search for similar image patches

% Use a smoothed version of the noisy image instead of the noisy image
% itself. This leads to computing weights that lead to results of better
% reconstruction quality. 
yG=imfilter(x_wiener,fspecial('gaussian',[5 5],1.25),'conv','symmetric');
[D_STV_NL,C_STV_NL]=NL_wdist(yG,[7 7],winsize,'K',9,'bc','symmetric','isgrad',false);
W_STV_NL=exp(-D_STV_NL/(0.25)^2);

lambda = 0.0004;
[STV_NL,obj_STV_NL,res_STV_NL,ISNR_STV_NL]=deconv_STVNL_AL(y,h,lambda,...
  W_STV_NL,C_STV_NL,options_STVNL_AL{:},'alpha',10*lambda,'x_init',...
  x_wiener,'img',f);
% The ground-truth image is used only for the compuation of the increase in
% the SNR (ISNR) per iteration, but it can also be removed from the input
% arguments.

fprintf('\n Output PSNR=%2.2f\n',psnr(STV_NL,f));

figure(4),imshow(STV_NL,[]);
