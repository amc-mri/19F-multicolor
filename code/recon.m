% Example code for "A novel iterative sparse deconvolution method for multicolor 19F-MRI"
% This code runs a multicolor 19F reconstruction
% Provided data contains four direction readout 19F FLASH scans + 1H proton FLASH
% Images show a mouse injected with PFCE and PFOB in the leg muscles
% For more details, please refer to paper 

% Jasper Schoormans - 2019
% This code relies on the SPOT linear operator toolbox (v1.2, http://www.cs.ubc.ca/labs/scl/spot/)

clear all; close all; clc; 
addpath('functions')
addpath(genpath('spot-master'))

%% reconstruction parameters
dim= 3;                     % shared phase-encoding dimension
phaseremoval=1;             % option to remove image phase

% CG OPTIONS
lambda=2e1;                 % regularization parameter
niter=25;                     % number of iterations for conjugate-gradient algorithm

%% Loading data
fprintf('--------------------------------------------------------------------\n')
fprintf('A novel iterative sparse deconvolution method for multicolor 19F-MRI\n')
fprintf('--------------------------------------------------------------------\n')
fprintf('Reconstruction of a multicolor 19F scan\n')
fprintf('Mouse injected with PFOB and PFCE\n\n')

fprintf('Loading k-space and image data\n')
load('kspaces.mat')
load('protonimage.mat')

%% shift correction - registers k-spaces 1,2,3 to k-space 4
% needs image processing toolbox 
% if not available, load the pre-registered images instead:
% load('kspaces_registered.mat')

fprintf('Calculating shift correction of direction 1...\n')
K1 = registration_correction_PFCE(K1,K4);
fprintf('Calculating shift correction of direction 2...\n')
K2 = registration_correction_PFCE(K2,K4);
fprintf('Calculating shift correction of direction 3...\n')
K3 = registration_correction_PFCE(K3,K4);

%% loop over slices

for sl=85
fprintf('Starting Reconstruction of slice %i\n',sl)
    
% define functions, operators and derive parameters from data
fprintf('define functions, operators and derive parameters from data \n')
nx=size(K1,1);      % imsize 1
ny=size(K1,2);      % imsize 2
N=nx*ny;            % number of pixels in one image

% defining local functions (reshaping ops vector <-> tensor)
rr = @(I) reshape(I,[nx,ny*4]);
vec= @(I) reshape(I,[numel(I), 1]);
rr1 = @(I) reshape(I,[nx,ny]);
rr2 = @(I) reshape(I,[nx,ny*2]);

normalize = @(I) I./max(abs(I(:)));
ifft_meas = @(I,dim) fftshift(ifft(ifftshift(I,dim),[],dim),dim);       % iFFT in measurement direction

% defining matrix operators 
F=opDFT2(nx,ny,1);                                      % Fourier operator
ShiftOp=opConvolve(nx,ny,1,[nx/2 ny/2],'cyclic');       % Fourier shift operator 
FS=F*ShiftOp;                                           % combined
F2=opBlockDiag(FS,FS,FS,FS);                            % 4 Fourier ops - one for each image 

% preprocess data to reconstruct (there is only one shared phase-encoding direction...)
fprintf('normalizes, ifft in shared phase-encoding direction... \n')

K1=normalize(K1);
K2=normalize(K2);
K3=normalize(K3);
K4=normalize(K4);

k1=extract_slice(ifft_meas(K1,dim),sl,3);
k2=extract_slice(ifft_meas(K2,dim),sl,3);
k3=extract_slice(ifft_meas(K3,dim),sl,3);
k4=extract_slice(ifft_meas(K4,dim),sl,3);

data=[vec(k1);vec(k2); vec(k3); vec(k4)]; % data matrix 

if phaseremoval
fprintf('removing phase of image data... \n')
data= F2*abs(opInverse(F2)*data); % 
end

k1=data(1:N);
k2=data(N+1:N*2);
k3=data(2*N+1:3*N);
k4=data(3*N+1:N*4);

% Linear recon
fprintf('linear (zero-filled) recon... \n')
linear_recon=opInverse(F2)*data;    %linear recon of data 

% Visualization of linear recon
figure(1);subplot(211); imshow(abs(rr(linear_recon)),[]); title('magnitude linear recon') 
figure(1);subplot(212); imshow(angle(rr(linear_recon)),[-pi, pi]); title('phase linear recon')


% Calculating the spectrum
fprintf('Calculating the spectrum... \n')
BW=4.4643e+04;        
ppm=282.5685;              
BWpix=BW/nx;
[PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(ppm,BWpix);
PFOB_alpha=[251 735 354 625 355];                %manually changing peak heights 
pixlocs=1+round(-PFOB); 
pixlocs(pixlocs<0)=nx+pixlocs(pixlocs<0);
Spectrum_BW=zeros(1,nx); Spectrum_BW(pixlocs)=PFOB_alpha./sum(PFOB_alpha(:));

% Visualizing the deconvolution spectrum
figure(2);
plot(abs(Spectrum_BW)); hold on; plot(abs(Spectrum_BW)); hold off


% make convolution operators (include fftshift in operators for simplicity)
Spectrum_flipped=flip(Spectrum_BW,2);
Spectrum_vert=Spectrum_BW.';
Spectrum_vert_flipped=flip(Spectrum_BW.',1);

A1=opConvolve(nx,ny,Spectrum_vert,[1 1],'cyclic');
A2=opConvolve(nx,ny,Spectrum_vert_flipped,[1 1],'cyclic');
A3=opConvolve(nx,ny,Spectrum_BW,[1 1],'cyclic');
A4=opConvolve(nx,ny,Spectrum_flipped,[1 1],'cyclic');
B=opDirac(nx*ny);  

% constructing measurement operator
M=[FS*A1,FS*B;...
    FS*A2,FS*B;...
    FS*A3,FS*B;...
    FS*A4,FS*B] 

% first guess (pseudo-inverse) solution
first_guess=pinv(M)*data;
figure(20); imshow(rr2(abs(first_guess)),[]); axis off; title('first guess'); colormap('jet')

% LASSO (conjugate gradient) solution
RCG=nl_conjgrad_fluor_test(M,data,first_guess,niter,zeros([2*N,1]),lambda,nx,ny*2,1);

% visualizing result
imPFCE(:,:,sl)=rr1(abs(RCG(N+1:2*N)));
imPFOB(:,:,sl)=rr1(abs(RCG(1:N)));

figure(3); 
subplot(121)
imshow(imPFOB(:,:,sl),[]); axis off;
title('PFOB IMAGE')
subplot(122)
imshow(imPFCE(:,:,sl),[]); axis off;
title('PFCE IMAGE')

end
%% overlay on proton image

params = CreateOverlayImageParams()
params.export=0
params.fignumber=4
params.threshold=0.1
params.Alpha=0.5;
params.n=2
params.slicerange=85
params.maxfactor=1;
params.minfactor=1;
params.orientation='cor'

CreateOverlayImage(PROTONIM,imPFCE,imPFOB,params)


 
 
 
 
 
 
 
 
 


