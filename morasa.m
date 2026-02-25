function [ img_x ] = morasa( k_data, mask, LR_rank, sparsity_threshold, tol)
% MORASA for accelerating T2 mapping
% data_k: undersampled k-space data of size Nx Ny ETL
% mask: unsersampling pattern
% LR_rank: rank of the spatiotemporal matrix, default value is 3
% sparsity_threshold: soft threshold of the sparse coefficient, default valaue is 0.01
% tol: stopping tol, a small value, default is 1e-4
% img_x: output T2-weighted image seires
% Author Xi Peng 
[Nx, Ny, ETL] = size(mask);
% POCS-based reconstruction:
% wavelet transform, need Lustig Sparse-MRI code
XFM = Wavelet('Daubechies',4,4);
appr = 1;
iter = 0;
img_x0 = ifft2c(k_data);
while (iter<100) && (appr>tol)
    iter=iter+1;
    img_x=img_x0;
    % distributed CS 
	ssx = 2^ceil(log2(Nx)); 
	ssy = 2^ceil(log2(Ny));
	ss = max(ssx, ssy);
    for i=1:ETL
        res_wave(:,:,i)=XFM*zpad(img_x(:,:,i),ss,ss);
    end
	res_wave = JointSoftThresh(res_wave, sparsity_threshold); %3D
	for i=1:ETL
        img_x(:,:,i)=crop(XFM'*res_wave(:,:,i),Nx,Ny);
    end
    % data_consistency
    img_x=reshape(img_x,Nx,Ny,ETL); 
	img_k=fft2c(img_x);
	img_k=img_k.*(~mask)+k_data;
	img_x=ifft2c(img_k);
    
	% Low rank
	img_x = reshape(img_x, Nx*Ny, ETL);
	[U S V] = svd(img_x,'econ');
	s=diag(S);
	s(LR_rank+1:ETL)=0; % strictly low rank
	img_x=U*diag(s)*V';
	img_x=reshape(img_x,Nx,Ny,ETL);
    % data_consistency
    img_x=reshape(img_x,Nx,Ny,ETL); 
	img_k=fft2c(img_x);
	img_k=img_k.*(~mask)+k_data;
	img_x=ifft2c(img_k);
    
    % Linear predictability
	for i=1:Nx
        parfor j=1:Ny
            img_x(i,j,:) = reshape(cadzow_filter(squeeze(img_x(i,j,:)), 1),[],1);
        end
    end
     
    % data_consistency
    img_x=reshape(img_x,Nx,Ny,ETL); 
	img_k=fft2c(img_x);
	img_k=img_k.*(~mask)+k_data;
	img_x=ifft2c(img_k);

	appr = norm(img_x(:)-img_x0(:))/norm(img_x0(:));
	img_x0 = img_x;
end