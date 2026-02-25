function [ img_comb ] =morasa_comb( k_data, mask, LR_rank, sparsity_threshold, tol)
% MORASA for accelerating T2 mapping
% data_k: undersampled k-space data of size Nx Ny ETL
% mask: unsersampling pattern
% LR_rank: rank of the spatiotemporal matrix, default value is 3
% sparsity_threshold: soft threshold of the sparse coefficient, default valaue is 0.01
% tol: stopping tol, a small value, default is 1e-4
% img_x: output T2-weighted image seires
% Author Xi Peng 
[Nx, Ny, coil,ETL] = size(k_data);
% [Nx, Ny, ETL] = size(mask);
% POCS-based reconstruction:
% wavelet transform, need Lustig Sparse-MRI code
XFM = Wavelet('Daubechies',4,4);
appr = 1;
iter = 0;
for i=1:ETL
	img_x0(:,:,:,i)=ifft2c(k_data(:,:,:,i));
end

% estimate sensitivity map from the first time point
sensi_data = k_data(:,:,:,1);
sensi_img = ifft2c(sensi_data);
csm = sensi_img./repmat(sqrt(sum(sensi_img.*conj(sensi_img),3)),[1,1,coil]);

% coil combine
    img_comb0= squeeze(sum(bsxfun(@times, img_x0, conj(csm)), 3));
    iterations=50;
    
while (iter<iterations) && (appr>tol)
    iter=iter+1;
    img_comb=img_comb0;
    % distributed CS 
	ssx = 2^ceil(log2(Nx)); 
	ssy = 2^ceil(log2(Ny));
	ss = max(ssx, ssy);
    for i=1:ETL
        res_wave(:,:,i)=XFM*zpad(img_comb(:,:,i),ss,ss);
    end
	res_wave = JointSoftThresh(res_wave, sparsity_threshold); %3D
	for i=1:ETL
        img_comb(:,:,i)=crop(XFM'*res_wave(:,:,i),Nx,Ny);
    end
    
    %%data consistency 1
        for i=1:coil
            img_x(:,:,i,:)=img_comb.*repmat(squeeze(csm(:,:,i)),[1,1,ETL]);
        end
        img_x=reshape(img_x,Nx,Ny,coil,ETL); 
        for i=1:ETL
            img_k(:,:,:,i)=fft2c(img_x(:,:,:,i));
        end
        for i=1:coil
            img_k(:,:,i,:)=squeeze(img_k(:,:,i,:)).*(~mask)+squeeze(k_data(:,:,i,:));
        end
        for i=1:ETL
            img_x(:,:,:,i)=ifft2c(img_k(:,:,:,i));
        end   
    img_comb= squeeze(sum(bsxfun(@times, img_x, conj(csm)), 3));
    
	% Low rank
	img_comb = reshape(img_comb, Nx*Ny, ETL);
	[U S V] = svd(img_comb,'econ');
	s=diag(S);
	s(LR_rank+1:ETL)=0; % strictly low rank
	img_comb=U*diag(s)*V';
	img_comb=reshape(img_comb,Nx,Ny,ETL);
    
    % data_consistency
        for i=1:coil
            img_x(:,:,i,:)=img_comb.*repmat(csm(:,:,i),[1,1,ETL]);
        end
        img_x=reshape(img_x,Nx,Ny,coil,ETL); 

        for i=1:ETL
            img_k(:,:,:,i)=fft2c(img_x(:,:,:,i));
        end
        for i=1:coil
            img_k(:,:,i,:)=squeeze(img_k(:,:,i,:)).*(~mask)+squeeze(k_data(:,:,i,:));
        end
        for i=1:ETL
            img_x(:,:,:,i)=ifft2c(img_k(:,:,:,i));
        end
     img_comb= squeeze(sum(bsxfun(@times, img_x, conj(csm)), 3));
    % Linear predictability
	for i=1:Nx
        parfor j=1:Ny
            img_comb(i,j,:) = reshape(cadzow_filter2(squeeze(img_comb(i,j,:)), 0.02),[],1);
        end
        
    end
    
    %%data consistency 
        for i=1:coil
            img_x(:,:,i,:)=img_comb.*repmat(csm(:,:,i),[1,1,ETL]);
        end
        img_x=reshape(img_x,Nx,Ny,coil,ETL); 
     
        for i=1:ETL
            img_k(:,:,:,i)=fft2c(img_x(:,:,:,i));
        end
        for i=1:coil
            img_k(:,:,i,:)=squeeze(img_k(:,:,i,:)).*(~mask)+squeeze(k_data(:,:,i,:));
        end
        for i=1:ETL
            img_x(:,:,:,i)=ifft2c(img_k(:,:,:,i));
        end
        img_comb= squeeze(sum(bsxfun(@times, img_x, conj(csm)), 3));

        appr=norm(img_comb(:)-img_comb0(:))/norm(img_comb0(:));
        img_comb0 = img_comb;
      
      sprintf('appr=%04f  for iter=%d',appr,iter) 
end
end
