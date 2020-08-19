% Load a pair of images
img_o=double(imread('diff_circ.png')) / 256;
img_c=double(imread('diff_circ2.png')) / 256;

% Smooth and resample to 128x128
sim_c = imresize(imfilter(img_c,fspecial('gaussian',32,4)),[128 128]);
sim_o = imresize(imfilter(img_o,fspecial('gaussian',32,4)),[128 128]);

% Create initial problem 
[vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 * size(sim_o,1) * size(sim_o,2), 1, 0.008);

% Run for a couple of iterations
tic;
[vxo2 vyo2] = lddmm_optimize_basic(vx,vy,p,1e-2,4);
toc;

%% Write deformation fields as NIFTI
defx = vxo2(:,:,6);
defy = vyo2(:,:,6);

vfout = zeros(128,128,1,1,2);
vfout(:,:,1,1,1) = defx';
vfout(:,:,1,1,2) = defy';
h = make_nii(vfout,[],[],16);
h.hdr.dime.intent_code = 1007;
save_nii(h, 'testdata/test_fft_input.nii');

%% Perform operation
resx = ifft2(fft2(defx) .* p.f_kernel_sq, 'symmetric');
resy = ifft2(fft2(defy) .* p.f_kernel_sq, 'symmetric');
vfout(:,:,1,1,1) = resx';
vfout(:,:,1,1,2) = resy';
h = make_nii(vfout,[],[],16);
h.hdr.dime.intent_code = 1007;
save_nii(h, 'testdata/test_fft_result.nii');

%% Now, save the velocity field and image for warp experiment
for i=1:size(vxo2,3)
    vfout(:,:,1,1,1) = vxo2(:,:,i)';
    vfout(:,:,1,1,2) = vyo2(:,:,i)';
    h = make_nii(vfout,[],[],64);
    h.hdr.dime.intent_code = 1007;
    save_nii(h, sprintf('testdata/test_warp_input_%02d.nii',i-1));
end

h = make_nii(sim_o', [], [], 64);
save_nii(h, 'testdata/test_warp_image.nii');

[ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(vxo2, vyo2, p);
J0 = lddmm_warp_scalar_field(sim_o, ft0x(:,:,end), ft0y(:,:,end), p);

h = make_nii(J0', [], [], 64);
save_nii(h, 'testdata/test_warp_output.nii');
        
%% Save some data for actuall LDDMM registration
h = make_nii(sim_o', [], [], 64);
save_nii(h, 'testdata/test_imreg_oshape.nii');

h = make_nii(sim_c', [], [], 64);
save_nii(h, 'testdata/test_imreg_cshape.nii');

%% Now, create a random field
[vxr vyr] = lddmm_random_field(p, 4);
for i=1:size(vxr,3)
    vfout(:,:,1,1,1) = vxr(:,:,i)';
    vfout(:,:,1,1,2) = vyr(:,:,i)';
    h = make_nii(vfout,[],[],64);
    h.hdr.dime.intent_code = 1007;
    save_nii(h, sprintf('testdata/test_randomfield_%02d.nii',i-1));
end

%%
[vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 * size(sim_o,1) * size(sim_o,2), 1, 0.008);
[vxo2 vyo2] = lddmm_optimize_basic(vx,vy,p,1e-2,4);

%% lddmm_test_gradient(vx, vy, p);
lddmm_test_gradient(vxo2, vyo2, p);
