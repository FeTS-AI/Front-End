function [dx1, dy1] = lddmm_register_multiscale(Is, Im, rlist, scale_list)


if length(scale_list) == 1
    scale_list = ones(length(rlist), 1) * scale_list;
end;

for k = 1:length(rlist);
    r = rlist(k);
    nb_iter_per_scale = scale_list(k);
    
        fprintf(2, '---> scale=%f, iter=%d\n', r, nb_iter_per_scale);
    
    Is1 = imresize(Is, r, 'bilinear');
    Im1 = imresize(Im, r, 'bilinear');

    if k == 1
        [dx1, dy1, vx1, vy1] = lddmm_basic(Is1, Im1, nb_iter_per_scale);
    else
        vx1 = zeros(size(Is1,1), size(Is1,2), size(vx0, 3));
        vy1 = zeros(size(vx1));
        for ch = 1:size(vx0, 3)
            vx1(:, :, ch) = imresize(vx0(:,:,ch) * r / rlist(k-1), size(Is1), 'bilinear');
            vy1(:, :, ch) = imresize(vy0(:,:,ch) * r / rlist(k-1), size(Is1), 'bilinear');
        end;
        [dx1, dy1, vx1, vy1] = lddmm_basic(Is1, Im1, nb_iter_per_scale, vx1, vy1);
    end;
    



    
    vx0 = vx1;
    vy0 = vy1;
    
%    pause;
end;