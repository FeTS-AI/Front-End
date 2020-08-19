function depth = depthmap(dm, lm)

i_inside = find(dm > 0);
lm_inside = lm(i_inside);
lm_unique = unique(lm_inside);

depth = zeros(size(dm));

for i = 1:length(lm_unique())
    v = lm_unique(i);
    iv = i_inside(find(lm_inside == v));
    div = dm(iv);
    dmax = max(div);
    depth(iv) = dm(iv) / dmax;
end
