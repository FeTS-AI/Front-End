img_cc1=double(imread('fh20-cc-bigstudy009.png') > 128);
img_cc2=double(imread('fh20-cc-bigstudy010.png') > 128);

[dm1 lm1] = bwdist(1-img_cc1);
