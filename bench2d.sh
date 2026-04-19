#!/bin/bash

touch bench2d.out

# 2D forward transforms
echo "2D forward transforms, t = 10" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 10 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 11" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 11 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 12" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 12 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 13" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 13 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 14" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 14 >> bench2d.out

# 2D backward transforms
echo "" >> bench2d.out
echo "2D forward transforms, t = 10" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 10 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 11" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 11 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 12" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 12 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 13" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 13 >> bench2d.out

echo "" >> bench2d.out
echo "2D forward transforms, t = 14" >> bench2d.out
stdbuf -oL ./bench2d_forward_fft bench -t 14 >> bench2d.out

# 2D isotropic convolutions

echo "" >> bench2d.out
echo "2D isotropic convolution, t = 10" >> bench2d.out
stdbuf -oL ./bench2d_convolution bench -t 10 >> bench2d.out

echo "" >> bench2d.out
echo "2D isotropic convolution, t = 11" >> bench2d.out
stdbuf -oL ./bench2d_convolution bench -t 11 >> bench2d.out

echo "" >> bench2d.out
echo "2D isotropic convolution, t = 12" >> bench2d.out
stdbuf -oL ./bench2d_convolution bench -t 12 >> bench2d.out

echo "" >> bench2d.out
echo "2D isotropic convolution, t = 13" >> bench2d.out
stdbuf -oL ./bench2d_convolution bench -t 13 >> bench2d.out

echo "" >> bench2d.out
echo "2D isotropic convolution, t = 14" >> bench2d.out
stdbuf -oL ./bench2d_convolution bench -t 14 >> bench2d.out

# 2D anisotropic convolutions

echo "" >> bench2d.out
echo "2D anisotropic convolution, t = 20" >> bench2d.out
stdbuf -oL ./bench2d_convolution_anisotropic bench -t 20 >> bench2d.out

echo "" >> bench2d.out
echo "2D anisotropic convolution, t = 22" >> bench2d.out
stdbuf -oL ./bench2d_convolution_anisotropic bench -t 22 >> bench2d.out

echo "" >> bench2d.out
echo "2D anisotropic convolution, t = 24" >> bench2d.out
stdbuf -oL ./bench2d_convolution_anisotropic bench -t 24 >> bench2d.out

echo "" >> bench2d.out
echo "2D anisotropic convolution, t = 26" >> bench2d.out
stdbuf -oL ./bench2d_convolution_anisotropic bench -t 26 >> bench2d.out

