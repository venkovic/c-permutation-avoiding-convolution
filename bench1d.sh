#!/bin/bash

touch bench1d.out

# Forward transforms

echo "1D forward transforms, t = 20" >> bench1d.out
stdbuf -oL ./bench1d_forward_fft bench -t 20 >> bench1d.out

echo "" >> bench1d.out  
echo "1D forward transforms, t = 22" >> bench1d.out
stdbuf -oL ./bench1d_forward_fft bench -t 22 >> bench1d.out

echo "" >> bench1d.out
echo "1D forward transforms, t = 24" >> bench1d.out
stdbuf -oL ./bench1d_forward_fft bench -t 24 >> bench1d.out

echo "" >> bench1d.out
echo "1D forward transforms, t = 26" >> bench1d.out
stdbuf -oL ./bench1d_forward_fft bench -t 26 >> bench1d.out

echo "" >> bench1d.out
echo "1D forward transforms, t = 28" >> bench1d.out
stdbuf -oL ./bench1d_forward_fft bench -t 28 >> bench1d.out

# Backward transforms

echo "" >> bench1d.out
echo "1D backward transforms, t = 20" >> bench1d.out
stdbuf -oL ./bench1d_backward_fft bench -t 20 >> bench1d.out

echo "" >> bench1d.out
echo "1D backward transforms, t = 22" >> bench1d.out
stdbuf -oL ./bench1d_backward_fft bench -t 22 >> bench1d.out

echo "" >> bench1d.out
echo "1D backward transforms, t = 24" >> bench1d.out
stdbuf -oL ./bench1d_backward_fft bench -t 24 >> bench1d.out

echo "" >> bench1d.out
echo "1D backward transforms, t = 26" >> bench1d.out
stdbuf -oL ./bench1d_backward_fft bench -t 26 >> bench1d.out

echo "" >> bench1d.out
echo "1D backward transforms, t = 28" >> bench1d.out
stdbuf -oL ./bench1d_backward_fft bench -t 28 >> bench1d.out

# Discrete convolutions

echo "" >> bench1d.out
echo "1D convolutions, t = 20" >> bench1d.out
stdbuf -oL ./bench1d_convolution bench -t 20 >> bench1d.out

echo "" >> bench1d.out
echo "1D convolutions, t = 22" >> bench1d.out
stdbuf -oL ./bench1d_convolution bench -t 22 >> bench1d.out

echo "" >> bench1d.out
echo "1D convolutions, t = 24" >> bench1d.out
stdbuf -oL ./bench1d_convolution bench -t 24 >> bench1d.out

echo "" >> bench1d.out
echo "1D convolutions, t = 26" >> bench1d.out
stdbuf -oL ./bench1d_convolution bench -t 26 >> bench1d.out

echo "" >> bench1d.out
echo "1D convolutions, t = 28" >> bench1d.out
stdbuf -oL ./bench1d_convolution bench -t 28 >> bench1d.out


