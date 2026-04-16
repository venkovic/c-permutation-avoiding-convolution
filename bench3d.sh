#!/bin/bash

touch bench3d.out

# 3D isotropic convolutions

echo "" >> bench3d.out
echo "3D isotropic convolution, t = 7" >> bench3d.out
./bench3d_convolution bench -t 7 >> bench3d.out

echo "" >> bench3d.out
echo "3D isotropic convolution, t = 8" >> bench3d.out
./bench3d_convolution bench -t 8 >> bench3d.out

echo "" >> bench3d.out
echo "3D isotropic convolution, t = 9" >> bench3d.out
./bench3d_convolution bench -t 9 >> bench3d.out

# 3D anisotropic convolutions

echo "" >> bench3d.out
echo "3D anisotropic convolution, t = 20" >> bench3d.out
./bench3d_convolution_anisotropic bench -t 20 >> bench3d.out

echo "" >> bench3d.out
echo "3D anisotropic convolution, t = 22" >> bench3d.out
./bench3d_convolution_anisotropic bench -t 22 >> bench3d.out

echo "" >> bench3d.out
echo "3D anisotropic convolution, t = 24" >> bench3d.out
./bench3d_convolution_anisotropic bench -t 24 >> bench3d.out

