In the algorithms presented in Section~4 of Venkovic and Anzt (2025), for any given $x\in\mathbb{C}^n$, we need to evaluate $A_{r,n}x$, $\overline{ A_{r,n}}x$ and $A_{r,n}^Tx$.
Irrespective of the radix, these kernels are implemented as described in Algos. 1 and 2.
In what follows, we show detailed implementations of these kernels for radices 2, 4, and 8.  

<div align="center">
  <img src="Printing-of-tex-algorithms/algo01_butterfly-related-kernels.png">
</div>

<div align="center">
  <img src="Printing-of-tex-algorithms/algo02_transposed-butterfly-kernel.png">
</div>

## Radix-2 butterfly-related kernels

For the radix-2 case, we have

$$
B_{2,k}=
(F_2\otimes I_{k/2})\text{diag}\left(I_{k/2},\Omega_{2,k/2}\right)
$$

where 

$$
F_2=\begin{bmatrix}
     1 & 1 \\ 
     1 & -1
    \end{bmatrix}
$$

so that

$$
B_{2,k}=
\begin{bmatrix}
I_{k/2} &  \Omega_{2,k/2} \\
I_{k/2} & -\Omega_{2,k/2}
\end{bmatrix}.
$$

Then, for all $x\in\mathbb{C}^k$, as we let $z_i:=x_{(i-1)k/2+1:ik/2}$ for $i=1,2$, we obtain

$$
\begin{align}
B_{2,k}x=
\begin{bmatrix}
z_{1}+\Omega_{2,k/2}z_{2}\\
z_{1}-\Omega_{2,k/2}z_{2}\\
\end{bmatrix}.
\end{align}
$$

We then let $\tau:=\Omega_{2,k/2}z_{2}$ so that

$$
\begin{align}
B_{2,k}x=
\begin{bmatrix}
z_{1}+\tau\\
z_{1}-\tau
\end{bmatrix}.
\end{align}
$$

This leads to Algo. 3 for the computation of $x\mapsto A_{2,n}x$.  

<div align="center">
  <img src="Printing-of-tex-algorithms/algo03_butterfly-kernel-radix-2.png">
</div>

The radix-2 conjugate butterfly kernel relies on the expression

$$
\begin{align}
\overline{B_{2,k}}x=
\begin{bmatrix}
z_{1}+\overline{\Omega_{2,k/2}}z_{2}\\
z_{1}-\overline{\Omega_{2,k/2}}z_{2}\\
\end{bmatrix}
\end{align}
$$

so that, as we let $\tau:=\overline{\Omega_{2,k/2}}z_{2}$, we still have

$$
\begin{align}
\overline{B_{2,k}}x=
\begin{bmatrix}
z_{1}+\tau\\
z_{1}-\tau
\end{bmatrix}.
\end{align}
$$

This leads to Algo. 4 for the computation of $x\mapsto\overline{A_{2,n}}x$.  

<div align="center">
  <img src="Printing-of-tex-algorithms/algo04_conjugate-butterfly-kernel-radix-2.png">
</div>

On the other hand, we also have

$$
\begin{align}
B_{2,k}^Tx=
\begin{bmatrix}
z_{1}+z_{2}\hfill\\
\Omega_{2,k/2}z_{1}-\Omega_{2,k/2}z_{2}
\end{bmatrix}
\end{align}
$$
so that if we let $\tau:=z_2$, we obtain
$$
\begin{align}
B_{2,k}^Tx=
\begin{bmatrix}
z_{1}+\tau\hfill\\
\Omega_{2,k/2}(z_{1}-\tau)
\end{bmatrix}.
\end{align}
$$

This leads to Algo. 5 for the computation of $x\mapsto A_{2,n}^Tx$.  

<div align="center">
  <img src="Printing-of-tex-algorithms/algo05_transposed-butterfly-kernel-radix-2.png">
</div>

## Radix-4 butterfly-related kernels
For the radix-4 case, we have

$$
\begin{align}
B_{4,k}=
(F_4\otimes I_{k/4})\text{diag}\left(I_{k/4},\Omega_{4,k/4},\Omega_{4,k/4}^2,\Omega_{4,k/4}^3\right)
\end{align}
$$

where

$$
\begin{align}
F_4=
\begin{bmatrix}
1& 1& 1& 1\\
1&-i&-1& i\\
1&-1& 1&-1\\
1& i&-1&-i
\end{bmatrix}
\end{align}
$$

so that

$$
\begin{align}
B_{4,k}=
\begin{bmatrix}
I_{k/4}&       \Omega_{4,k/4}& \Omega_{4,k/4}^2&       \Omega_{4,k/4}^3\\
I_{k/4}&-i\cdot\Omega_{4,k/4}&-\Omega_{4,k/4}^2& i\cdot\Omega_{4,k/4}^3\\
I_{k/4}&      -\Omega_{4,k/4}& \Omega_{4,k/4}^2&      -\Omega_{4,k/4}^3\\
I_{k/4}& i\cdot\Omega_{4,k/4}&-\Omega_{4,k/4}^2&-i\cdot\Omega_{4,k/4}^3
\end{bmatrix}.
\end{align}
$$

Then, for all $x\in\mathbb{C}^k$, as we let $z_i:=x_{(i-1)k/4+1:ik/4}$ for $i=1,\dots,4$, we obtain

$$
\begin{align}
B_{4,k}x=
\begin{bmatrix}
z_{1}+\Omega_{4,k/4}z_{2}+\Omega_{4,k/4}^2z_{3}+\Omega_{4,k/4}^3z_{4}\hfill\\
z_{1}-i\cdot\Omega_{4,k/4}z_{2}-\Omega_{4,k/4}^2z_{3}+i\cdot\Omega_{4,k/4}^3z_{4}\\
z_{1}-\Omega_{4,k/4}z_{2}+\Omega_{4,k/4}^2z_{3}-\Omega_{4,k/4}^3z_{4}\hfill\\
z_{1}+i\cdot\Omega_{4,k/4}z_{2}-\Omega_{4,k/4}^2z_{3}-i\cdot\Omega_{4,k/4}^3z_{4}\hfill
\end{bmatrix}.
\end{align}
$$

Once we introduce

$$
\begin{align}
&\tau_1:=z_1+\Omega_{4,k/4}^2z_3\\
&\tau_2:=z_1-\Omega_{4,k/4}^2z_3\\
&\tau_3:=\Omega_{4,k/4}z_2+\Omega_{4,k/4}^3z_4\\
&\tau_4:=\Omega_{4,k/4}z_2-\Omega_{4,k/4}^3z_4\\
\end{align}
$$

we get

$$
\begin{align}
B_{4,k}x=
\begin{bmatrix}
\tau_1+\tau_3\\
\tau_2-i\cdot\tau_4\\
\tau_1-\tau_3\hfill\\
\tau_2+i\cdot\tau_4\hfill
\end{bmatrix}.
\end{align}
$$

This leads to Algo. 6 for the computation of $x\mapsto A_{4,n}x$.  

<div align="center">
  <img src="Printing-of-tex-algorithms/algo06_butterfly-kernel-radix-4.png">
</div>


The radix-4 conjugate butterfly kernel relies on the expression

$$
\begin{align}
\overline{B_{4,k}}x=
\begin{bmatrix}
z_{1}+\overline{\Omega_{4,k/4}}z_{2}+\overline{\Omega_{4,k/4}^2}z_{3}+\overline{\Omega_{4,k/4}^3}z_{4}\hfill\\
z_{1}+i\cdot\overline{\Omega_{4,k/4}}z_{2}-\overline{\Omega_{4,k/4}^2}z_{3}-i\cdot\overline{\Omega_{4,k/4}^3}z_{4}\\
z_{1}-\overline{\Omega_{4,k/4}}z_{2}+\overline{\Omega_{4,k/4}^2}z_{3}-\overline{\Omega_{4,k/4}^3}z_{4}\hfill\\
z_{1}-i\cdot\overline{\Omega_{4,k/4}}z_{2}-\overline{\Omega_{4,k/4}^2}z_{3}+i\cdot\overline{\Omega_{4,k/4}^3}z_{4}\hfill
\end{bmatrix}
\end{align}
$$

so that, as we introduce

$$
\begin{align}
&\tau_1:=z_1+\overline{\Omega_{4,k/4}^2}z_3\\
&\tau_2:=z_1-\overline{\Omega_{4,k/4}^2}z_3\\
&\tau_3:=\overline{\Omega_{4,k/4}}z_2+\overline{\Omega_{4,k/4}^3}z_4\\
&\tau_4:=\overline{\Omega_{4,k/4}}z_2-\overline{\Omega_{4,k/4}^3}z_4
\end{align}
$$

we get

$$
\begin{align}
\overline{B_{4,k}}x=
\begin{bmatrix}
\tau_1+\tau_3\hfill\\
\tau_2+i\cdot\tau_4\hfill\\
\tau_1-\tau_3\hfill\\
\tau_2-i\cdot\tau_4\hfill
\end{bmatrix}.
\end{align}
$$

This leads to Algo. 7 for the computation of $x\mapsto \overline{A_{4,n}}x$.

<div align="center">
  <img src="Printing-of-tex-algorithms/algo07_conjugate-butterfly-kernel-radix-4.png">
</div>

The radix-4 transposed butterfly kernel relies on the expression

$$
\begin{align}
B_{4,k}^Tx=
\begin{bmatrix}
z_{1}+z_{2}+z_{3}+z_{4}\hfill\\
\Omega_{4,k/4}z_{1}-i\cdot\Omega_{4,k/4}z_{2}-\Omega_{4,k/4}z_{3}+i\cdot\Omega_{4,k/4}z_{4}\\
\Omega_{4,k/4}^2z_{1}-\Omega_{4,k/4}^2z_{2}+\Omega_{4,k/4}^2z_{3}-\Omega_{4,k/4}^2z_{4}\hfill\\
\Omega_{4,k/4}^3z_{1}+i\cdot\Omega_{4,k/4}^3z_{2}-\Omega_{4,k/4}^3z_{3}-i\cdot\Omega_{4,k/4}^3z_{4}\hfill
\end{bmatrix}.
\end{align}
$$

As we introduce

$$
\begin{align}
&\tau_1:=z_1+z_3\\
&\tau_2:=\Omega_{4,k/4}(z_1-z_3)\\
&\tau_3:=z_2+z_4\\
&\tau_4:=\Omega_{4,k/4}(z_2-z_4)
\end{align}
$$

we obtain

$$
\begin{align}
B_{4,k}^Tx=
\begin{bmatrix}
\tau_1+\tau_3\hfill\\
\tau_2-i\cdot\tau_4\hfill\\
\Omega_{4,k/4}^2(\tau_1-\tau_3)\hfill\\
\Omega_{4,k/4}^2(\tau_2+i\cdot\tau_4)\hfill
\end{bmatrix}.
\end{align}
$$

This leads to Algo. 8 for the computation of $x\mapsto A_{4,n}^Tx$.

<div align="center">
  <img src="Printing-of-tex-algorithms/algo08_transposed-butterfly-kernel-radix-4.png">
</div>

## Radix-8 butterfly-related kernels

For the radix-8 case, we have

$$
\begin{align}
B_{8,k}=
(F_8\otimes I_{k/8})\text{diag}\left(I_{k/8},\Omega_{8,k/8},\Omega_{8,k/8}^2,\Omega_{8,k/8}^3,\Omega_{8,k/8}^4,\Omega_{8,k/8}^5,\Omega_{8,k/8}^6,\Omega_{8,k/8}^7\right)
\end{align}
$$

where

$$
\begin{align}
F_8=
\begin{bmatrix}
1& 1& 1& 1& 1& 1& 1& 1\\
1& a&-i& b&-1&-a& i&-b\\
1&-i&-1& i& 1&-i&-1& i\\
1& b& i& a&-1&-b&-i&-a\\
1&-1& 1&-1& 1&-1& 1&-1\\
1&-a&-i&-b&-1& a& i& b\\
1& i&-1&-i& 1& i&-1&-i\\
1&-b& i&-a&-1& b&-i& a
\end{bmatrix}
\end{align}
$$

in which $a=(1-i)/\sqrt{2}$ and $b=-(1+i)/\sqrt{2}=-ia$, so that, proceeding similarly as for previous radices, if we introduce:

$$
\begin{align}
&\tau_1:=z_1+\Omega_{8,k/8}^4z_5\\
&\tau_2:=z_1-\Omega_{8,k/8}^4z_5\\
&\tau_3:=\Omega_{8,k/8}z_2+\Omega_{8,k/8}^5z_6\\
&\tau_4:=\Omega_{8,k/8}z_2-\Omega_{8,k/8}^5z_6\\
&\tau_5:=\Omega_{8,k/8}^2z_3+\Omega_{8,k/8}^6z_7\\
&\tau_6:=\Omega_{8,k/8}^2z_3-\Omega_{8,k/8}^6z_7\\
&\tau_7:=\Omega_{8,k/8}^3z_4+\Omega_{8,k/8}^7z_8\\
&\tau_8:=\Omega_{8,k/8}^3z_4-\Omega_{8,k/8}^7z_8
\end{align}
$$

we obtain

$$
\begin{align}
B_{8,k}x=
\begin{bmatrix}
\tau_1+\tau_3+\tau_5+\tau_7\hfill\\
\tau_2+a\cdot\tau_4-i\cdot\tau_6+b\cdot\tau_8\hfill\\
\tau_1-i\cdot\tau_3-\tau_5+i\cdot\tau_7\hfill\\
\tau_2+b\cdot\tau_4+i\cdot\tau_6+a\cdot\tau_8\hfill\\
\tau_1-\tau_3+\tau_5-\tau_7\hfill\\
\tau_2-a\cdot\tau_4-i\cdot\tau_6-b\cdot\tau_8\hfill\\
\tau_1+i\cdot\tau_3-\tau_5-i\cdot\tau_7\hfill\\
\tau_2-b\cdot\tau_4+i\cdot\tau_6-a\cdot\tau_8\hfill
\end{bmatrix}.
\end{align}
$$

This leads to Algo. 9 for the computation of $x\mapsto A_{8,n}x$.

<div align="center">
  <img src="Printing-of-tex-algorithms/algo09_butterfly-kernel-radix-8.png">
</div>

For the radix-8 conjugate butterfly kernel, we introduce

$$
\begin{align}
&\tau_1:=z_1+\overline{\Omega_{8,k/8}^4}z_5\\
&\tau_2:=z_1-\overline{\Omega_{8,k/8}^4}z_5\\
&\tau_3:=\overline{\Omega_{8,k/8}}z_2+\overline{\Omega_{8,k/8}^5}z_6\\
&\tau_4:=\overline{\Omega_{8,k/8}}z_2-\overline{\Omega_{8,k/8}^5}z_6\\
&\tau_5:=\overline{\Omega_{8,k/8}^2}z_3+\overline{\Omega_{8,k/8}^6}z_7\\
&\tau_6:=\overline{\Omega_{8,k/8}^2}z_3-\overline{\Omega_{8,k/8}^6}z_7\\
&\tau_7:=\overline{\Omega_{8,k/8}^3}z_4+\overline{\Omega_{8,k/8}^7}z_8\\
&\tau_8:=\overline{\Omega_{8,k/8}^3}z_4-\overline{\Omega_{8,k/8}^7}z_8
\end{align}
$$

and we obtain

$$
\begin{align}
\overline{B_{8,k}}x=
\begin{bmatrix}
\tau_1+\tau_3+\tau_5+\tau_7\hfill\\
\tau_2+\overline{a}\cdot\tau_4+i\cdot\tau_6+\overline{b}\cdot\tau_8\hfill\\
\tau_1+i\cdot\tau_3-\tau_5-i\cdot\tau_7\hfill\\
\tau_2+\overline{b}\cdot\tau_4-i\cdot\tau_6+\overline{a}\cdot\tau_8\hfill\\
\tau_1-\tau_3+\tau_5-\tau_7\hfill\\
\tau_2-\overline{a}\cdot\tau_4+i\cdot\tau_6-\overline{b}\cdot\tau_8\hfill\\
\tau_1-i\cdot\tau_3-\tau_5+i\cdot\tau_7\hfill\\
\tau_2-\overline{b}\cdot\tau_4-i\cdot\tau_6-\overline{a}\cdot\tau_8\hfill
\end{bmatrix}
\end{align}
$$

which leads to Algo. 10 for the computation of $x\mapsto \overline{A_{8,n}}x$.

<div align="center">
  <img src="Printing-of-tex-algorithms/algo10_conjugate-butterfly-kernel-radix-8.png">
</div>

For the transposed butterfly kernel, we introduce

$$
\begin{align}
&\tau_1:=z_1+z_5\\
&\tau_2:=\Omega_{8,k/8}(z_1-z_5)\\
&\tau_3:=z_2+z_6\\
&\tau_4:=\Omega_{8,k/8}(z_2-z_6)\\
&\tau_5:=z_3+z_7\\
&\tau_6:=\Omega_{8,k/8}(z_3-z_7)\\
&\tau_7:=z_4+z_8\\
&\tau_8:=\Omega_{8,k/8}(z_4-z_8)
\end{align}
$$

to obtain

$$
\begin{align}
B_{8,k}^Tx=
\begin{bmatrix}
\tau_1+\tau_3+\tau_5+\tau_7\hfill\\
\tau_2+a\cdot\tau_4-i\cdot\tau_6+b\cdot\tau_8\hfill\\
\Omega_{8,k/8}^2(\tau_1-i\cdot\tau_3-\tau_5+i\cdot\tau_7)\hfill\\
\Omega_{8,k/8}^2(\tau_2+b\cdot\tau_4+i\cdot\tau_6+a\cdot\tau_8)\hfill\\
\Omega_{8,k/8}^4(\tau_1-\tau_3+\tau_5-\tau_7)\hfill\\
\Omega_{8,k/8}^4(\tau_2-a\cdot\tau_4-i\cdot\tau_6-b\cdot\tau_8)\hfill\\
\Omega_{8,k/8}^6(\tau_1+i\cdot\tau_3-\tau_5-i\cdot\tau_7)\hfill\\
\Omega_{8,k/8}^6(\tau_2-b\cdot\tau_4+i\cdot\tau_6-a\cdot\tau_8)\hfill
\end{bmatrix}
\end{align}
$$

and that leads to Algo. 11 for the computation of $x\mapsto A_{8,n}^Tx$.

<div align="center">
  <img src="Printing-of-tex-algorithms/algo11_transposed-butterfly-kernel-radix-8.png">
</div>

N. Venkovic and H. Anzt (2025). Permutation-avoiding FFT-based convolution. arXiv preprint, arXiv:12506.12718.