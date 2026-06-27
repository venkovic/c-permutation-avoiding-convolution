***
**Algorithm:** Butterfly kernel of radix-8  
***
**Input:**  $x\in\mathbb{C}^n$, $n=8^t$  
**Output:** $x:=A_{8,n}^Tx$
***
**for** $q=1,\dots,t$  **do**  
$k:=8^q$     
$z_1:=x_{(b-1)k+j}$  
$z_2:=x_{(b-1)k+\ell+j}$  
$z_3:=x_{(b-1)k+2\ell+j}$  
$z_4:=x_{(b-1)k+3\ell+j}$  
$z_5:=x_{(b-1)k+4\ell+j}$  
$z_6:=x_{(b-1)k+5\ell+j}$  
$z_7:=x_{(b-1)k+6\ell+j}$  
$z_8:=x_{(b-1)k+7\ell+j}$  
$\tau_1:=z_1+z_5$, $\tau_2:=\omega_k^{j-1}(z_1-z_5)$  
$\tau_3:=z_2+z_6$, $\tau_4:=\omega_k^{j-1}(z_2-z_6)$  
$\tau_5:=z_3+z_7$, $\tau_6:=\omega_k^{j-1}(z_3-z_7)$  
$\tau_7:=z_4+z_8$, $\tau_8:=\omega_k^{j-1}(z_4-z_8)$  
$x_{(b-1)k+j}:=\tau_1+\tau_3+\tau_5+\tau_7$  
$x_{(b-1)k+\ell+j}:=\tau_2+a\cdot\tau_4-i\cdot\tau_6+b\cdot\tau_8$  
$x_{(b-1)k+2\ell+j}:=\omega_{k}^{2(j-1)}(\tau_1-i\cdot\tau_3-\tau_5+i\cdot\tau_7)$  
$x_{(b-1)k+3\ell+j}:=\omega_{k}^{2(j-1)}(\tau_2+b\cdot\tau_4+i\cdot\tau_6+a\cdot\tau_8)$  
$x_{(b-1)k+4\ell+j}:=\omega_{k}^{4(j-1)}(\tau_1-\tau_3+\tau_5-\tau_7)$  
$x_{(b-1)k+5\ell+j}:=\omega_{k}^{4(j-1)}(\tau_2-a\cdot\tau_4-i\cdot\tau_6-b\cdot\tau_8)$  
$x_{(b-1)k+6\ell+j}:=\omega_{k}^{6(j-1)}(\tau_1+i\cdot\tau_3-\tau_5-i\cdot\tau_7)$  
$x_{(b-1)k+7\ell+j}:=\omega_{k}^{6(j-1)}(\tau_2-b\cdot\tau_4+i\cdot\tau_6-a\cdot\tau_8)$  
***

N. Venkovic and H. Anzt ,(2025). Permutation-avoiding FFT-based convolution. arXiv preprint, arXiv:12506.12718.