# Quasi-Newton method

## L-BFGS algorithm

Based on BFGS algorithm, $\boldsymbol{H}_k$ is defined as bellow.

$$\boldsymbol{H}_{k+1}=\left(\boldsymbol{I}-\frac{\boldsymbol{s}_k\boldsymbol{y}_k^T}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\right)^T\boldsymbol{H}_k\left(\boldsymbol{I}-\frac{\boldsymbol{s}_k\boldsymbol{y}_k^T}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\right)+\frac{\boldsymbol{s}_k\boldsymbol{s}_k^T}{\boldsymbol{y}_k^T\boldsymbol{s}_k}$$

Define $\boldsymbol{V}_k$ and $\rho_k$ as below respectively.

$$\begin{aligned}
    \boldsymbol{V}_k&=\left(\boldsymbol{I}-\frac{\boldsymbol{s}_k\boldsymbol{y}_k^T}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\right) \\
    \rho_k&=\frac{1}{\boldsymbol{y}_k^T\boldsymbol{s}_k}
\end{aligned}$$

Then,

$$\begin{aligned}
    \boldsymbol{H}_{k+1}=&\boldsymbol{V}_k^T\boldsymbol{H}_k\boldsymbol{V}
    +\rho_k\boldsymbol{s}_k\boldsymbol{s}_k^T \\
    =&\boldsymbol{V}_k^T\boldsymbol{V}_{k-1}^T\boldsymbol{H}_{k-1}\boldsymbol{V}_{k-1}\boldsymbol{V}_k+\rho_{k-1}\boldsymbol{V}_k\boldsymbol{s}_{k-1}\boldsymbol{s}_{k-1}^T\boldsymbol{V}_k+\rho_k\boldsymbol{s}_k\boldsymbol{s}_k^T \\
    &\vdots \\
    =&\left(\boldsymbol{V}_0\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right)^T\boldsymbol{H}_0\left(\boldsymbol{V}_0\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right) \\
    &+\rho_0\left(\boldsymbol{V}_1\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right)^T\boldsymbol{s}_0\boldsymbol{s}_0^T\left(\boldsymbol{V}_1\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right) \\
    &+\rho_1\left(\boldsymbol{V}_2\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right)^T\boldsymbol{s}_1\boldsymbol{s}_1^T\left(\boldsymbol{V}_2\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right) \\
    &+\cdots+\rho_m\left(\boldsymbol{V}_{m+1}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right)^T\boldsymbol{s}_m\boldsymbol{s}_m^T\left(\boldsymbol{V}_{m+1}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right) \\
    &+\cdots+\rho_{k-1}\boldsymbol{V}_k^T\boldsymbol{s}_{k-1}\boldsymbol{s}_{k-1}^T\boldsymbol{V}_k+\rho_k\boldsymbol{s}_k\boldsymbol{s}_k^T
\end{aligned}$$

To reduce memory, approximate $\boldsymbol{H}_k$ as bellow.

$$\begin{aligned}
    \boldsymbol{H}_{k+1}=&\left(\boldsymbol{V}_{k-m}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right)^T\boldsymbol{H}_0\left(\boldsymbol{V}_{k-m}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right) \\
    &+\rho_{k-m}\left(\boldsymbol{V}_{k-m+1}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right)^T\boldsymbol{s}_{k-m}\boldsymbol{s}_{k-m}^T\left(\boldsymbol{V}_{k-m+1}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right) \\
    &+\rho_{k-m+1}\left(\boldsymbol{V}_{k-m+2}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right)^T\boldsymbol{s}_{k-m+1}\boldsymbol{s}_{k-m+1}^T\left(\boldsymbol{V}_{k-m+2}\cdots\boldsymbol{V}_{k-1}\boldsymbol{V}_k\right) \\
    &+\cdots+\rho_{k-1}\boldsymbol{V}_k^T\boldsymbol{s}_{k-1}\boldsymbol{s}_{k-1}^T\boldsymbol{V}_k+\rho_k\boldsymbol{s}_k\boldsymbol{s}_k^T
\end{aligned}$$

Search direction is $\boldsymbol{d}=\boldsymbol{H}_k\boldsymbol{g}_k$, and

$$\begin{aligned}
\boldsymbol{V}_k\boldsymbol{a}=&\left(\boldsymbol{I}-\frac{\boldsymbol{s}_k\boldsymbol{y}_k^T}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\right)\boldsymbol{a} \\
=&\boldsymbol{a}-\frac{\boldsymbol{y}_k^T\boldsymbol{a}}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\boldsymbol{s}_k
\end{aligned}$$

$$\begin{aligned}
\boldsymbol{V}_k^T\boldsymbol{a}=&\left(\boldsymbol{I}-\frac{\boldsymbol{s}_k\boldsymbol{y}_k^T}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\right)^T\boldsymbol{a} \\
=&\left(\boldsymbol{I}-\frac{\boldsymbol{y}_k\boldsymbol{s}_k^T}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\right)^T\boldsymbol{a} \\
=&\boldsymbol{a}-\frac{\boldsymbol{s}_k^T\boldsymbol{a}}{\boldsymbol{y}_k^T\boldsymbol{s}_k}\boldsymbol{y}_k
\end{aligned}$$


[ref](https://ysk24ok.github.io/2017/03/31/lbfgs.html)
