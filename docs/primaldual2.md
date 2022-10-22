# Primal-Dual Interior Point Method（with Box constraint）

$\text{min}\,\bm{c}\cdot\bm{x}$  
$\text{s.t.}\,\begin{cases}
        {\bm{A}_e\bm{x}=\bm{b}_e}\\
        {\bm{A}_i\bm{x}+\bm{s}=\bm{b}_i}\\
        {\bm{x}+\bm{t}=\bm{u}}\\
        {x_i\ge 0,s_i\ge 0,t_i\ge 0}
    \end{cases}$

Lagrangian is

$L=\bm{c}\cdot\bm{x}-\bm{y}_e\cdot\left(\bm{A}_e\bm{x}-\bm{b}_e\right)-\bm{y}_i\cdot\left(\bm{A}_i\bm{x}+\bm{s}-\bm{b}_i\right)-\bm{z}_x\cdot\bm{x}-\bm{z}_s\cdot\bm{s}-\bm{z}_t\cdot\left(\bm{u}-\bm{x}-\bm{t}\right)$ .

KKT condition is

$\begin{cases}
    {\bm{c}-\bm{A}_e^T\bm{y}_e-\bm{A}_i^T\bm{y}_i-\bm{z}_x+\bm{z}_t=\bm{0}}\\
    {\bm{A}_e\bm{x}-\bm{b}_e=\bm{0}}\\
    {\bm{A}_i\bm{x}+\bm{s}-\bm{b}_i=\bm{0}}\\
    {-\bm{y}_i-\bm{z}_s=\bm{0}}\\
    {\bm{Z}_x\bm{x}=\bm{0}}\\
    {\bm{Z}_s\bm{s}=\bm{0}}\\
    {\bm{Z}_t\bm{t}=\bm{0}}\\
    {x_i\ge 0,s_i\ge 0,t_i\ge 0,z_{xi}\ge 0,z_{si}\ge 0,z_{ti}\ge 0}
\end{cases}$ .

Set initial value of $\bm{x}$, $\bm{s}$, $\bm{t}$, $\bm{y}_e$, $\bm{y}_i$, $\bm{z}_x$, $\bm{z}_s$ and $\bm{z}_t$ as

$\begin{cases}
    {\bm{x}=s_0\bm{e}}\\
    {\bm{s}=s_0\bm{e}}\\
    {\bm{t}=s_0\bm{e}}\\
    {\bm{y}_e=\bm{0}}\\
    {\bm{y}_i=\bm{0}}\\
    {\bm{z}_x=s_0\bm{e}}\\
    {\bm{z}_s=s_0\bm{e}}\\
    {\bm{z}_t=s_0\bm{e}}
\end{cases}$

then, $s_0$ is a large number like $1000$ .

Residual vector is

$\begin{cases}
    {\bm{r}_{dx}=\bm{c}-\bm{A}_e^T\bm{y}_e-\bm{A}_i^T\bm{y}_i-\bm{z}_x+\bm{z}_t}\\
    {\bm{r}_{ds}=-\bm{y}_i-\bm{z}_s}\\
    {\bm{r}_{dt}=\bm{u}-\bm{x}-\bm{t}}\\
    {\bm{r}_{pe}=\bm{A}_e\bm{x}-\bm{b}_e}\\
    {\bm{r}_{pi}=\bm{A}_i\bm{x}+\bm{s}-\bm{b}_i}\\
    {\bm{r}_{cx}=\sigma\mu\bm{e}-\bm{Z}_x\bm{x}}\\
    {\bm{r}_{cs}=\sigma\mu\bm{e}-\bm{Z}_s\bm{s}}\\
    {\bm{r}_{ct}=\sigma\mu\bm{e}-\bm{Z}_t\bm{t}}
\end{cases}$

and solve

$\begin{cases}
    {\bm{A}_e^Td\bm{y}_e+\bm{A}_i^Td\bm{y}_i+d\bm{z}_x-d\bm{z}_t}=\bm{r}_{dx}\\
    {d\bm{y}_i+d\bm{z}_s=\bm{r}_{ds}}\\
    {d\bm{x}+d\bm{t}=\bm{r}_{dt}}\\
    {-\bm{A}_ed\bm{x}=\bm{r}_{pe}}\\
    {-\bm{A}_id\bm{x}-d\bm{s}=\bm{r}_{pi}}\\
    {\bm{Z}_xd\bm{x}+\bm{X}d\bm{z}_x=\bm{r}_{cx}}\\
    {\bm{Z}_sd\bm{s}+\bm{S}d\bm{z}_s=\bm{r}_{cs}}\\
    {\bm{Z}_td\bm{t}+\bm{T}d\bm{z}_t=\bm{r}_{ct}}
\end{cases}$ .

First, solve equation of $\bm{y}_e$ and $\bm{y}_i$ .

$\left(\begin{array}{c|c}
\bm{A}_e\left(\bm{Z}_x+\bm{X}\bm{T}^{-1}\bm{Z}_t\right)^{-1}\bm{A}_e^T &
\bm{A}_e\left(\bm{Z}_x+\bm{X}\bm{T}^{-1}\bm{Z}_t\right)^{-1}\bm{A}_i^T \\ \hline  
\bm{A}_i\left(\bm{Z}_x+\bm{X}\bm{T}^{-1}\bm{Z}_t\right)^{-1}\bm{A}_e^T &
\bm{A}_i\left(\bm{Z}_x+\bm{X}\bm{T}^{-1}\bm{Z}_t\right)^{-1}\bm{A}_i^T+\bm{Z}_s\bm{S}
\end{array}\right)
\left(\begin{array}{c}
d\bm{y}_e \\ \hline d\bm{y}_i
\end{array}\right)=
\left(\begin{array}{c}
-\bm{r}_{pe}-\bm{A}_e\left(\bm{Z}_x+\bm{X}\bm{T}^{-1}\bm{Z}_t\right)^{-1}\left[\bm{r}_{cx}-\bm{X}\left\{\bm{r}_{dx}+\bm{T}^{-1}\left(\bm{r}_{ct}-\bm{Z}_t\bm{r}_{dt}\right)\right\}\right] \\ \hline
-\bm{r}_{pi}-\bm{A}_i\left(\bm{Z}_x+\bm{X}\bm{T}^{-1}\bm{Z}_t\right)^{-1}\left[\bm{r}_{cx}-\bm{X}\left\{\bm{r}_{dx}+\bm{T}^{-1}\left(\bm{r}_{ct}-\bm{Z}_t\bm{r}_{dt}\right)\right\}\right]-\bm{Z}_s^{-1}\left(\bm{r}_{cs}-\bm{S}\bm{r}_{ds}\right)
\end{array}\right)$

Next, solve equations of $\bm{x}$, $\bm{s}$, $\bm{t}$, $\bm{z}_t$, $\bm{z}_s$ and $\bm{z}_x$ .

$\begin{cases}
    {d\bm{x}=\left(\bm{Z}_x+\bm{X}\bm{T}^{-1}\bm{Z}_t\right)^{-1}\left[\bm{r}_{cx}-\bm{X}\left\{\bm{r}_{dx}-\bm{A}_e^Td\bm{y}_e-\bm{A}_i^Td\bm{y}_i+\bm{T}^{-1}\left(\bm{r}_{ct}-\bm{W}\bm{r}_{dt}\right)\right\}\right]}\\
    {d\bm{s}=\bm{Z}_s^{-1}\left\{\bm{r}_{cs}-\bm{S}\left(\bm{r}_{ds}-d\bm{y}_i\right)\right\}}\\
    {d\bm{t}=\bm{r}_{dt}-d\bm{x}}\\
    {d\bm{z}_t=\bm{T}^{-1}\left(\bm{r}_{ct}-\bm{Z}_td\bm{t}\right)}\\
    {d\bm{z}_s=\bm{r}_{ds}-d\bm{y}_i}\\
    {d\bm{z}_x=\bm{X}^{-1}\left(\bm{r}_{cx}-\bm{Z}_xd\bm{x}\right)}
\end{cases}$

Get step size $\alpha$ as follows.

$\alpha=\text{min}\left\{\tau\alpha_x,\tau\alpha_s,\tau\alpha_t,\tau\alpha_{zx},\tau\alpha_{zs},\tau\alpha_{zt},1\right\}$

and

$\alpha_*=\underset{d*_i<0}{\text{min}}\left\{-\frac{*_i}{d*_i}\right\}$ .

Update $\bm{x}$, $\bm{s}$, $\bm{t}$, $\bm{y}_e$, $\bm{y}_i$, $\bm{z}_x$, $\bm{z}_s$ and $\bm{z}_t$ with

$\bm{*}+=\alpha d\bm{*}$ .

Calculate above flow from residual vector again, until match optimal criteria.
