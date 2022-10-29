# Primal-Dual Interior Point Method（with Box constraint）

$$\begin{aligned}
\text{min}&\,\boldsymbol{c}\cdot\boldsymbol{x}\\  
\text{s.t.}&\,\begin{cases}
        {\boldsymbol{A}_e\boldsymbol{x}=\boldsymbol{b}_e}\\
        {\boldsymbol{A}_i\boldsymbol{x}+\boldsymbol{s}=\boldsymbol{b}_i}\\
        {\boldsymbol{x}+\boldsymbol{t}=\boldsymbol{u}}\\
        {x_i\ge 0,s_i\ge 0,t_i\ge 0}
    \end{cases}
\end{aligned}$$

Lagrangian is as follows.

$$L=\boldsymbol{c}\cdot\boldsymbol{x}-\boldsymbol{y}_e\cdot\left(\boldsymbol{A}_e\boldsymbol{x}-\boldsymbol{b}_e\right)-\boldsymbol{y}_i\cdot\left(\boldsymbol{A}_i\boldsymbol{x}+\boldsymbol{s}-\boldsymbol{b}_i\right)-\boldsymbol{z}_x\cdot\boldsymbol{x}-\boldsymbol{z}_s\cdot\boldsymbol{s}-\boldsymbol{z}_t\cdot\left(\boldsymbol{u}-\boldsymbol{x}-\boldsymbol{t}\right)$$

KKT condition is as follows.

$$\begin{cases}
    {\boldsymbol{c}-\boldsymbol{A}_e^T\boldsymbol{y}_e-\boldsymbol{A}_i^T\boldsymbol{y}_i-\boldsymbol{z}_x+\boldsymbol{z}_t=\boldsymbol{0}}\\
    {\boldsymbol{A}_e\boldsymbol{x}-\boldsymbol{b}_e=\boldsymbol{0}}\\
    {\boldsymbol{A}_i\boldsymbol{x}+\boldsymbol{s}-\boldsymbol{b}_i=\boldsymbol{0}}\\
    {-\boldsymbol{y}_i-\boldsymbol{z}_s=\boldsymbol{0}}\\
    {\boldsymbol{Z}_x\boldsymbol{x}=\boldsymbol{0}}\\
    {\boldsymbol{Z}_s\boldsymbol{s}=\boldsymbol{0}}\\
    {\boldsymbol{Z}_t\boldsymbol{t}=\boldsymbol{0}}\\
    {x_i\ge 0,s_i\ge 0,t_i\ge 0,z_{xi}\ge 0,z_{si}\ge 0,z_{ti}\ge 0}
\end{cases}$$

Set initial value of $\boldsymbol{x}$, $\boldsymbol{s}$, $\boldsymbol{t}$, $\boldsymbol{y}_e$, $\boldsymbol{y}_i$, $\boldsymbol{z}_x$, $\boldsymbol{z}_s$ and $\boldsymbol{z}_t$ as

$$\begin{cases}
    {\boldsymbol{x}=s_0\boldsymbol{e}}\\
    {\boldsymbol{s}=s_0\boldsymbol{e}}\\
    {\boldsymbol{t}=s_0\boldsymbol{e}}\\
    {\boldsymbol{y}_e=\boldsymbol{0}}\\
    {\boldsymbol{y}_i=\boldsymbol{0}}\\
    {\boldsymbol{z}_x=s_0\boldsymbol{e}}\\
    {\boldsymbol{z}_s=s_0\boldsymbol{e}}\\
    {\boldsymbol{z}_t=s_0\boldsymbol{e}}
\end{cases}$$

then, $s_0$ is a large number like $1000$ .

Residual vector is

$$\begin{cases}
    {\boldsymbol{r}_{dx}=\boldsymbol{c}-\boldsymbol{A}_e^T\boldsymbol{y}_e-\boldsymbol{A}_i^T\boldsymbol{y}_i-\boldsymbol{z}_x+\boldsymbol{z}_t}\\
    {\boldsymbol{r}_{ds}=-\boldsymbol{y}_i-\boldsymbol{z}_s}\\
    {\boldsymbol{r}_{dt}=\boldsymbol{u}-\boldsymbol{x}-\boldsymbol{t}}\\
    {\boldsymbol{r}_{pe}=\boldsymbol{A}_e\boldsymbol{x}-\boldsymbol{b}_e}\\
    {\boldsymbol{r}_{pi}=\boldsymbol{A}_i\boldsymbol{x}+\boldsymbol{s}-\boldsymbol{b}_i}\\
    {\boldsymbol{r}_{cx}=\sigma\mu\boldsymbol{e}-\boldsymbol{Z}_x\boldsymbol{x}}\\
    {\boldsymbol{r}_{cs}=\sigma\mu\boldsymbol{e}-\boldsymbol{Z}_s\boldsymbol{s}}\\
    {\boldsymbol{r}_{ct}=\sigma\mu\boldsymbol{e}-\boldsymbol{Z}_t\boldsymbol{t}}
\end{cases}$$

and solve bellow.

$$\begin{cases}
    {\boldsymbol{A}_e^Td\boldsymbol{y}_e+\boldsymbol{A}_i^Td\boldsymbol{y}_i+d\boldsymbol{z}_x-d\boldsymbol{z}_t}=\boldsymbol{r}_{dx}\\
    {d\boldsymbol{y}_i+d\boldsymbol{z}_s=\boldsymbol{r}_{ds}}\\
    {d\boldsymbol{x}+d\boldsymbol{t}=\boldsymbol{r}_{dt}}\\
    {-\boldsymbol{A}_ed\boldsymbol{x}=\boldsymbol{r}_{pe}}\\
    {-\boldsymbol{A}_id\boldsymbol{x}-d\boldsymbol{s}=\boldsymbol{r}_{pi}}\\
    {\boldsymbol{Z}_xd\boldsymbol{x}+\boldsymbol{X}d\boldsymbol{z}_x=\boldsymbol{r}_{cx}}\\
    {\boldsymbol{Z}_sd\boldsymbol{s}+\boldsymbol{S}d\boldsymbol{z}_s=\boldsymbol{r}_{cs}}\\
    {\boldsymbol{Z}_td\boldsymbol{t}+\boldsymbol{T}d\boldsymbol{z}_t=\boldsymbol{r}_{ct}}
\end{cases}$$

First, solve equation of $\boldsymbol{y}_e$ and $\boldsymbol{y}_i$ .

$$\left\(\begin{array}{c|c}
\boldsymbol{A}_e\left\(\boldsymbol{Z}_x+\boldsymbol{X}\boldsymbol{T}^{-1}\boldsymbol{Z}_t\right\)^{-1}\boldsymbol{A}_e^T &
\boldsymbol{A}_e\left\(\boldsymbol{Z}_x+\boldsymbol{X}\boldsymbol{T}^{-1}\boldsymbol{Z}_t\right\)^{-1}\boldsymbol{A}_i^T \\ 
\\hline  
\boldsymbol{A}_i\left\(\boldsymbol{Z}_x+\boldsymbol{X}\boldsymbol{T}^{-1}\boldsymbol{Z}_t\right\)^{-1}\boldsymbol{A}_e^T &
\boldsymbol{A}_i\left\(\boldsymbol{Z}_x+\boldsymbol{X}\boldsymbol{T}^{-1}\boldsymbol{Z}_t\right\)^{-1}\boldsymbol{A}_i^T+\boldsymbol{Z}_s\boldsymbol{S}
\end{array}\right\)
\left\(\begin{array}{c}
d\boldsymbol{y}_e \\ 
\\hline d\boldsymbol{y}_i
\end{array}\right\)=
\left\(\begin{array}{c}
-\boldsymbol{r}_{pe}-\boldsymbol{A}_e\left\(\boldsymbol{Z}_x+\boldsymbol{X}\boldsymbol{T}^{-1}\boldsymbol{Z}_t\right\)^{-1}\left[\boldsymbol{r}_{cx}-\boldsymbol{X}\left\\{\boldsymbol{r}_{dx}+\boldsymbol{T}^{-1}\left\(\boldsymbol{r}_{ct}-\boldsymbol{Z}_t\boldsymbol{r}_{dt}\right\)\right\\}\right] \\ 
\\hline
-\boldsymbol{r}_{pi}-\boldsymbol{A}_i\left\(\boldsymbol{Z}_x+\boldsymbol{X}\boldsymbol{T}^{-1}\boldsymbol{Z}_t\right\)^{-1}\left[\boldsymbol{r}_{cx}-\boldsymbol{X}\left\\{\boldsymbol{r}_{dx}+\boldsymbol{T}^{-1}\left\(\boldsymbol{r}_{ct}-\boldsymbol{Z}_t\boldsymbol{r}_{dt}\right\)\right\\}\right]-\boldsymbol{Z}_s^{-1}\left\(\boldsymbol{r}_{cs}-\boldsymbol{S}\boldsymbol{r}_{ds}\right\)
\end{array}\right\)$$

Next, solve equations of $\boldsymbol{x}$, $\boldsymbol{s}$, $\boldsymbol{t}$, $\boldsymbol{z}_t$, $\boldsymbol{z}_s$ and $\boldsymbol{z}_x$ .

$$\begin{cases}
    {d\boldsymbol{x}=\left\(\boldsymbol{Z}_x+\boldsymbol{X}\boldsymbol{T}^{-1}\boldsymbol{Z}_t\right\)^{-1}\left[\boldsymbol{r}_{cx}-\boldsymbol{X}\left\\{\boldsymbol{r}_{dx}-\boldsymbol{A}_e^Td\boldsymbol{y}_e-\boldsymbol{A}_i^Td\boldsymbol{y}_i+\boldsymbol{T}^{-1}\left\(\boldsymbol{r}_{ct}-\boldsymbol{W}\boldsymbol{r}_{dt}\right\)\right\\}\right]}\\
    {d\boldsymbol{s}=\boldsymbol{Z}_s^{-1}\left\\{\boldsymbol{r}_{cs}-\boldsymbol{S}\left\(\boldsymbol{r}_{ds}-d\boldsymbol{y}_i\right\)\right\\}}\\
    {d\boldsymbol{t}=\boldsymbol{r}_{dt}-d\boldsymbol{x}}\\
    {d\boldsymbol{z}_t=\boldsymbol{T}^{-1}\left\(\boldsymbol{r}_{ct}-\boldsymbol{Z}_td\boldsymbol{t}\right\)}\\
    {d\boldsymbol{z}_s=\boldsymbol{r}_{ds}-d\boldsymbol{y}_i}\\
    {d\boldsymbol{z}_x=\boldsymbol{X}^{-1}\left\(\boldsymbol{r}_{cx}-\boldsymbol{Z}_xd\boldsymbol{x}\right\)}
\end{cases}$$

Get step size $\alpha$ as follows.

$$\alpha=\text{min}\left\\{\tau\alpha_x,\tau\alpha_s,\tau\alpha_t,\tau\alpha_{zx},\tau\alpha_{zs},\tau\alpha_{zt},1\right\\}$$

and

$$\alpha_\*=\underset{d\*_i<0}{\text{min}}\left\\{-\frac{\*_i}{d\*_i}\right\\}.$$

Update $\boldsymbol{x}$, $\boldsymbol{s}$, $\boldsymbol{t}$, $\boldsymbol{y}_e$, $\boldsymbol{y}_i$, $\boldsymbol{z}_x$, $\boldsymbol{z}_s$ and $\boldsymbol{z}_t$ with

$$\boldsymbol{\*}+=\alpha d\boldsymbol{\*}.$$

Calculate above flow from residual vector again, until match optimal criteria.
