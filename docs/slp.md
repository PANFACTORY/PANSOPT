# Sequential Linear Programming

$$\begin{aligned}
\text{min}&\,F\left(\boldsymbol{x}\right)\\
\text{s.t.}&\,\begin{cases}
    {G_j\left(\boldsymbol{x}\right)\le 0}\\
    {L_i\le x_i\le U_i}
\end{cases}
\end{aligned}$$

Generate subproblem below with Taylor expansion to the first order term, and solve it until match convergence criteria.

$$\begin{aligned}
\text{min}&\,F\left(\boldsymbol{x}^k\right)+\nabla F\left(\boldsymbol{x}^k\right)\cdot\Delta \boldsymbol{x}\\
\text{s.t.}&\,\begin{cases}
    {G_j\left(\boldsymbol{x}^k\right)+\nabla G_j\left(\boldsymbol{x}^k\right)\cdot\Delta \boldsymbol{x}\le 0}\\
    {\text{max}\left(L_i-x^k_i,-\text{move}\right)\le\Delta x_i\le\text{min}\left(U_i-x^k_i,\text{move}\right)}
\end{cases}
\end{aligned}$$
