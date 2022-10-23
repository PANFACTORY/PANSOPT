# Sequential Linear Programming

$\text{min}\,F\left(\bm{x}\right)$  
$\text{s.t.}\,\begin{cases}
    {G_j\left(\bm{x}\right)\le 0}\\
    {L_i\le x_i\le U_i}
\end{cases}$

Generate subproblem below with Taylor expansion to the first order term, and solve it until match convergence criteria.

$\text{min}\,F\left(\bm{x}^k\right)+\nabla F\left(\bm{x}^k\right)\cdot\Delta \bm{x}$  
$\text{s.t.}\,\begin{cases}
    {G_j\left(\bm{x}^k\right)+\nabla G_j\left(\bm{x}^k\right)\cdot\Delta \bm{x}\le 0}\\
    {\text{max}\left(L_i-x^k_i,-\text{move}\right)\le\Delta x_i\le\text{min}\left(U_i-x^k_i,\text{move}\right)}
\end{cases}$
