# 2026.09.06
## Cholesky decomposition of SSD design
To obtain a zero-mean vector, we define $\mathbf{u} = -\frac{1}{\sqrt{n}}\mathbf{Z}\_{\mathcal{A}}\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{F}\_{\mathcal{A}}^T\mathbf{e}$ . Then, $\mathbf{u}$ follows a k-dimensional normal distribution, denoted as $\mathbf{u} \sim N\_k(\mathbf{0} , \mathbf{Z}\_{\mathcal{A}}\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{Z}\_{\mathcal{A}})$ , and the probability $P(S\_{\lambda})$ can be expressed as

$$
P(S\_{\lambda}) = P\left(\mathbf{u} < \sqrt{n}\mathbf{Z}\_{\mathcal{A}}\left(\mathbf{V}\_{\mathcal{A}}^{\frac{1}{2}}\boldsymbol{\beta}_{\mathcal{A}}-\lambda\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{z}\_{\mathcal{A}}\right)\right).
$$

Denote $\mathbf{Z}\_{\mathcal{A}}\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{Z}\_{\mathcal{A}}$ by $\mathbf{\Sigma}$ ,and denote $\sqrt{n}\mathbf{Z}\_{\mathcal{A}}(\mathbf{V}\_{\mathcal{A}}^{\frac{1}{2}}\boldsymbol{\beta}\_{\mathcal{A}}-\lambda\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{z}\_{\mathcal{A}})$ by $\mathbf{v} = (v\_1, v\_2, \ldots, v\_k)^T$ . Then 

$$
P(S\_{\lambda}) = P(\mathbf{u} < \mathbf{v}) = \frac{1}{\sqrt{|\mathbf{\Sigma}|(2\pi)^k}}\int\_{-\infty}^{v\_1}\int\_{-\infty}^{v\_2}\cdots\int\_{-\infty}^{v\_k}e^{-\frac{1}{2}\mathbf{u}^T\mathbf{\Sigma}^{-1}\mathbf{u}}d\mathbf{u}.
$$

Applying the Cholesky decomposition to $\mathbf{\Sigma}$ , we have $\mathbf{\Sigma} = \mathbf{\Gamma}^{-1}(\mathbf{\Gamma}^{-1})^T$ , where both $\mathbf{\Gamma}^{-1}$ and $\mathbf{\Gamma}$ are lower triangular matrices and $\mathbf{\Gamma} = (\gamma\_{ij})\_{k \times k}$ . Let $\boldsymbol{\xi} = \mathbf{\Gamma}\mathbf{u}$ , then $\boldsymbol{\xi}$ follows a k-dimensional standard normal distribution, i.e., $\boldsymbol{\xi} \sim N\_k(\mathbf{0},\mathbf{I}\_{k}).$ Now $\mathbf{u}^T\Sigma^{-1}\mathbf{u} = \boldsymbol{\xi}^T(\mathbf{\Gamma}^T)^{-1}\mathbf{\Gamma}^T\mathbf{\Gamma}\mathbf{\Gamma}^{-1}\boldsymbol{\xi} = \boldsymbol{\xi}^T\boldsymbol{\xi}$ , and $d\boldsymbol{\xi} = |\mathbf{\Gamma}|d\mathbf{u} = |\mathbf{\Sigma}|^{-\frac{1}{2}}d\mathbf{u}$ .Because $\quad\mathbf{u} = \mathbf{\Gamma}^{-1}\boldsymbol{\xi}\le\mathbf{v}$ implies $\xi\_i\le\left(v\_i-\sum\_{j=1}^{i-1}\gamma\_{ij}\xi\_j\right)/\gamma\_{ii}\quad$ for $\quad i = 1,2,\dots,m\quad$ , we have

$$
P(S\_{\lambda}) = \frac{1}{\sqrt{(2\pi)^k}}\int\_{-\infty}^{v'\_1}\int\_{-\infty}^{v'\_2}\cdots\int\_{-\infty}^{v'\_k}e^{-\frac{1}{2}\boldsymbol{\xi}^T\boldsymbol{\xi}}d\boldsymbol{\xi} = \frac{1}{\sqrt{(2\pi)^k}} \int_{-\infty}^{v'_1} e^{-\frac{\xi_1^2}{2}} \int_{-\infty}^{v'_2(\xi_1)} e^{-\frac{\xi_2^2}{2}} \dots \int_{-\infty}^{v'_k(\xi_1, \dots, \xi_{k-1})} e^{-\frac{\xi_k^2}{2}} d\boldsymbol{\xi},
$$

with $\quad v'\_i(\xi\_1,\dots,\xi\_{i-1}) = \left(v\_i-\sum\_{j=1}^{i-1}\gamma\_{ij}\xi\_j\right)/\gamma\_{ii}.$

## Subsequent treatment
Since the density has been reduced to the standard normal distribution, we focus on the integration region. The $i$ -th integration upper limit $v_i$ depends on $\xi_1, \xi_2, \dots, \xi_{i-1}$ , and the elements of $\Gamma$ are not independent of the components of $\mathbf{v}$ , which makes the integration upper limits rather complicated. In the nested integration of Genz (1992), we have $a'\_i(y_1,\dots,y_{i-1}) = \left(a_i-\sum_{j=1}^{i-1}c_{ij}y_j\right)/c_{ii}\quad$ and $\quad b'\_i(y_1,\dots, y_{i-1}) = \left(b_i-\sum_{j=1}^{i-1}c_{ij}y_j\right)/c_{ii}.$ We note that the length of the $i$ -th integration interval is $b'\_i(y_1,\dots, y_{i-1})-a'\_i(y_1,\dots,y_{i-1}) = (b\_i-a\_i)/c\_{ii}$ , which is independent of the first $i-1$ random variables, thus convenient to handle. In our problem, $b\_i = v\_i, c\_{ij} = \gamma\_{ij},$ and all $a\_i$ are identically $-\infty$ ; therefore, we initially thought that we might only consider the quantities $\frac{b\_i}{\gamma\_{ii}}$ . In fact, the integration interval length for each component here tends to $\infty$ , so our idea of examining the interval length is merely an intuition. We vaguely felt that when the lower bound is very small and approaches $-\infty$ , the interval length almost determines the integral value, but the numerical simulation results do not support this. Due to the complexity of the integration upper limits, a theoretically strictly equivalent objective function may not be a problem that can be tackled in the short term. At present, we have a simple optimization idea:

(1) First, we intend to prove that, with other $v_j$ 's fixed, the probability (integral value) $P(S_{\lambda})$ is monotonically increasing with respect to $\frac{v_i}{\gamma_{ii}}$ .

(2) A heuristic optimization objective we plan to examine is $\sum_{i=1}^k{\frac{v_i}{\gamma_{ii}}}$ .

(3) We could then apply Genz's algorithm within a certain range to further compare the performance of various designs. For example, when a design satisfies $\sum_{i=1}^k{\frac{v_i}{\gamma_{ii}}} > M$ , where $M$ is a prespecified threshold, we would retain this design. The class of all designs meeting this threshold is denoted by $\mathcal{M}$ .

(4) For the design class $\mathcal{M}$ , we would apply Genz's algorithm to select the design that maximizes the integral value, and take it as the final optimal design.  

We will continue to explore this idea while also attempting to find other more analytical and effective optimization methods. 
