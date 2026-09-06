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
