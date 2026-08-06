# 2026.08.06
## Standard procedure of transformation
Genz (1992) presented a method for computing integrals, in which the change‑of‑variables step is particularly insightful. Applying this matrix decomposition and transformation approach to our design problem can simplify the integral to some extent. We first introduce Genz's standard transformation procedure.

A problem that arises in many statistics applications is that of computing the multi-variate normal distribution function

$$
F(\mathbf{a},\mathbf{b}) = \frac{1}{\sqrt{|\Sigma|(2\pi)^m}}\int_{a_1}^{b_1}\int_{a_2}^{b_2}\cdots\int_{a_m}^{b_m}e^{-\frac{1}{2}\boldsymbol{\theta}^T\Sigma^{-1}\boldsymbol{\theta}}d\boldsymbol{\theta},
$$

where $\boldsymbol{\theta} = (\boldsymbol{\theta_1}, \boldsymbol{\theta_2}, \ldots, \boldsymbol{\theta_m})^T$ and $\boldsymbol{\Sigma}$ is an $m \times m$ symmetric positive definite covariance matrix.

According to the Cholesky decomposition, a positive definite covariance matrix can be decomposed as $\mathbf{\Sigma} = \mathbf{C}\mathbf{C^T}$ ,where $\mathbf{C}$ is a lower triangular matrix. Applying the Cholesky transformation $\boldsymbol{\theta} = \mathbf{C}\mathbf{y}$ ,we have $\boldsymbol{\theta}^T \Sigma^{-1} \boldsymbol{\theta} = \mathbf{y}^T\mathbf{C}^T\(\mathbf{C}^T\)^{-1}\mathbf{C}^{-1}\mathbf{C}\mathbf{y} = \mathbf{y}^T\mathbf{y}$ ,and $d\boldsymbol{\theta} = |\mathbf{C}|d\mathbf{y} = |\mathbf{\Sigma}|^{-\frac{1}{2}}d\mathbf{y}$ .Because $\quad\mathbf{a} \le \theta = \mathbf{C}\mathbf{y} \le \mathbf{b}\quad$ implies $\quad\left(a_i-\sum_{j=1}^{i-1}c_{ij}y_j\right)/c_{ii}\le y_i\le\left(b_i-\sum_{j=1}^{i-1}c_{ij}y_j\right)/c_{ii}\quad$ for $\quad i = 1, 2, \dots, m\quad$ ,we have

$$
F(\mathbf{a}, \mathbf{b}) = \frac{1}{\sqrt{(2\pi)^m}} \int_{a'_1}^{b'_1} e^{-\frac{y_1^2}{2}} \int_{a'_2(y_1)}^{b'_2(y_1)} e^{-\frac{y_2^2}{2}} \dots \int_{a'_m(y_1, \dots, y_{m-1})}^{b'_m(y_1, \dots, y_{m-1})} e^{-\frac{y_m^2}{2}} d\mathbf{y},
$$

with $a'\_i(y_1,\dots,y_{i-1}) = \left(a_i-\sum_{j=1}^{i-1}c_{ij}y_j\right)/c_{ii}\quad$ and $\quad b'\_i(y_1,\dots, y_{i-1}) = \left(b_i-\sum_{j=1}^{i-1}c_{ij}y_j\right)/c_{ii}.$

## Application of the transformation to SSD design
To obtain a zero-mean vector, we define $\mathbf{u} = -\frac{1}{\sqrt{n}}\mathbf{Z}\_{\mathcal{A}}\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{F}\_{\mathcal{A}}^T\mathbf{e}$ .Then, $\mathbf{u}$ follows a k-dimensional normal distribution, denoted as $\mathbf{u} \sim N\_k(\mathbf{0} , \mathbf{Z}\_{\mathcal{A}}\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{Z}\_{\mathcal{A}})$ , and the probability $P(S\_{\lambda})$ can be expressed as

$$
P(S\_{\lambda}) = P(\mathbf{u} < \sqrt{n}\mathbf{Z}\_{\mathcal{A}}(\mathbf{V}\_{\mathcal{A}}^{\frac{1}{2}}\boldsymbol{\beta}_{\mathcal{A}}-\lambda\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{z}\_{\mathcal{A}})).
$$

Denote $\mathbf{Z}\_{\mathcal{A}}\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{Z}\_{\mathcal{A}}$ by $\mathbf{\Sigma}$ ,and denote $\sqrt{n}\mathbf{Z}\_{\mathcal{A}}(\mathbf{V}\_{\mathcal{A}}^{\frac{1}{2}}\boldsymbol{\beta}\_{\mathcal{A}}-\lambda\mathbf{C}\_{\mathcal{A}}^{-1}\mathbf{z}\_{\mathcal{A}})$ by $\mathbf{v} = (v\_1, v\_2, \ldots, v\_k)^T$ .Then 

$$
P(S\_{\lambda}) = P(\mathbf{u} < \mathbf{v}) = \frac{1}{\sqrt{|\Sigma|(2\pi)^k}}\int\_{-\infty}^{v\_1}\int\_{-\infty}^{v\_2}\cdots\int\_{-\infty}^{v\_k}e^{-\frac{1}{2}\mathbf{u}^T\Sigma^{-1}\mathbf{u}}d\mathbf{u}.
$$

Applying the Cholesky decomposition to $\mathbf{\Sigma}$ ,we have $\mathbf{\Sigma} = \mathbf{\Gamma}^{-1}(\mathbf{\Gamma}^{-1})^T$ ,where both $\mathbf{\Gamma}^{-1}$ and $\mathbf{\Gamma}$ are lower triangular matrices and $\mathbf{\Gamma} = (\gamma\_{ij})\_{k \times k}$ . Let $\boldsymbol{\xi} = \mathbf{\Gamma}\mathbf{u}$ ,then $\boldsymbol{\xi}$ follows a k-dimensional standard normal distribution, i.e., $\boldsymbol{\xi} \sim N\_k(\mathbf{0},\mathbf{I}\_{k}).$ Now $\mathbf{u}^T\Sigma^{-1}\mathbf{u} = \boldsymbol{\xi}^T(\mathbf{\Gamma}^T)^{-1}\mathbf{\Gamma}^T\mathbf{\Gamma}\mathbf{\Gamma}^{-1}\boldsymbol{\xi} = \boldsymbol{\xi}^T\boldsymbol{\xi}$ ,and $d\boldsymbol{\xi} = |\mathbf{\Gamma}|d\mathbf{u} = |\mathbf{\Sigma}|^{-\frac{1}{2}}d\mathbf{u}$ .Because $\quad\mathbf{u} = \mathbf{\Gamma}^{-1}\boldsymbol{\xi}\le\mathbf{v}$ implies $\xi\_i\le\left(v\_i-\sum\_{j=1}^{i-1}\gamma\_{ij}\xi\_j\right)/\gamma\_{ii}\quad$ for $\quad i = 1,2,\dots,m\quad$ ,we have

$$
P(S\_{\lambda}) = \frac{1}{\sqrt{(2\pi)^k}}\int\_{-\infty}^{v'\_1}\int\_{-\infty}^{v'\_2}\cdots\int\_{-\infty}^{v'\_k}e^{-\frac{1}{2}\boldsymbol{\xi}^T\boldsymbol{\xi}}d\boldsymbol{\xi},
$$

with $\quad v'\_i(\xi\_1,\dots,\xi\_{i-1}) = \left(v\_i-\sum\_{j=1}^{i-1}\gamma\_{ij}\xi\_j\right)/\gamma\_{ii}.$
