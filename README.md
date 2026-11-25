# Sensor placement for a simply supported beam under one sensor loss

A finite-grid computational study of modal recovery on a uniform Euler–Bernoulli beam. It compares the best layout with every sensor working to the best layout when any one sensor can fail. Both optimization problems are solved by complete enumeration with integer interval bounds.

## Results

Select $m+1$ displacement sensors from the 31 interior positions $x/L=i/32$, for $i=1,\ldots,31$, to estimate the first $m$ bending-mode amplitudes.

| Modes | Layouts examined | Failure-aware optimum, positions i/32 | All-sensors optimum, positions i/32 | Worst-case determinant ratio* | Retained all-sensors determinant |
| ---: | ---: | --- | --- | ---: | ---: |
| 2 | 4,495 | 8, 16, 24 | 9, 10, 22 | 97.40× | 84.92% |
| 3 | 31,465 | 6, 12, 20, 26 | 7, 15, 17, 25 | 7.46× | 76.85% |
| 4 | 169,911 | 5, 10, 16, 22, 27 | 5, 11, 14, 20, 27 | 6.38× | 80.72% |
| 5 | 736,281 | 4, 8, 13, 19, 24, 28 | 4, 10, 15, 17, 22, 28 | 8.30× | 65.42% |

\* Determinant of the remaining information matrix, minimized over the failed sensor, for the failure-aware layout divided by the same worst case for the all-sensors optimum. Mirror-equivalent layouts appear once in the table. The program identifies every optimum. All failure-aware optima are unique and symmetric about midspan in this study.

In the four-mode case, the failure-aware layout retains about 81% of the all-sensors information determinant while increasing the worst-case post-failure information determinant by 6.38 times. The determinant measures the inverse squared volume of the joint uncertainty ellipsoid, with common factors canceled. Uncertainty in individual modes depends on the eigenvalues and eigenvectors of the information matrix.

## Model and certificate

Let $A_{ik}=\sin(k\pi x_i/L)$, the displacement of mode $k$ at sensor $i$ with unit modal amplitude. Assume independent measurements with equal noise variance, known sensor locations, and a known missing channel. With common noise scaling omitted, the all-sensors information matrix is $A^\mathsf{T}A$. After deletion of sensor $j$, $A_{-j}$ is square and the determinant is $\det(A_{-j})^2$. The failure-aware objective maximizes the minimum of these determinants. The all-sensors objective maximizes their sum by Cauchy–Binet.

For any square set of mode samples, write $\theta_i=\pi x_i/L$. Factoring $\sin\theta_i$ from each row leaves the Chebyshev polynomials $U_{k-1}(\cos\theta_i)$. Their leading coefficients are $2^{k-1}$, so the absolute determinant is

$$
2^{m(m-1)/2}
\prod_i\sin\theta_i
\prod_{i<j}(\cos\theta_i-\cos\theta_j)
$$

for positions ordered from left to right. Half-angle identities, the cosine recurrence, and outward-rounded integer products bound every factor and minor at scale $10^{28}$. At each optimum the best lower bound exceeds the largest upper bound of any excluded layout. Any remaining mirror pair has equal objectives since reflection only permutes sensor rows and changes signs of mode columns. The guarantee applies to the specified discrete grid and idealized measurement model.

## Reproduce

Python 3.9+ with the standard library.

```sh
python design.py --self-test > results.json
python design.py --modes 4
```

The self-test compares the product identity with direct determinant expansion. The JSON reports both candidate layouts and the exact scaled integers `best_score_lower_scaled` and `excluded_layouts_upper_scaled`. Their difference divided by `score_scale` is a certified positive lower bound on the gap to any excluded layout. The printed decimal minor and Gram intervals are for display.

Background references are [Dynamics of Structures, Université libre de Bruxelles](https://structuraldynamics.ulb.be/Vibrations/0-lecture_notes.pdf) for the beam modes and [Optimal Sensor Placement in Structural Dynamics](https://opencourses.uth.gr/modules/document/file.php/MHXB124/Unit%202/OED_2.pdf) for information-matrix criteria.
