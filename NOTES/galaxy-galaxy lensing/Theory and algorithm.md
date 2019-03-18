<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

Hubbule parameter:

$$ \begin{aligned}H^2 &= H_0^2[\frac{\Omega_r}{a^4}+\frac{\Omega_m}{a^3}-\frac{Kc^2}{a^2H_0^2}+\Omega_{\Lambda}] \\\\
&=H_0^2[\frac{\Omega_r}{a^4}+\frac{\Omega_m}{a^3}+\frac{1-\Omega_0}{a^2}+\Omega_{\Lambda}] \end{aligned}$$

Comving distance:

$$ dt = \frac{da}{\dot{a}} \Rightarrow -dw = \frac{cdt}{a} = \frac{cda}{a\dot{a}}=\frac{cda}{a^2H}$$

$$ \begin{aligned}w(z_1,z_2) &= \frac{c}{H_0}\int_{a(z_2)}^{a(z_1)} \frac{da}{\sqrt{a\Omega_m + a^2(1-\Omega_m -\Omega_\Lambda) + a^4\Omega_\Lambda}}, z_1 < z_2 \\\\ 
&= \frac{c}{H_0}\int_{z_1}^{z_2} \frac{dz}{\sqrt{(1+z)^3\Omega_m + (1+z)^2(1-\Omega_m -\Omega_\Lambda) + \Omega_\Lambda}}, z_1 < z_2 \\\\
&= \frac{c}{H_0}\alpha(z_1, z_2) \end{aligned}$$

The search radius is \\( R h^{-1} Mpc\\). Then, the search radius in arcmin is

$$ \begin{aligned} &w\theta = \frac{c}{H_0}\theta\alpha(z_1, z_2) = \frac{c\times10^5 Km\cdot s^{-1}}{100 h Km\cdot s^{-1} {Mpc}^{-1}}\theta \alpha(z_1,z_2) = R h^{-1} Mpc \\\\ &\Rightarrow \theta = \frac{R}{1000c\alpha(z_1,z_2)}\frac{180\times60}{\pi} = \frac{10.8R}{c\pi\alpha(z_1, z_2)}\end{aligned}$$

The \\(\alpha(z_1, z_2)\\) is typically \\(10^{-1} ~ 10 ^{2}\\).







# <center>Process

1). Run Prepare_data.py collect data from each field. The data will be selected by some cutoffs. The name of the result file is "cata_result_ext.hdf5".

"mpirun -np ....  prepare_data.py collect"

2). Run the sym_mc_plot_cfht.py for cutoffs. Then determine the cutoff threshold (flux_alt or ..) according to the results (multiplicative bias and additive bias).

3). Run prepare_data.py to select the data needed. The name of the result file is "cata_result_ext_cut.hdf5". 

"mpirun -np ....  prepare_data.py select"

4). Run the C++ program to build the grid and assign the source to each grid for final calculation.

"mpirun -n ....  "
