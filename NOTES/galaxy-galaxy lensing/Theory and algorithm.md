<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>


Hubbule parameter:

$$ \begin{aligned}H^2 &= H_0^2[\frac{\Omega_r}{a^4}+\frac{\Omega_m}{a^3}-\frac{Kc^2}{a^2H_0^2}+\Omega_{\Lambda}] \\\\
&=H_0^2[\frac{\Omega_r}{a^4}+\frac{\Omega_m}{a^3}+\frac{1-\Omega_0}{a^2}+\Omega_{\Lambda}] \end{aligned}$$

Comving distance:

$$ dt = \frac{da}{\dot{a}} \Rightarrow -dw = \frac{cdt}{a} = \frac{cda}{a\dot{a}}=\frac{cda}{a^2H}$$

$$ \begin{aligned}\omega(z_1,z_2) &= \frac{c}{H_0}\int_{a(z_2)}^{a(z_1)} \frac{da}{\sqrt{a\Omega_m + a^2(1-\Omega_m -\Omega_\Lambda) + a^4\Omega_\Lambda}}, z_1 < z_2 \\\\ 
&= \frac{c}{H_0}\int_{z_1}^{z_2} \frac{dz}{\sqrt{(1+z)^3\Omega_m + (1+z)^2(1-\Omega_m -\Omega_\Lambda) + \Omega_\Lambda}}, z_1 < z_2 \\\\
&= \frac{c}{H_0}\alpha(z_1, z_2) \end{aligned}$$

The physical distance relates the angle to the size of a distant object is the angular diameter distance \\(D_A\\). 

$$ D_A(z) = a(z)f_k(\omega(0,z))$$

The search radius is \\( R h^{-1} \rm{Mpc}\\). In the flat universe, the search radius in unit of degree is

$$ \begin{aligned} \omega(0,z)\theta &= \frac{c}{H_0}\alpha(0, z)\theta = \theta\alpha(0,z) \frac{c\times10^5 \rm{Km\cdot s^{-1}}}{100 h \rm{Km \cdot s^{-1}\cdot {Mpc}^{-1}}}  = R h^{-1} \rm{Mpc} \\\\ &\Rightarrow \theta = \frac{R}{1000c\alpha(0,z)}\frac{180}{\pi} = \frac{0.18}{c\pi}\frac{R}{\alpha(0, z)}\end{aligned} $$

The angular distance between \\(z_1\\) and \\(z_2\\) is

$$D_A(z_1, z_2) = a(z_2)f_k(\omega(z_1,z_2)) = \frac{c}{H_0}a(z_2)[\alpha(0, z_2) - \alpha(0, z_1)]$$

It is valid only for \\(\omega_k \geq 0 \\) (see astro-ph/9905116 or Principles of Physical Cosmology pp 336–337, Peebles)

The critial surface density in comoving distance:

$$ \begin{aligned}\Sigma_{crit}(z_1, z_2) &= \frac{c^2}{4\pi G}\frac{D_A (z_2)}{D_A (z_1) D_A (z_1, z_2) (1+z_1)^2}, z_1<z_2 \\\\
&= \frac{c^2}{4 \pi G} \frac{\omega(0,z_2)}{\omega(0,z_1)[\omega(0,z_2) - \omega(0,z_1)] (1+z_1)} \\\\
&=\frac{c^2}{4\pi G}\frac{H_0}{c}\frac{\alpha(0,z_2)}{\alpha(0,z_1)[\alpha(0,z_2) - \alpha(0,z_1)] (1+z_1)} \\\\
&= \frac{c\times 10^{8} \rm{m\cdot s^{-1}}100h\rm{Km\cdot s^{-1}{Mpc}^{-1}}}{4\pi G\times 10^{-11} \rm{m^3\cdot s^{-2} \cdot Kg^{-1}}}\frac{\alpha(0,z_2)}{[\alpha(0,z_2) - \alpha(0,z_1)]\alpha(0,z_1)(1+z_1)} \\\\
& = \frac{\alpha(0,z_2)}{[\alpha(0,z_2) - \alpha(0,z_1)]\alpha(0,z_1)(1+z_1)}\frac{c \cdot h}{4\pi G}\times 10^{24}\rm{Kg \cdot m^{-1}\cdot Mpc^{-1}} \\\\ 
& = \frac{\alpha(0,z_2)}{[\alpha(0,z_2) - \alpha(0,z_1)]\alpha(0,z_1)(1+z_1)}\frac{c \cdot h \cdot m_{pc} }{4\pi G\cdot m_s}\times 10^{4}\rm{M_{sum} \cdot pc^{-2}} \\\\
& = \frac{\alpha(0,z_2)}{[\alpha(0,z_2) - \alpha(0,z_1)]\alpha(0,z_1)(1+z_1)}3.88283351\times 10^{2}\rm{M_{sum} \cdot pc^{-2}}
\end{aligned}$$

$$  h = 0.7, \quad c = 2.99792458 \quad  m_{pc} = 3.085677581 \quad m_s = 1.98847 \quad G = 6.6740831313$$

Tangential shear and cross shear:

$$ \begin{aligned} g e^{-2i\phi} &=(g_1+ig_2)e^{-2i\phi} \\\\
&= g_1\cos(2\phi)+g_2\sin(2\phi) + i(g_2 \cos(2\phi) - g_1\sin(2\phi)) \\\\ 
g_t &= -Re[g e^{-2i\phi}] \\\\ &= -g_1\cos(2\phi) - g_2\sin(2\phi) \\\\ 
g_{\times} & = -Im[g e^{-2i\phi}] \\\\
&=g_1\sin(2\phi) - g_2 \cos(2\phi)\end{aligned}$$



# <center>Process

1). Run Prepare_data.py collect data from each field. The data will be selected by some cutoffs. The name of the result file is "cata_result_ext.hdf5".

"mpirun -np ....  prepare_data.py collect"

2). Run the sym_mc_plot_cfht.py for cutoffs. Then determine the cutoff threshold (flux_alt or ..) according to the results (multiplicative bias and additive bias).

3). Run prepare_data.py to select the data needed. The name of the result file is "cata_result_ext_cut.hdf5". 

"mpirun -np ....  prepare_data.py select"

4). Run the C++ program to build the grid and assign the source to each grid for final calculation.

"mpirun -n ....  "

# <center> Code structure

Loop the radius bin (in unit of degree):


1). Loop the all foreground galaxy to find the background source galaxies in radius bin ([radius_s, radius_e]) and label them.
~~~
The source galaxy: Z > Z(len) + 0.3 (will change with the new redshift catalog).

The searching radius is calculated by the former formula.

The seperation radius (in unit of degree) is calculated in the Cartesian coordinate because of the small seperation.

Each thread label the found source galaxy in the mask array. It is possible that one galaxy may be identified as the source galaxy of the two or more foreground galaxy due to the small seperation of foreground galaxies. Once it is identified as source galaxy, the corresponding mask will increase by 1.

Then they will assign the mask to the final mask that can seen by each thread.
~~~

2). The rank 0 thread calculates the shear in this radius bin,  [radius_s, radius_e].

~~~
Each source galaxy is just used once even though it may be identified as the source galaxy by many foregroudn galaxies.
~~~

1). Loop the foreground galaxy to calculate the tangential shear and the \\(\Delta \Sigma(R)\\)