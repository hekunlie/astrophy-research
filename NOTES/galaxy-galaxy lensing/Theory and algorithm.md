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
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>

The critial surface density in comoving distance:

<table border="0">
 <tr>
    <td>$$\Sigma_{crit}(z_L, z_S)$$</b></td>
    <td>$$\Sigma_{crit}(z_L, z_S)$$</b></td>
 </tr>
 <tr>
   <td>
   $$ \begin{aligned}&= \frac{c^2}{4\pi G}\frac{D_A (z_S)}{D_A (z_L) D_A (z_L, z_S) (1+z_L)^2}, z_L < z_S \\\\
   &= \frac{c^2}{4 \pi G} \frac{\omega(0,z_S)}{\omega(0,z_L)[\omega(0,z_S) - \omega(0,z_L)] (1+z_L)} \\\\
   &= L(z_L, z_s)\frac{c_0^2\times 10^{16}\ \rm{m^2\cdot s^{-2} \cdot h \cdot Mpc^{-1}}}{4\pi G_0\times 10^{-11}\  \rm{m^3\cdot s^{-2} \cdot Kg^{-1}}} \\\\
   & = L(z_L, z_s)\frac{c_0^2}{4\pi G_0}\times 10^{21}\ \rm{h \cdot Kg \cdot m^{-1} \cdot pc^{-1}} \\\\
   & = L(z_L, z_s)\frac{c_0^2}{4\pi G_0 l_{pc} m_{\odot}}\times 10^{7}\ \rm{h \cdot M_{\odot} \cdot pc^{-2}} \\\\
   & = L(z_L, z_s)\times 174648.0703379\ \rm{h \cdot M_{\odot} \cdot pc^{-2}}, \\\\
   &L(z_L, z_s) = \frac{\omega(0,z_S)}{\omega(0,z_L)[\omega(0,z_S) - \omega(0,z_L)] (1+z_L)}
   \end{aligned}$$
   </td>
   <td>
   $$ \begin{aligned} &= \frac{c^2}{4\pi G}\frac{D_A (z_S)}{D_A (z_L) D_A (z_L, z_S) (1+z_L)^2}, z_L < z_S \\\\
   &= \frac{c^2}{4 \pi G}\frac{H_0}{c} \frac{\alpha(0,z_S)}{\alpha(0,z_L)[\alpha(0,z_S) - \alpha(0,z_L)] (1+z_L)} \\\\
   &= L^{\prime}(z_L, z_S)\frac{c_0\times 10^{16}\ \rm{m^2\cdot s^{-2}} \times 10^{-3} \ \rm{h \cdot Mpc^{-1}}}{4\pi G_0\times 10^{-11}\  \rm{m^3\cdot s^{-2} \cdot Kg^{-1}}}  \\\\
   &= L^{\prime}(z_L, z_S)\frac{c_0}{4\pi G_0}\times 10^{24} \ \rm{h \cdot Kg \cdot m^{-1} \cdot Mpc^{-1}} \\\\
   &= L^{\prime}(z_L, z_S)\frac{c_0 l_{pc}}{4\pi G_0 m_{\odot}}\times 10^{4} \ \rm{h \cdot M_{\odot} \cdot pc^{-2}} \\\\
   & = L^{\prime}(z_L, z_S)\times 554.682\ 135\ 528\ \rm{h \cdot M_{\odot} \cdot pc^{-2}} \\\\
   &L^{\prime}(z_L, z_S) = \frac{\alpha(0,z_S)}{\alpha(0,z_L)[\alpha(0,z_S) - \alpha(0,z_L)] (1+z_L)}
   \end{aligned}$$
   </td>
 </tr>
</table>


$$ c_0 = 2.99\ 792\ 458; \quad  l_{pc} = 3.085\ 677\ 581; \quad m_{\odot} = 1.9885, \quad M_{\odot} = 1.9885\times 10^{30}\  \rm{Kg}; \quad G_0 = 6.67\ 408\ 313$$


The rotation of the shear eastimators: \\(G_1\\), \\(G_2\\), \\(N\\), \\(U\\), \\(V\\).

<table border="0">
 <tr>
    <td>Spin-0</b></td>
    <td>Spin-2</b></td>
    <td>Spin-4</b></td>
 </tr>
 <tr>
 <td>$$N \rightarrow N$$</td>
    <td>$$\begin{aligned} G_1^{\prime}+iG_2^{\prime} &= (G_1+iG_2)\exp(2i\phi) \\\\
   G_1^{\prime} &= G_1\cos(2\phi)  - G_2\sin(2\phi) \\\\
   G_2^{\prime} &= G_1\sin(2\phi)  + G_2\cos(2\phi) \end{aligned}$$</td>
    <td>$$\begin{aligned} U^{\prime}+iV^{\prime} &= (U+iV)\exp(4i\phi) \\\\
   U^{\prime} &= U\cos(4\phi)  - V\sin(4\phi) \\\\
   V^{\prime} &= U\sin(4\phi)  + V\cos(4\phi) \end{aligned}$$</td>
 </tr>
</table>

# <center>NFW:
Profile:

 $$ \rho(r)= \frac{\rho_s}{[\frac{r}{r_s}][1+\frac{r}{r_s}]^2}$$

Mass:

 $$ \begin{aligned} M &= \int_0^{2\pi}d\phi\int_0^{\pi}d\theta\ \sin\theta \int_0^{r_{\Delta}}dr\ r^2\frac{\rho_s}{ [ \frac{r}{r_s} ][ 1+\frac{r}{r_s} ]^2 } \\\\ &= 4 \pi \rho_s r_s^3 \left[ \ln (1 + \frac{r}{r_s}) + \frac{1}{1+\frac{r}{r_s}} \right], r=0, r_{\Delta} \\\\ &= 4\pi\rho_s r_s^3 \left[\ln(1+\frac{r_\Delta}{r_s}) - \frac{r_\Delta}{r_s +r_\Delta}\right] \end{aligned}$$

Projected surface density:

$$ \begin{aligned} \Sigma(\theta) &= \int_{\alpha_1(\theta)}^{\alpha_2(\theta)}d\alpha\ \sin\alpha \int_{r_m}^{r_{\Delta}}dr \frac{\rho_s}{[\frac{r}{r_s}][1+\frac{r}{r_s}]^2} \\\\ &= [\cos\alpha_1(\theta) - \cos\alpha_2(\theta)]\ \rho_s r_s^3 \left[\ln(1+\frac{r}{r_s}) + \frac{1}{1+\frac{r}{r_s}}\right], r=r_m, r_{\Delta} \\\\
& = [\cos\alpha_1(\theta) - \cos\alpha_2(\theta)]\ \rho_s r_s^3 \left[\ln(\frac{r_s + r_{\Delta}}{r_s + r_m}) + \frac{r_s}{r_s + r_{\Delta}}  - \frac{r_s}{r_s + r_m}\right]\end{aligned}$$

<br>
<br>

<table border="0">
 <tr>
    <td>g</b></td>
    <td>G</b></td>
    <td>U & V</b></td>
 </tr>
 <tr>
    <td>$$ \begin{aligned} g &=(g_1+ig_2)e^{2i\phi} \\\\
    &= g_1\cos(2\phi) - g_2\sin(2\phi) \\\\
    &+ i( g_1\sin(2\phi) + g_2 \cos(2\phi)) \\\\ 
    g_t &= -Re[g e^{-2i\phi}] \\\\ &= -g_1\cos(2\phi) - g_2\sin(2\phi) \\\\ 
    g_{\times} & = -Im[g e^{-2i\phi}] \\\\
    & = g_1\sin(2\phi) - g_2 \cos(2\phi)\end{aligned}$$</td>
    <td>Rotate as above $$\begin{aligned} G &= (G_1+iG_2) e^{2i\phi} \\\\
    G_t &= G_1\cos(2\phi)  - G_2\sin(2\phi) \\\\
    G_{\times} &= G_1\sin(2\phi) + G_2\cos(2\phi) \end{aligned}$$</td>
    <td>Rotate as above $$\begin{aligned} U^{\prime} &= (U+iV) e^{4i\phi} \\\\
    U_t &= U\cos(4\phi) - V\sin(4\phi) \\\\ \end{aligned}$$</td>
 </tr>
</table>
<img src="rotation.jpg" width = 100% height = 100% div align=left/>

# <center>Process

1). Run Prepare_data.py collect data from each field. The data will be selected by some cutoffs. The name of the result file is "cata_result_ext.hdf5".

"mpirun -np ....  prepare_data.py collect"

2). Run the sym_mc_plot_cfht.py for cutoffs. Then determine the cutoff threshold (flux_alt or ..) according to the results (multiplicative bias and additive bias).

3). Run prepare_data.py to select the data needed. The name of the result file is "cata_result_ext_cut.hdf5". 

"mpirun -np 4 python prepare_data.py select"

4). Run the ggl_com_dist.cpp to assign the comoving distance (only the integrate part) to each galaxy in the file "cata_result_ext_cut.hdf5". The comoving distances have been calculated for 0 to 10 with an interval \\(\delta z = 0.0001\\).

"./ggl_com_dist"

5). Run ggl_grid.cpp to build the grid for background galaxies.

"mpirun -n 30 ./ggl_grid 0.15" ...

6). Run prepare_foreground.py to prepare the foreground data for measurement.

"python prepare_foreground.py"

7). Run the C++ program to build the grid and assign the source to each grid for final calculation.



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


