# polyLJ

polydisperse LJ potential for LAMMPS

LAMMPS doesn't provide native implementation of polydisperse system. One tedious way of doing that is defining hundreds of type to approaching the desired size distribution, which also need to set inter-type parameter.

This potential use charge as atom size in LJ potential. General form of poly-lj is:

<img src="https://latex.codecogs.com/svg.image?\phi(r/\sigma_{\alpha\beta})=\left\{\begin{align*}&v_0[(\frac{\sigma_{\alpha\beta}}{r})^m-(\frac{\sigma_{\alpha\beta}}{r})^n&space;]&plus;\sum_{k=0}^qc_k(\frac{r^{\alpha\beta}}{\sigma_{\alpha\beta}})^{2k}&space;\quad&space;&r/\sigma_{\alpha\beta}\le&space;\tilde&space;r_c&space;\\&space;&0&space;\quad&space;&otherwise\end{align*}&space;\right&space;\}" title="\phi(r/\sigma_{\alpha\beta})=\left\{\begin{align*}&v_0[(\frac{\sigma_{\alpha\beta}}{r})^m-(\frac{\sigma_{\alpha\beta}}{r})^n ]+\sum_{k=0}^qc_k(\frac{r^{\alpha\beta}}{\sigma_{\alpha\beta}})^{2k} \quad &r/\sigma_{\alpha\beta}\le \tilde r_c \\ &0 \quad &otherwise\end{align*} \right \}" />

and 

<img src="https://latex.codecogs.com/svg.image?\sigma_{\alpha\beta}=\frac{1}{2}(\sigma_\alpha&plus;\sigma_\beta)(1-\epsilon|\sigma_\alpha-\sigma_\beta|)" title="\sigma_{\alpha\beta}=\frac{1}{2}(\sigma_\alpha+\sigma_\beta)(1-\epsilon|\sigma_\alpha-\sigma_\beta|)" />

the polynomial in the end is one way of flatten the cutoff. 

The current version don't introduce the polynomial term, only normal cutoff version is implemented.

## Notes on using the potential:

- How to set neighbor correct?

LAMMPS uses cutoff to build neighbor list. Here we use a overall cutoff for all size, means the skin distance must be long enough to cover the largest particle. For example, when cutoff is 2.5 <img src="https://latex.codecogs.com/svg.image?\sigma_{ij}" title="\sigma_{ij}" /> and largest size is 1.55 <img src="https://latex.codecogs.com/svg.image?\sigma" title="\sigma" />:


```
neighbor        1.5 bin 
#1.375 is minimum skin distance
neigh_modify    every 20 delay 0 check no 
#check no to make sure build neighbor list correctly 
```

- How to set pair coefficient?

```
#global cutoff
pair_style      lj/cut/poly 2.5
# args are type type  e non-addtive-para
pair_coeff      1 1  1 0.2 
```

The non-additive parameter is a global variable. It will be set when first occurs.

- How to set size distribution?

Use `generate.py`  to generate a atomfile for LAMMPS. The file contains id and size of particle

```
atom_style	    charge
variable        size atomfile number_of_atom.xyz
set             atom * charge v_size
```

use `python generata.py -h` to get help. Now only a bimodal size distribution function is implemented.

