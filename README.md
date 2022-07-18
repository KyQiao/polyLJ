# polyLJ

polydisperse LJ potential for LAMMPS

LAMMPS doesn't provide a native implementation of the polydisperse system. One tedious way of doing that is defining hundreds of types to approach the desired size distribution, which also needs to set inter-type parameters.

This potential uses charge as atom size in LJ potential. The general form of poly-LJ is:

<img src="https://latex.codecogs.com/svg.image?\phi(r/\sigma_{\alpha\beta})=\left\{\begin{align*}&v_0[(\frac{\sigma_{\alpha\beta}}{r})^m-(\frac{\sigma_{\alpha\beta}}{r})^n&space;]&plus;\sum_{k=0}^qc_k(\frac{r^{\alpha\beta}}{\sigma_{\alpha\beta}})^{2k}&space;\quad&space;&r/\sigma_{\alpha\beta}\le&space;\tilde&space;r_c&space;\\&space;&0&space;\quad&space;&otherwise\end{align*}&space;\right&space;\}" title="\phi(r/\sigma_{\alpha\beta})=\left\{\begin{align*}&v_0[(\frac{\sigma_{\alpha\beta}}{r})^m-(\frac{\sigma_{\alpha\beta}}{r})^n ]+\sum_{k=0}^qc_k(\frac{r^{\alpha\beta}}{\sigma_{\alpha\beta}})^{2k} \quad &r/\sigma_{\alpha\beta}\le \tilde r_c \\ &0 \quad &otherwise\end{align*} \right \}" />

and 

<img src="https://latex.codecogs.com/svg.image?\sigma_{\alpha\beta}=\frac{1}{2}(\sigma_\alpha&plus;\sigma_\beta)(1-\epsilon|\sigma_\alpha-\sigma_\beta|)" title="\sigma_{\alpha\beta}=\frac{1}{2}(\sigma_\alpha+\sigma_\beta)(1-\epsilon|\sigma_\alpha-\sigma_\beta|)" />

the polynomial in the end is one way of flattening the cutoff. 

The current version doesn’t introduce the polynomial term, only the normal cutoff version is implemented.

## Installation

Move all c++ files (sources and headers) to LAMMPS source folder (e.g. `/lammps/src`) and compile again.

## Notes on using the potential:

- How to set neighbor correct?

LAMMPS uses cutoff to build the neighbor list. Here we use an overall cutoff for all sizes, which means the skin distance must be long enough to cover the largest particle. For example, when the cutoff is 2.5 <img src="https://latex.codecogs.com/svg.image?\sigma_{ij}" title="\sigma_{ij}" /> and the largest size is 1.55 <img src="https://latex.codecogs.com/svg.image?\sigma" title="\sigma" />:


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

A smallest case to use `swap` is:

```
units           lj
#must be the charge style
atom_style	    charge
dimension       2 
boundary        p p p

#to import size from outside
atom_modify     map yes

lattice         sq 0.3 spacing 0.5 0.5 0.5
region          box block 0 32  0 32 -0.5 0.5  units box
create_box      1 box
create_atoms    1 region box units box
mass            1 1.0


# python generate.py 324 --shuffle --addfig --binary
variable        size atomfile 324atomShuffle.xyz

set             atom * charge v_size


# set global cut off (rc) to 2.5
pair_style      lj/cut/poly 2.5
# sigma here are dumb variable
# args are type type  e non-addtive-para
# set non-addtive para to 0.2 provent xtal
pair_coeff      1 1  1 0.2


velocity        all create 0.01 87287
velocity        all zero linear
velocity        all zero angular


neighbor        1.4 bin 
#1.375 is minimum skin distance
neigh_modify    every 20 delay 0 check no 
#check no to make sure build neighbor list correctly 


fix             1 all nvt temp 0.01 0.05 0.25
# 	swap 10000 times every 40000 steps with a random seed 111821 with accept Temp threshold 0.1
# 	last parameter means start swap when T<0.1, 
fix             swap all swap 100 400 111821 0.1
fix             2d all enforce2d

thermo          1000


run             10000
unfix           2d
unfix           swap
```

