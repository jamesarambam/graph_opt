# Graph_Opt

This repository implementations paper, Graph Based Optimization For Multiagent Cooperation (AAMAS-2019).



# Dependencies

1. Python 2.7
2. GNU Scientific Library (GSL) - 2.3



# GSL installation on Ubuntu 16.04

Create a folder name "gsl" in your home directory

```
$mkdir gsl
```

Download GSL-2.3 package from ftp://ftp.gnu.org/gnu/gsl

```
$tar -zxvf gsl-2.3.tar.gz
$cd gsl-2.3/
$./configure --prefix=/home/<username>/gsl
$make
$make check
$make install
```

Add line  `` export LD_LIBRARY_PATH="/home/<username>/gsl/lib:$LD_LIBRARY_PATH"`` in .bashrc, file(located at /home/username/.bashrc)

 

