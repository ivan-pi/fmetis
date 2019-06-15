---
project: Fortran METIS Interface
src_dir: ./src
output_dir: ./docs
project_github: https://github.com/ivan-pi/fmetis
author: Ivan Pribec
summary: Fortran METIS Interface
github: https://github.com/ivan-pi
email: ivan.pribec@gmail.com
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
source: true
graph: true
preprocess: false
exclude_dir: ./src/tests
exclude: parmetis_interface.f90
---

# Brief description

This is a Fortran interface to the [METIS software package](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) 
for partitioning unstructured graphs, partitioning meshes, and computing fill-reducing orderings
of sparse matrices. The interface makes use of the C interoperability features available in modern Fortran 
(i.e., Fortran 2003+) and provides a simple and safe way to call the original serial routines.

# License



# Further information

* [METIS Home page](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
* [METIS Manual (PDF)](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf)