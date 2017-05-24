# CME 342 - Homework 2

## Installation

### Prerequisites -- METIS and PARMETIS
First, install both the metis and parmetis libraries:

```bash
$ wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
$ wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
$ tar xzf metis-5.1.0.tar.gz
$ tar xzf parmetis-4.0.3.tar.gz
$ cd metis-5.1.0
$ make config
$ make
$ cd ../parmetis-4.0.3
$ make config
$ make
```

### Clone repository and build this code

```bash
$ git clone https://github.com/rehnd/cme342.git
$ cd cme342/hw2
```

Now compile the serial code

```bash
$ cd 1-serial
$ make
```

and next, the parallel

```bash
$ cd ~/cme342/hw2/2-parallel
$ make
```

### Executing the code after compilation

The executables require at least three additional command-line options.
An example for the `oneram6` mesh file is the following:

```bash
$ ./serial 2 ../meshes/oneram6.conn ../meshes/oneram6.xyz
```

The correspoding parallel call is

```bash
$ mpirun -np 2 ./parallel 2 ../meshes/oneram6.conn ../meshes/oneram6.xyz
```

where we have run on 2 processors.
