# Definition of the correlation function

The displacement correlation function for pair of residues is defined as

<img src="https://latex.codecogs.com/gif.latex?C(i,j):=\langle\vec{u}_i(t)\cdot\vec{u}_j(t)\rangle_t" />. with <img src="https://latex.codecogs.com/gif.latex?\vec{u}_i(t):=\vec{p}_i(t+\Delta)-\vec{p}_i(t)." />

Here, <img src="https://latex.codecogs.com/gif.latex?\vec{p}_i" /> is center of atom positions of 
<img src="https://latex.codecogs.com/gif.latex?i" />
-th residue (without "CA" atom).
<img src="https://latex.codecogs.com/gif.latex?\Delta" /> is time interval in the input file.
Time average <img src="https://latex.codecogs.com/gif.latex?C(i,j):=\langle\cdots\rangle_t" /> 
is calculated using all samples.
If you need the center of mass or diffrent time interval, you can *write* next DCBR!


Our DCBR program provides 
<img src="https://latex.codecogs.com/gif.latex?\{C(i,j)\}" />
, where index of residue <img src="https://latex.codecogs.com/gif.latex?i" /> is fixed (this residue is reference).

## Using DCBR
dcbr in=input_file.pdb out=output_file_name step=0-1000 target=index_of_reference_residue

## Building
example: clang++-5.0 -std=c++1z -stdlib=libc++ main.cpp -O2 -o dcbr
