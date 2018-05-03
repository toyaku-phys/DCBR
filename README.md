# Definition of the function

The displacement correlation function for pair of residues is defined as

<img src="https://latex.codecogs.com/gif.latex?C(i,j):=\langle\vec{u}_i(t)\cdot\vec{u}_j(t)\rangle_t." />

Here,
<img src="https://latex.codecogs.com/gif.latex?\vec{u}(t):=p(t+\Delta)-p(t)." />

Our DCBR program provides 
<img src="https://latex.codecogs.com/gif.latex?\{C(i,j)\}" />
, where residue index <img src="https://latex.codecogs.com/gif.latex?i" /> is fixed.

## Using DCBR
dcbr in=input_file.pdb out=output_file_name step=0-1000 target=index_of_reference_residue
