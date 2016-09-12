%{
Reaction Diffusion solver (Forward time, center space)
Matt Bovyn
September 6 2016

Function to create the matrix which encodes the FTCS scheme. 

It takes in:
mu1 and mu2, the time step over the space step in each direction squared
n, the number of points on the space grid

It puts out:
a sparse matrix which, when multiplied by the last time point converted 
into a vector, generates the next time point.

%note: this creates full matricies first, then converts into sparse. This
is not the best way to do this because the full memory must be allocated
for the big matricies and this can run out quickly (can't do more than
about 100 points in space). To do this correctly, everything should be
sparse from the beginning. Should rewrite at some point if we want to do
higher spatial resolution.
%}

function mat=FTCS_matrix(n,mu1,mu2)
%Creates the matrix for use with the Forward Time Center Space scheme

%vectors to put on our diagonals
center=(1-2*mu1-2*mu2)*ones(n,1);
top=mu1*ones(n-1,1);
bottom=mu1*ones(n-1,1);

%create a small matrix that is our block
block=diag(center,0)+diag(top,1)+diag(bottom,-1);
block(1,n)=mu1;
block(n,1)=mu1;

%create a block diagonal matrix by repeating the block we've made n times
blocks_only=kron(eye(n),block);

%create a matrix with just the mu2's that wrap around the outside
col=zeros(n^2,1);
col(n+1)=mu2;
col(end-(n-1))=mu2;
row=zeros(n^2,1);
row(n+1)=mu2;
row(end-(n-1))=mu2;
wrappers_only=toeplitz(col,row);

%put them together into a sparse matrix
mat=sparse(wrappers_only+blocks_only);

end
