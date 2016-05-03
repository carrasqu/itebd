clear;
rng(1);
chimax=25; % bond dimension
d=2; % Dimension of the local hilbert space
Jz=4.0; % Jz coupling
h=1; % transverse field
eps=0.005; % time step
tau=eps; % imaginary time evolution
terror=0.00000000001;

 
% randomly initializing tensors defining the wave function
LA=rand(chimax,d,chimax);
lambA=rand(chimax,chimax);
LB=rand(chimax,d,chimax);
lambB=rand(chimax,chimax);
 
% pauli matrices
x=0.5*[0 1;1 0];
y=-0.5*[0 1i;-1i 0]; 
z=0.5*[1 0; 0 -1];
e=[1 0; 0 1];


% s+ s-
sp=x+1i*y;
sm=x-1i*y;

% --- gates ----

% two-site operators defining the Hamiltonian gates to be applied duuring the time
% evolution
szsz=reshape(kron(z,z),2,2,2,2);
sx1=reshape(kron(x,e),2,2,2,2);
sx2=reshape(kron(e,x),2,2,2,2);
id2=reshape(kron(e,e),2,2,2,2);


hl=reshape(Jz*kron(z,z)+h*kron(x,e)+h*kron(e,x),2,2,2,2); % local hamiltonian
Jz*kron(z,z)+h*kron(x,e)+h*kron(e,x);

U=reshape(expm(-tau*reshape(hl,4,4)),2,2,2,2 );


time=200; % maximum time
tsteps=ceil(100/eps);


% imaginary time evolution using U gate
 for t=1:4000;
     
 [ LA,lambA,LB,lambB,U,ee ] = itebdA( LA,lambA,LB,lambB,U,chimax,d,terror,tau );    
    
 [ LB,lambB,LA,lambA,U,ee ] = itebdA( LB,lambB,LA,lambA,U,chimax,d,terror,tau ); 
 
%  [energy1]=tedben(LB,lambB,LA,lambA,hl,chimax);
%  [energy2]=tedben(LA,lambA,LB,lambB,hl,chimax);
%  energy=(energy1+energy2)/2.0
%  ee
 
 end
 
%energy measurements
[energy1]=tedben(LB,lambB,LA,lambA,hl,chimax);
[energy2]=tedben(LA,lambA,LB,lambB,hl,chimax);
 energy=(energy1+energy2)/2.0
 ee
 -1.27323954473516 % exact energy at the critical point in the thermodynamic limit which is the case I am considering 
 % at criticality e=-4/pi =             -1.27323954473516
 

return

