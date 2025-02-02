%% Projective instrument
%%Let outcome 0 correspond to N and 1 to S. Given a success probability,
%%this code gives the optimal fidelity using a projective instrument
%%Check next section for code for arbitrary instrument
clc 
clearvars

pmax = 3/4;
phip = [1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]/2; %maximally entangled state for qubits
pwin1=2/3; %Specify success probability

cvx_begin sdp
% Define CVX variables for decomposition
variable x0(4, 4) hermitian semidefinite
variable x1(4, 4) hermitian semidefinite
variable y0(4, 4) hermitian semidefinite
variable y1(4, 4) hermitian semidefinite
   
%Defining etas
eta0 = x0+y0;
eta1 = x1+y1;   

F = 1/3+ (2/3)*trace(phip*(eta0+eta1)); %Forumula for average fidelity


  
%probability constraint 
pwin = 1/2+(1/4)*trace(2*PartialTrace(eta0,1)*[1 0;0 -1]);%Formula for success probability
pwin == pwin1;



%Separability constraint on Chois corresponding to rank 1 projectors
PartialTranspose(x0) >= 0;
PartialTranspose(x1) >= 0;

% Trace constraint on Chois:
trace(x0+x1+y0+y1)==1;
PartialTrace(x0+x1,1)==trace(x0)*eye(2);
PartialTrace(y0,1)==trace(y0)*eye(2)/2;
PartialTrace(y1,1)==trace(y1)*eye(2)/2;

%Objectiv: given success probability, maximize F
maximize(F)
  
cvx_end

%% Arbitrary instrument
%%Let outcome 0 correspond to N and 1 to S. Given a success probability,
%%this code gives the optimal fidelity using any quantum instrument

phip = [1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]/2; %Maximally entangled state for qubits
pwin1=2/3;

cvx_begin sdp
cvx_solver sedumi
cvx_precision best

% Define CVX variables for decomposition
variable eta0(4, 4) hermitian semidefinite
variable eta1(4, 4) hermitian semidefinite

   

F = 1/3+ (2/3)*trace(phip*(eta0+eta1)); %Formula for average fidelity

% Trace constraint:
PartialTrace(eta0+eta1,1)==eye(2)/2;     

  
%probability constraint 
pwin = 1/2+(1/4)*trace(2*PartialTrace(eta0,1)*[1 0;0 -1]); %Formula for success probability
pwin==pwin1;

    maximize(F)
  
cvx_end

