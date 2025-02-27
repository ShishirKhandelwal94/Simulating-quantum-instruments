%% For a specific noise model (see next section for worst-case noise)

clc
clearvars

N=2; %Specify number of outcomes
d=2; %Specify dimension
plist=genpart(N,d); %Generates a list of all unique rank vectors, see genpart function
Nr=size(plist,1); %Gives the total number of unique rank vectors

%Instrument Choi matrices need to be specified as a list eta(:,:,a), where a is an outcome
%Example: below we give the Choi matrices for weak Z in three dimensions

gamma = 0.9; %Degree of sharpness

for i=1:4
    eta(:,:,i) = eta_a(i,d,gamma); %Defines weak Z Chois using the eta_a function
end

% Select type type of noise

%noise = eye(d^2)/(d^2*N); %White noise
noise = dephasing_noise(d,N);  %Dephasing noise, see dephasing_noise function



cvx_begin sdp
cvx_solver sedumi
cvx_precision best

    variable X(d^2,d^2,Nr,N) hermitian semidefinite %For each element of all rank vectors, we have one X
    %corresponding to a subnormalised Choi for the corresponding element
    variable v %visibility

%setting all zero rank Xs to zero and applying reduction constraint for Schmidt number for all non-zero Xs
%with rank less than the dimension
for i=1:Nr
    for j=1:N
        if plist(i,j)==0
            X(:,:,i,j)==0;
        elseif plist(i,j)>0 && plist(i,j)<=d-1
             kron(eye(d),PartialTrace(X(:,:,i,j),1)) - X(:,:,i,j)/plist(i,j) >= 0; %Schmidt number constraint
            %PartialTranspose(X(:,:,i,j),1)>=0; %If d=2, use this PPT contraint instead
        end
    end
end

%Choi constraint on total Choi states for each rank vector and trace
%constraint on each X
for i=1:Nr
    sumX=0;
    for j=1:N
        sumX = sumX + X(:,:,i,j);
    end
    PartialTrace(sumX,1) == trace(sumX)*eye(d)/d;
    for k=1:N
        trace(X(:,:,i,k)) == trace(sumX)*plist(i,k)/d;
    end
end

        
%Main decomposition constraint

for a=1:N
    sum=0;
    for j=1:Nr
        sum = sum+X(:,:,j,a);
    end 
    sum == v*eta(:,:,a)+(1-v)*noise;
end

%objective
maximize(v)

         
cvx_end



%% For worst-case noise. Here noise is a variable, and trace constraints need to be applied to it
clc
clearvars

N=2; %Specify number of outcomes
d=2; %Specify dimension
plist=genpart(N,d); %Generates a list of all unique rank vectors, see genpart function
Nr=size(plist,1); %Gives the total number of unique rank vectors

%Instrument Choi matrices need to be specified as a list eta(:,:,a), where a is an outcome
%Example: below we give the Choi matrices for weak Z in three dimensions

gamma = 0.1; %Degree of sharpness



for a=1:N
    eta(:,:,a) = eta_a(a,d,gamma);%Defines weak Z Chois using eta_a function
end



cvx_begin sdp
cvx_solver sedumi
cvx_precision best

    variable X(d^2,d^2,Nr,N) hermitian semidefinite %For each element of all rank vectors, we have one X
    %corresponding to a subnormalised Choi for the corresponding element
    variable noise(d^2,d^2,N) hermitian semidefinite %arbitrary noise
    variable v %visibility

%setting all zero rank Xs equal to and applying reduction constraint for Schmidt number for all non-zero Xs
%with rank less than the dimension
for i=1:Nr
    for j=1:N
        if plist(i,j)==0
            X(:,:,i,j)==0;
        elseif plist(i,j)>0 && plist(i,j)<=d-1
             kron(eye(d),PartialTrace(X(:,:,i,j),1)) - X(:,:,i,j)/plist(i,j) >= 0; %Schmidt number constraint
            %PartialTranspose(X(:,:,i,j),1)>=0; %If d=2, use this PPT contraint instead
        end
    end
end

%Choi constraint on total Choi states for each rank vector and trace
%constraint on each X
for i=1:Nr
    sumX=0;
    for j=1:N
        sumX = sumX + X(:,:,i,j);
    end
    PartialTrace(sumX,1) == trace(sumX)*eye(d)/d;
    for k=1:N
        trace(X(:,:,i,k)) == trace(sumX)*plist(i,k)/d;
    end
end

%Noise constraints

snoise=0;
for a=1:N
    snoise=snoise+noise(:,:,a);
end
PartialTrace(snoise,1)==(1-v)*eye(d)/d;

        
%Main decomposition constraint

for a=1:N
    sum=0;
    for j=1:Nr
        sum = sum+X(:,:,j,a);
    end 
    sum == v*eta(:,:,a)+noise(:,:,a);
end

%objective
maximize(v)
        

cvx_end


