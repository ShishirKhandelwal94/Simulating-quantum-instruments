%% For a specific noise model (see next section for worst-case noise)

clc
clearvars

N=3; %Specify number of outcomes
d=3; %Specify dimension
plist=genpart(N,d); %Generates a list of all unique rank vectors, see genpart function
Nr=size(plist,1); %Gives the total number of unique rank vectors

%Instrument Choi matrices need to be specified as a list eta(:,:,a), where a is an outcome
%Example: below we give the Choi matrices for weak Z in three dimensions

gamma = 0.5; %Degree of sharpness
eta(:,:,1) = [
    (1 + gamma) / 6, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12
];

eta(:,:,2) = [
    (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 + gamma) / 6, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12
];

eta(:,:,3) = [
    (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 + gamma) / 6
];

% Select type type of noise

noise = eye(d^2)/(d^2*N); %White noise
%noise = dephasing_noise(d,N);  %Dephasing noise, see dephasing_noise function



cvx_begin sdp

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

%Choi constraint on total Choi states for each rank vector
for i=1:Nr
    sumX=0;
    for j=1:N
        sumX = sumX + X(:,:,i,j);
    end
    PartialTrace(sumX,1) == trace(sumX)*eye(d)/d;
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

N=3; %Specify number of outcomes
d=3; %Specify dimension
plist=genpart(N,d); %Generates a list of all unique rank vectors, see genpart function
Nr=size(plist,1); %Gives the total number of unique rank vectors

%Instrument Choi matrices need to be specified as a list eta(:,:,a), where a is an outcome
%Example: below we give the Choi matrices for weak Z in three dimensions

gamma = 0.5; %Degree of sharpness
eta(:,:,1) = [
    (1 + gamma) / 6, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12
];

eta(:,:,2) = [
    (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 + gamma) / 6, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 - gamma) / 12
];

eta(:,:,3) = [
    (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    (1 - gamma) / 12, 0, 0, 0, (1 - gamma) / 12, 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2));
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, sqrt(1 - gamma^2) / (6 * sqrt(2)), 0, 0, 0, (1 + gamma) / 6
];






cvx_begin sdp

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

%Choi constraint on total Choi states for each rank vector
for i=1:Nr
    sumX=0;
    for j=1:N
        sumX = sumX + X(:,:,i,j);
    end
    PartialTrace(sumX,1) == trace(sumX)*eye(d)/d;
end

%Noise constraints
for a=1:N
    PartialTrace(noise(:,:,a),1)==trace(noise(:,:,a))*eye(d)/d;
end
snoise=0;
for a=1:N
    snoise=snoise+noise(:,:,a);
end
trace(snoise)==1-v;
        
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
