% This code performs a seesaw search in the trade-off between two 
% sequential CHSH parameter.  It returns an explicit quantum model for the
% double violation.
% Input: alpha (S_AB >= alpha)
% Needs: CVX and qetlab


clear all
clc



X=[0 1;1 0];Y=[0 -i;i 0]; Z=[1 0;0 -1]; id=eye(2);



dat=[];
T=[1 1; 1 -1];

for alpha=[2.13:0.01:2.13];

    ll=[];
for repeat=1:50
    
for x=1:2
        sigma(:,:,1,x)=RandomDensityMatrix(2,1,1)/2;
        sigma(:,:,2,x)=id/2-sigma(:,:,1,x);
end


M(:,:,1)=X;%RandomDensityMatrix(2,0,1);%(X+Z)/sqrt(2);
M(:,:,2)=Z;%RandomDensityMatrix(2,0,1);%(X-Z)/sqrt(2);






for K=1:100


    
    clear eta rhopost
cvx_begin sdp quiet
cvx_solver mosek
variable x1(4,4,2) hermitian semidefinite
variable y1(4,4,2) hermitian semidefinite
variable x2(4,4,2) hermitian semidefinite
variable y2(4,4,2) hermitian semidefinite


eta(:,:,1,1)=x1(:,:,1)+y1(:,:,1);
eta(:,:,2,1)=x1(:,:,2)+y1(:,:,2);
eta(:,:,1,2)=x2(:,:,1)+y2(:,:,1);
eta(:,:,2,2)=x2(:,:,2)+y2(:,:,2);


PartialTranspose(x1(:,:,1))>=0;
PartialTranspose(x1(:,:,2))>=0;
PartialTranspose(x2(:,:,1))>=0;
PartialTranspose(x2(:,:,2))>=0;


PartialTrace(x1(:,:,1)+x1(:,:,2),1)==trace(x1(:,:,1)+x1(:,:,2))*eye(2)/2;
PartialTrace(x2(:,:,1)+x2(:,:,2),1)==trace(x2(:,:,1)+x2(:,:,2))*eye(2)/2;
PartialTrace(y1(:,:,1),1)==trace(y1(:,:,1))*eye(2)/2;
PartialTrace(y1(:,:,2),1)==trace(y1(:,:,2))*eye(2)/2;
PartialTrace(y2(:,:,1),1)==trace(y2(:,:,1))*eye(2)/2;
PartialTrace(y2(:,:,2),1)==trace(y2(:,:,2))*eye(2)/2;

trace(x1(:,:,1)+x1(:,:,2)+y1(:,:,1)+y1(:,:,2))==1;
trace(x2(:,:,1)+x2(:,:,2)+y2(:,:,1)+y2(:,:,2))==1;

trace(x1(:,:,1))==trace(x1(:,:,1)+x1(:,:,2))*1/2;
trace(x1(:,:,2))==trace(x1(:,:,1)+x1(:,:,2))*1/2;
trace(x2(:,:,1))==trace(x2(:,:,1)+x2(:,:,2))*1/2;
trace(x2(:,:,2))==trace(x2(:,:,1)+x2(:,:,2))*1/2;




for x=1:2
    for a=1:2
        s=0;
        for y=1:2
            for b=1:2
                s=s+2*PartialTrace(Tensor(id,transpose(sigma(:,:,a,x)))*eta(:,:,b,y),2)/2;
            end
        end
        rhopost(:,:,a,x)=s;
    end
end

    
s=0;
for x=1:2
    for a=1:2
        for y=1:2
            s=s+T(x,y)*2*(-1)^(a+1)*trace(Tensor(id,transpose(sigma(:,:,a,x)))*(eta(:,:,1,y)-eta(:,:,2,y)));
        end
    end
end
S1=s;


obj=0;
for x=1:2
    for a=1:2
        for z=1:2
            obj=obj+T(x,z)*(-1)^(a+1)*trace(rhopost(:,:,a,x)*M(:,:,z));
        end
    end
end
S2=obj;
S1>=alpha;

maximise(real(S2))

cvx_end
cvx_optval;


clear rhopost

cvx_begin sdp quiet
cvx_solver mosek
variable sigma(2,2,2,2) symmetric semidefinite

T0=sigma(:,:,1,1)+sigma(:,:,2,1);
sigma(:,:,1,1)+sigma(:,:,2,1)==sigma(:,:,1,2)+sigma(:,:,2,2);
trace(sigma(:,:,1,1)+sigma(:,:,2,1))==1;





for x=1:2
    for a=1:2
        s=0;
        for y=1:2
            for b=1:2
                s=s+2*PartialTrace(Tensor(id,transpose(sigma(:,:,a,x)))*eta(:,:,b,y),2)/2;
            end
        end
        rhopost(:,:,a,x)=s;
    end
end

    
s=0;
for x=1:2
    for a=1:2
        for y=1:2
            s=s+T(x,y)*2*(-1)^(a+1)*trace(Tensor(id,transpose(sigma(:,:,a,x)))*(eta(:,:,1,y)-eta(:,:,2,y)));
        end
    end
end
S1=s;


obj=0;
for x=1:2
    for a=1:2
        for z=1:2
            obj=obj+T(x,z)*(-1)^(a+1)*trace(rhopost(:,:,a,x)*M(:,:,z));
        end
    end
end
S2=obj;
S1>=alpha;

maximise(real(S2))
cvx_end
cvx_optval;



clear M

cvx_begin sdp quiet
cvx_solver mosek
variable MM(2,2,2,2) hermitian semidefinite


M(:,:,1)=MM(:,:,1,1)-MM(:,:,2,1);
M(:,:,2)=MM(:,:,1,2)-MM(:,:,2,2);

for y=1:2
    s=0;
    for b=1:2
        s=s+MM(:,:,b,y);
    end
    s==id;
end



obj=0;
for x=1:2
    for a=1:2
        for z=1:2
            obj=obj+T(x,z)*(-1)^(a+1)*trace(rhopost(:,:,a,x)*M(:,:,z));
        end
    end
end
S2=obj;
S1>=alpha;

maximise(real(S2))
cvx_end
cvx_optval;
end


ll=[ll, cvx_optval]
end


dat=[dat; alpha max(ll)]
end





%% 

tau111=x2(:,:,1);
tau211=x2(:,:,2);
tau120=y2(:,:,1);
tau202=y2(:,:,2);
q11=trace(x2(:,:,1)+x2(:,:,2));
q20=trace(y2(:,:,1));
q02=trace(y2(:,:,2));





