function [results] = dsge_practical_forward_errors(inputs)
%Matrix quadratic backward errors and conditioning number follow Higham and Kim (2001)
%0=A*X^2+B*X+C
%X is the n by n solvent
%A, B, and C are n by n matrices
%Alexander Meyer-Gohde
%24/03/2022
if inputs.nstatic>0
A=[inputs.A_static;inputs.AA];%sparse(inputs.A_full);
B=[inputs.B_static inputs.B_rest;zeros(inputs.ndynamic,inputs.nstatic) inputs.BB];%sparse(inputs.B_full);
C=[inputs.C_static;inputs.CC];%sparse(inputs.C_full);
else
    A=inputs.AA;
    B=inputs.BB;
    C=inputs.CC;
end
A=[zeros(inputs.endo_nbr,inputs.npred+inputs.nstatic) A];
C=[zeros(inputs.endo_nbr,inputs.nstatic) C zeros(inputs.endo_nbr,inputs.nfwrd)];

P=inputs.X;
F=A*P+B;
results=NaN(8,5);
P_F=norm(P,'fro');
R_P=F*P+C;
R_P_F=norm(R_P,'fro');
[ny,~]=size(P);
if ny<51
V=kron(eye(ny),F)+kron(P',A);
results(4,1)=1/min(svd(V));
[temp_P,~]=linsolve(V,R_P(:));
results(7,1)=norm(temp_P)/P_F;
results(8,1)=results(4,1)*(R_P_F/P_F); %CHECK
%results1=results;
else
[A1,D1,P,Q] = qz(full(F),full(A));
[U,B1] = schur(full(-P),'complex');

C1=zeros(size(A1,1),size(B1,1));
F1=zeros(size(A1,1),size(B1,1));
E1=eye(size(B1));
[~, ~, dif, ~, ~, ~] = ztgsyl('N',3,complex(A1),complex(B1),complex(C1),complex(D1),complex(E1),complex(F1));
results(4,1)=1/dif;

C1=P*R_P*U;
[R, ~, ~, ~, ~, ~] = ztgsyl('N',0,complex(A1),complex(B1),complex(C1),complex(D1),complex(E1),complex(F1));
temp2=Q*R*U';
%temp2=temp2(:);
%results(7,1)=normest(temp2)/P_F;
results(7,1)=norm(temp2,'fro')/P_F;
results(8,1)=results(4,1)*R_P_F/P_F;
end
%[results1(4,1) results1(7,1) results1(8,1) ;results(4,1) results(7,1) results(8,1)]



