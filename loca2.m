function[lambda,ACC_LAMB] = loca2(E,W,wid,H00,H01,NMAX,EPS,NORTHO )
%NMAX the maximum steps and EPS is the precision
%renomalization in every NORTHO steps 
%W is the disorder strength
%wid is the width of the quasi-1d system
%H00 onsite matrix and H01 is the hopping matrix
invH01=inv(H01);
GAMMA=0;GAMMA2=0;%intermediate parameter
Q0=eye(2*wid); %initialized unitary matrix
T12=-invH01*H01';%intermeidate parameter
for N=1:NMAX/NORTHO
    %introduce onsite disorder
    Wr=5;
    HW = (Wr*rand(wid, NORTHO) - Wr/2) + 1i*(W*rand(wid, NORTHO) - W/2);
    %interant method to calculate the localization length
    for j=1:NORTHO
        %onsite matrix including  disorder
        H11=H00+diag( HW(:,j) );
        % the transfer matrix 
        TT=[H01\(eye(1*wid)*E-H11),T12;eye(1*wid),zeros(1*wid)];
        Q0=TT*Q0;
    end
    [Q,R]=qr(Q0);%QR decomposition
    Q=Q*diag(sign(diag(R)));%sign(diag(R)) is used to fix the sign 
    R=diag(sign(diag(R)))*R;%diag(R) is fixed to be positive   将Q和R固定为正
    Q0=Q;
    D=log(R(1*wid,1*wid));   %取最小对角元(中间)
    
    GAMMA=GAMMA+D;
    GAMMA2=GAMMA2+D*D;
    lambda=N*NORTHO/GAMMA;   %这里的lambda还并未约化
    ACC_LAMB=abs(GAMMA2/GAMMA^2-1/N);
    ACC_LAMB=sqrt(ACC_LAMB);%variance     
    xtt=mod(N,10000);
    if (N>10^4) && ACC_LAMB<EPS % EPS is the precision
        break
    end
    if xtt == 3
        N
        lambda
        ACC_LAMB
        toc
    end
end
end