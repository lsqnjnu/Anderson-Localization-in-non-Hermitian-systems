clear
delete(gcp)
parpool('local',48)

%define the hopping matrix on each lattice
V = 0;
Tx = -1;
Ty = -1;
Tz = -1;

NMAX=10^5;%maximum steps
NORTHO=8;%renomalization in every NORTHO steps 
EPS=0.01;%precision
E=0;%Fermi energy
%the tunable parameter
%Wn=1.5:0.5:5.5;   %disorder strength
Wn=6.2:0.01:6.6;
% Nn=[4,5]; 
Nn=[10,12,14,16];%size of quasi-1d sample
%localization length and covariance
LAMBDAE1=zeros(length(Nn),length(Wn)); %localization length
LAMBDAE2=zeros(length(Nn),length(Wn)); %variance of localiztion length
for nn=1:length(Nn)
Nx=Nn(nn);
Ny=Nx; %height and width to be the same 
%h00 hopping matrix of intra-line in y-direction;
%h01x hopping matrix of inter-line in x-direction;
h00=kron(eye(Ny),V)+kron(diag(ones(1,Ny-1),1),Ty')+kron(diag(ones(1,Ny-1),-1),Ty)+kron(diag(1,-Ny+1),Ty')+kron(diag(1,Ny-1),Ty);
h01x=kron(eye(Ny),Tx);
h01z=kron(eye(Ny),Tz);

H00=kron(eye(Nx),h00)+kron(diag(ones(1,Nx-1),1),h01x)+kron(diag(ones(1,Nx-1),-1),h01x')+kron(diag(1,-Nx+1),h01x)+kron(diag(1,Nx-1),h01x');
H01=kron(eye(Nx),h01z);


wid=Nx*Ny;
%H00 onsite matrix; H01 is the hopping matrix

parfor en=1:length(Wn)
        tic;
        W=Wn(en)
        [lambda,ACC_LAMB] = loca2(E,W,wid,H00,H01,NMAX,EPS,NORTHO );        
        LAMBDAE1(nn,en)=lambda/Ny;%localization length 
        LAMBDAE2(nn,en)=ACC_LAMB;%variance of localization length
        
end
save localizationscale LAMBDAE1 LAMBDAE2  Wn Nn 
end

