%% Generate a N-dimensional matrix which is used for recursion formula of edge states
clc
clear all

% Define parameter
N = 25; %number of layers of the lattice
M = zeros(4*N); %create an all-zero N*N matrix
t1 = 1.0;
t2 = 1.0;
t3 = 0.3; %constants from the Hamiltonian
%t2=0
%t3=0

step=100;%total footsteps of kx value
delta_kx=2*pi/step;%intervals of kx value
%kx=-pi/2; %initial value of kx value

kx = -pi;

res=zeros(4*N,step);% the result recording matrix
% Generate the complete matrix
%The first layer, i.e. the boundary condition 



%% Calculate the eigenvalue of the matrix

% Display setting for the screen print
%name = 'Alice';
%age = 12;
%ponsn = ['Calculating the number ',num2str(l),' eigenvalue now.'];
%disp(ponsn)

% Loop for discretized kx value from -pi to +pi and record the data
%l=1;
for l=1:step
kx=-pi+delta_kx*l;
%M=newmat(kx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Renewing the M matrix %%%%%%%%%%%%%%%%%%%%%%%%
M(1,2)=t1+1*1i*t3*exp(-1i*kx);
M(1,4)=t1-1*1i*t3*exp(-1i*kx);
M(1,3)=t2*exp(-1i*kx);
M(1,6)=1i*t3;

M(2,1)=t1-1*1i*t3*exp(1i*kx);
M(2,3)=t1+1*1i*t3*exp(-1i*kx);

M(3,2)=t1-1*1i*t3*exp(1i*kx);
M(3,1)=t2*exp(1i*kx);
M(3,4)=t1+1*1i*t3*exp(1i*kx);
M(3,6)=-1i*t3;

M(4,1)=t1+1*1i*t3*exp(1i*kx);
M(4,5)=1i*t3;
M(4,7)=-1i*t3;
M(4,3)=t1-1*1i*t3*exp(-1i*kx);
M(4,6)=t2;

%The last layer, i.e. the upward boundary
M(4*(N-1)+1,4*(N-1)+4)=t1-1*1i*t3*exp(-1i*kx);
M(4*(N-1)+1,4*(N-1)+3)=t2*exp(-1i*kx);
M(4*(N-1)+1,4*(N-1)+2)=t1+1*1i*t3*exp(-1i*kx);
M(4*(N-1)+1,4*(N-2)+4)=-1i*t3;

M(4*(N-1)+2,4*(N-1)+1)=t1-1*1i*t3*exp(1i*kx);
M(4*(N-1)+2,4*(N-1)+3)=t1+1*1i*t3*exp(-1i*kx);
M(4*(N-1)+2,4*(N-2)+4)=t2;
M(4*(N-1)+2,4*(N-2)+1)=-1i*t3;
M(4*(N-1)+2,4*(N-2)+3)=1i*t3;

M(4*(N-1)+3,4*(N-1)+4)=t1+1*1i*t3*exp(1i*kx);
M(4*(N-1)+3,4*(N-1)+1)=t2*exp(1i*kx);
M(4*(N-1)+3,4*(N-1)+2)=t1-1*1i*t3*exp(1i*kx);
M(4*(N-1)+3,4*(N-2)+4)=1i*t3;

M(4*(N-1)+4,4*(N-1)+1)=t1+1*1i*t3*exp(1i*kx);
M(4*(N-1)+4,4*(N-1)+3)=t1-1*1i*t3*exp(-1i*kx);

%The rest entries of the bulk

for n=2:N-1
M(4*(n-1)+1,4*(n-1)+4)=t1-1*1i*t3*exp(-1i*kx);
M(4*(n-1)+1,4*(n-1)+3)=t2*exp(-1i*kx);
M(4*(n-1)+1,4*(n-1)+2)=t1+1*1i*t3*exp(-1i*kx);
M(4*(n-1)+1,4*(n-2)+4)=-1i*t3;
M(4*(n-1)+1,4*n+2)=1i*t3;

M(4*(n-1)+2,4*(n-1)+1)=t1-1*1i*t3*exp(1i*kx);
M(4*(n-1)+2,4*(n-1)+3)=t1+1*1i*t3*exp(-1i*kx);
M(4*(n-1)+2,4*(n-2)+4)=t2;
M(4*(n-1)+2,4*(n-2)+1)=-1i*t3;
M(4*(n-1)+2,4*(n-2)+3)=1i*t3;

M(4*(n-1)+3,4*(n-1)+4)=t1+1*1i*t3*exp(1i*kx);
M(4*(n-1)+3,4*(n-1)+1)=t2*exp(1i*kx);
M(4*(n-1)+3,4*(n-1)+2)=t1-1*1i*t3*exp(1i*kx);
M(4*(n-1)+3,4*(n-2)+4)=1i*t3;
M(4*(n-1)+3,4*n+2)=-1i*t3;

M(4*(n-1)+4,4*(n-1)+1)=t1+1*1i*t3*exp(1i*kx);
M(4*(n-1)+4,4*(n-1)+3)=t1-1*1i*t3*exp(-1i*kx);
M(4*(n-1)+4,4*n+2)=t2;
M(4*(n-1)+4,4*n+1)=1i*t3;
M(4*(n-1)+4,4*n+3)=-1i*t3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res(:,l)=eig(M);
% Display setting for the screen print
%name = 'Alice';
%ll = l;
ponsn = ['Calculating the number ',num2str(l),' eigenvalue now.'];
disp(ponsn)
end

disp('End of the eigenvalue calculation and recording.')
% Calculatee the absolute value of the eigenvalues
%Abs_res=abs(res);
disp('End of calculating eigenvalue data.')

%% Extract the desired data for edge states to be plotted in Mathematica
temp_res=res';
temp_data=temp_res(:,17:34);
x=-pi+2*pi/100:2*pi/100:pi;

data17=zeros(100,2);
data18=zeros(100,2);
data19=zeros(100,2);
data20=zeros(100,2);
data21=zeros(100,2);
data22=zeros(100,2);
data23=zeros(100,2);
data24=zeros(100,2);
data25=zeros(100,2);
data26=zeros(100,2);
data27=zeros(100,2);
data28=zeros(100,2);
data29=zeros(100,2);
data30=zeros(100,2);
data31=zeros(100,2);
data32=zeros(100,2);
data33=zeros(100,2);
data34=zeros(100,2);
%data35=zeros(100,2);

data17(:,1)=x';
data18(:,1)=x';
data19(:,1)=x';
data20(:,1)=x';
data21(:,1)=x';
data22(:,1)=x';
data23(:,1)=x';
data24(:,1)=x';
data25(:,1)=x';
data26(:,1)=x';
data27(:,1)=x';
data28(:,1)=x';
data29(:,1)=x';
data30(:,1)=x';
data31(:,1)=x';
data32(:,1)=x';
data33(:,1)=x';
data34(:,1)=x';

data17(:,2)=temp_data(:,1);
data18(:,2)=temp_data(:,2);
data19(:,2)=temp_data(:,3);
data20(:,2)=temp_data(:,4);
data21(:,2)=temp_data(:,5);
data22(:,2)=temp_data(:,6);
data23(:,2)=temp_data(:,7);
data24(:,2)=temp_data(:,8);
data25(:,2)=temp_data(:,9);
data26(:,2)=temp_data(:,10);
data27(:,2)=temp_data(:,11);
data28(:,2)=temp_data(:,12);
data29(:,2)=temp_data(:,13);
data30(:,2)=temp_data(:,14);
data31(:,2)=temp_data(:,15);
data32(:,2)=temp_data(:,16);
data33(:,2)=temp_data(:,17);
data34(:,2)=temp_data(:,18);

data_test=temp_data(:,3)';

%save('mdata.mat','temp_data'); % save the data as a mat file
save('mdata17.mat','data17'); 
save('mdata18.mat','data18'); 
save('mdata19.mat','data19'); 
save('mdata20.mat','data20'); 
save('mdata21.mat','data21'); 
save('mdata22.mat','data22'); 
save('mdata23.mat','data23'); 
save('mdata24.mat','data24'); 
save('mdata25.mat','data25'); 
save('mdata26.mat','data26'); 
save('mdata27.mat','data27'); 
save('mdata28.mat','data28'); 
save('mdata29.mat','data29'); 
save('mdata30.mat','data30'); 
save('mdata31.mat','data31'); 
save('mdata32.mat','data32'); 
save('mdata33.mat','data33'); 
save('mdata34.mat','data34'); 


save('test_data.mat','data_test');
%% Plot the edge states

figure,plot(res(17,:))
hold on;
plot(res(18,:))
hold on;
plot(res(19,:))
hold on;
plot(res(20,:))
hold on;
plot(res(21,:))
hold on;
plot(res(22,:))
hold on;
plot(res(23,:))
hold on;
plot(res(24,:))
hold on;
plot(res(25,:))
hold on;
plot(res(26,:))
hold on;
plot(res(27,:))
hold on;
plot(res(28,:))
hold on;
plot(res(29,:))
hold on;
plot(res(30,:))
hold on;
plot(res(31,:))
hold on;
plot(res(32,:))
hold on;
plot(res(33,:))
hold on;
plot(res(34,:))

%  figure, plot(res(25,:))