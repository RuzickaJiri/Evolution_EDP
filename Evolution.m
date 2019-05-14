
%% First equation

% parameters
h = x_length/points ; % space step
x_length = 10 ; % space dimension
points = 100 ; % iteration number space
Dt = 0.01 ; % time step
tfinal = 200000 ; %final time
nt = tfinal/Dt ; % iteration number time

%Diffusion constants
phenotypes_number = 3 ; %number of phenotypes
phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
phe_diff(1) = 0.17 ;
phe_diff(2) = 0.2 ;
phe_diff(3) = 0.25 ;

%Environment matrix
E = zeros(phenotypes_number, points) ;
E(1,2) = 1 ;
E(2,9) = 1 ;
E(3,5) = 1 ;

%Sources vectors
a = ones(points,1) ;
A = zeros(points, 1) ;

for n=1:points
    A(n) = a(n) - sum(E(:,n)) ;
end

% Discretization matrix
M = zeros(points, points, phenotypes_number) ;

for i=1:phenotypes_number
    for n=2:points-1 
        M(n,n,i) = 1 + Dt*(phe_diff(i)*(-2)/(h*h) + A(n)) ;
        M(n+1,n,i) = Dt*phe_diff(i)/(h*h) ;
        M(n-1,n,i) = Dt*phe_diff(i)/(h*h) ;
    end
    M(1,1,i) = 1 + A(1);
    M(points,1,i) = 1 + A(1);
end

% Evolution 
t = 0 ;
while t < tfinal 
    for n=1:points
    A(n) = a(n) - sum(E(:,n)) ;
    end


for i=1:phenotypes_number
    for n=2:points-1
        M(n,n,i) = 1 + Dt*(phe_diff(i)*(-2)/(h*h) + A(n)) ;
    end
    M(1,1,i) = 1 + A(1);
    M(points,points,i) = 1 + A(points);
end

for i=1:phenotypes_number
    E(i,:) = E(i,:)*M(:,:,i);
    
end
if mod(t,10)==0
figure(1) ;


plot(1:points, E(1,1:points)) ;
hold on
plot(1:points, E(2,1:points)) ;
plot(1:points, E(3,1:points)) ;
hold off
end

t=t+1;
end


%% Equation with mutations


% parameters
h = x_length/points ; % space step
x_length = 10 ; % space dimension
points = 100 ; % iteration number space
Dt = 0.01 ; % time step
tfinal = 200000 ; %final time
nt = tfinal/Dt ; % iteration number time

%Diffusion constants
phenotypes_number = 3 ; %number of phenotypes
phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
phe_diff(1) = 0.17 ;
phe_diff(2) = 0.2 ;
phe_diff(3) = 0.25 ;

%Environment matrix
E = zeros(phenotypes_number, points) ;
E(1,2) = 1 ;
E(2,9) = 1 ;
E(3,5) = 1 ;

%Sources vectors
a = ones(points,1) ;
A = zeros(points, 1) ;

for n=1:points
    A(n) = a(n) - sum(E(:,n)) ;
end

% Discretization matrix
M = zeros(points, points, phenotypes_number) ;

for i=1:phenotypes_number
    for n=2:points-1 
        M(n,n,i) = 1 + Dt*phe_diff(i)*(-2)/(h*h) + A(n) ;
        M(n+1,n,i) = Dt*phe_diff(i)/(h*h) ;
        M(n-1,n,i) = Dt*phe_diff(i)/(h*h) ;
    end
    M(1,1,i) = 1 + A(1);
    M(points,points,i) = 1 + A(points);
end

%Mutation matrix
Mutation = zeros(phenotypes_number, phenotypes_number) ;
Mutation(1,1) = -0.5  ;
Mutation(1,2) = 0.0 ;
Mutation(1,3) = 0.0 ;
Mutation(2,1) = 0.0;
Mutation(2,2) = -0.5 ;
Mutation(2,3) = 0 ;
Mutation(3,1) = 0.5 ;
Mutation(3,2) = 0.5 ;
Mutation(3,3) = 0 ;

epsilon = 0.1 ;




% Evolution 
t = 0 ;
while t < tfinal 
    for n=1:points
    A(n) = a(n) - sum(E(:,n)) ;
    end
    
    for i=1:phenotypes_number
        for n=2:points-1 
            MM=0;
            for k=1:phenotypes_number
                MM = MM+Mutation(i,k)*E(k,n) ;
            end 
            M(n,n,i) = 1 + Dt*(phe_diff(i)*(-2)/(h*h) + A(n) +  epsilon*MM);
        end
        MM=0;
        for k=1:phenotypes_number
             MM = MM+Mutation(i,k)*E(k,1) ;
        end
        M(1,1,i) = 1 + A(1)+epsilon*MM;
        MM=0;
        for k=1:phenotypes_number
             MM = MM+Mutation(i,k)*E(k,points) ;
        end
        M(points,points,i) = 1 + A(points)+epsilon*MM;
    end 

    

    for i=1:phenotypes_number
        E(i,:) = E(i,:)*M(:,:,i);
    end
    if mod(t,100)==0
    figure(1) ;

    plot(1:points, E(1,1:points)) ;
    hold on
    plot(1:points, E(2,1:points)) ;
    plot(1:points, E(3,1:points)) ;
    hold off
    end

    t=t+1;
end

