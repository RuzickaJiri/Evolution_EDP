
%% First equation

% définition de taille
x_length = 10 ;
points = 100 ;
h = x_length/points ;
Dt = 0.01 ;

tfinal = 200000 ;
nt = tfinal/Dt ;

phenotypes_number = 3 ;
phe_diff = zeros(phenotypes_number, 1) ;
phe_diff(1) = 0.1 ;
phe_diff(2) = 0.2 ;
phe_diff(3) = 0.3 ;

E = zeros(phenotypes_number, points) ;
 
E(1,2) = 1 ;
E(2,9) = 1 ;
E(3,5) = 1 ;

a = ones(points,1) ;
A = zeros(points, 1) ;

for n=1:points
    A(n) = a(n) - sum(E(:,n)) ;
end

M = zeros(points, points, phenotypes_number) ;

for n=2:points-1
    for i=1:phenotypes_number
        M(n,n,i) = 1 + Dt*phe_diff(i)*(-2)/(h*h) + A(n) ;
        M(n+1,n,i) = Dt*phe_diff(i)/(h*h) ;
        M(n-1,n,i) = Dt*phe_diff(i)/(h*h) ;
    end
    M(1,1,i) = 1 + A(1);
    M(points,1,i) = 1 + A(1);
end

t = 0 ;
while t < tfinal 
    for n=1:points
    A(n) = a(n) - sum(E(:,n)) ;
    end


for n=2:points-1
    for i=1:phenotypes_number
        M(n,n,i) = 1 + Dt*phe_diff(i)*(-2)/(h*h) + A(n) ;
    end
    M(1,1,i) = 1 + A(1);
    M(points,1,i) = 1 + A(1);
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

E

