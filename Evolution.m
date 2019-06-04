
% %% First equation (faux)
% 
% % parameters
% x_length = 10 ; % space dimension
% points = 300 ; % iteration number space
% Dt = 0.0001 ; % time step
% tfinal = 20000 ; %final time
% nt = tfinal/Dt ; % iteration number time
% h = x_length/points ; % space step
% 
% %Diffusion constants
% phenotypes_number = 3 ; %number of phenotypes
% phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
% phe_diff(1) = 0.17 ;
% phe_diff(2) = 0.2 ;
% phe_diff(3) = 0.25 ;
% 
% %Environment matrix
% E = zeros(phenotypes_number, points) ;
% E(1,1) = 1 ;
% E(2,1) = 1 ;
% E(3,1) = 1 ;
% 
% %Sources vectors
% a = ones(points,1) ;
% A = zeros(points, 1) ;
% 
% for n=1:points
%     A(n) = a(n) - sum(E(:,n)) ;
% end
% 
% % Discretization matrix
% M = zeros(points, points, phenotypes_number) ;
% 
% for i=1:phenotypes_number
%     for n=2:points-1 
%         M(n,n,i) = 1 + Dt*(phe_diff(i)*(-2)/(h*h) + A(n)) ;
%         M(n+1,n,i) = Dt*phe_diff(i)/(h*h) ;
%         M(n-1,n,i) = Dt*phe_diff(i)/(h*h) ;
%     end
%     M(1,1,i) = 1 + A(1);
%     M(points,1,i) = 1 + A(1);
% end
% 
% % Evolution 
% t = 0 ;
% while t < tfinal 
%     for n=1:points
%     A(n) = a(n) - sum(E(:,n)) ;
%     end
% 
% 
% for i=1:phenotypes_number
%     for n=2:points-1
%         M(n,n,i) = 1 + Dt*(phe_diff(i)*(-2)/(h*h) + A(n)) ;
%     end
%     M(1,1,i) = 1 + A(1);
%     M(points,points,i) = 1 + A(points);
% end
% 
% for i=1:phenotypes_number
%     E(i,:) = E(i,:)*M(:,:,i);
%     
% end
% if mod(t,10)==0
% figure(1) ;
% 
% 
% plot(1:points, E(1,1:points)) ;
% hold on
% plot(1:points, E(2,1:points)) ;
% plot(1:points, E(3,1:points)) ;
% hold off
% end
% 
% t=t+1;
% end
% 
% 
% %% Equation with mutations (faux)
% 
% 
% % parameters
% x_length = 10 ; % space dimension
% points = 100 ; % iteration number space
% Dt = 0.001 ; % time step
% tfinal = 200000 ; %final time
% nt = tfinal/Dt ; % iteration number time
% h = x_length/points ; % space step
% 
% %Diffusion constants
% phenotypes_number = 3 ; %number of phenotypes
% phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
% phe_diff(1) = 0.17 ;
% phe_diff(2) = 0.2 ;
% phe_diff(3) = 0.25 ;
% 
% %Environment matrix
% E = zeros(phenotypes_number, points) ;
% E(1,2) = 1 ;
% E(2,9) = 1 ;
% E(3,5) = 1 ;
% 
% %Sources vectors
% a = ones(points,1) ;
% A = zeros(points, 1) ;
% 
% for n=1:points
%     A(n) = a(n) - sum(E(:,n)) ;
% end
% 
% % Discretization matrix
% M = zeros(points, points, phenotypes_number) ;
% 
% for i=1:phenotypes_number
%     for n=2:points-1 
%         M(n,n,i) = 1 + Dt*(phe_diff(i)*(-2)/(h*h) + A(n)) ;
%         M(n+1,n,i) = Dt*phe_diff(i)/(h*h) ;
%         M(n-1,n,i) = Dt*phe_diff(i)/(h*h) ;
%     end
%     M(1,1,i) = 1 + A(1);
%     M(points,points,i) = 1 + A(points);
% end
% 
% %Mutation matrix
% Mutation = zeros(phenotypes_number, phenotypes_number) ;
% Mutation(1,1) = -1;
% Mutation(1,2) = 1 ;
% Mutation(1,3) = 0;
% Mutation(2,1) = 0.0;
% Mutation(2,2) = -1 ;
% Mutation(2,3) = 1 ;
% Mutation(3,1) = 1 ;
% Mutation(3,2) = 0.0 ;
% Mutation(3,3) = -1 ;
% 
% epsilon = 0.01 ;
% 
% 
% 
% 
% % Evolution 
% t = 0 ;
% while t < tfinal 
%     for n=1:points
%         A(n) = a(n) - sum(E(:,n)) ;
%     end
%     
%     for i=1:phenotypes_number
%         for n=2:points-1 
%             MM=0;
%             for k=1:phenotypes_number
%                 MM = MM+Mutation(i,k)*E(k,n) ;
%             end 
%             M(n,n,i) = 1 + Dt*(phe_diff(i)*(-2)/(h*h) + A(n) +  epsilon*MM);
%         end
%         MM=0;
%         for k=1:phenotypes_number
%              MM = MM+Mutation(i,k)*E(k,1) ;
%         end
%         M(1,1,i) = 1 + A(1)+epsilon*MM;
%         MM=0;
%         for k=1:phenotypes_number
%              MM = MM+Mutation(i,k)*E(k,points) ;
%         end
%         M(points,points,i) = 1 + A(points)+epsilon*MM;
%     end 
% 
%     
% 
%     for i=1:phenotypes_number
%         E(i,:) = E(i,:)*M(:,:,i);
%     end
%     if mod(t,100)==0
%     figure(1) ;
% 
%     plot(1:points, E(1,1:points)) ;
%     hold on
%     plot(1:points, E(2,1:points)) ;
%     plot(1:points, E(3,1:points)) ;
%     hold off
%     end
% 
%     t=t+1;
% end
% 
% %% Equation 2D.(faux)
% 
% % parameters
% x_length = 10 ; % space dimension
% points = 20 ; % iteration number space
% Dt = 0.001 ; % time step
% tfinal = 20000 ; %final time
% nt = tfinal/Dt ; % iteration number time
% h = x_length/points ; % space step
% 
% %Diffusion constants
% phenotypes_number = 3 ; %number of phenotypes
% phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
% phe_diff(1) = 0.17 ;
% phe_diff(2) = 0.2 ;
% phe_diff(3) = 0.25 ;
% 
% %Environment matrix
% E = zeros(phenotypes_number, points^2) ;
% for i=1:phenotypes_number
%     for n =1:points^2
%         E(i,n)= rand;
%     end
% end
% 
% 
% %Sources vectors
% [X,Y]= meshgrid(1:points,1:points);
% a = ones(1,points^2);
% A = zeros(points^2, 1) ;
% 
% for n=1:points^2
%     A(n) = a(n) - sum(E(:,n)) ;
% end
% 
% % Discretization matrix
% M = zeros(points^2, points^2, phenotypes_number) ;
% 
% for i=1:phenotypes_numberquation de 
%     for n=1:points^2
%         M(n,n,i) = 1 + Dt*phe_diff(i)*(-4)/(h*h) + A(n);
%         
%     end
%     for n=0:points-1
%         for m=1:points
%             if m<points
%                 M(n*points+m+1,n*points+m,i) = Dt*phe_diff(i)/(h*h) ;
%             end
%             if m>1
%                 M(n*points+m-1,n*points+m,i) = Dt*phe_diff(i)/(h*h) ;
%             end
%             if n<points-1
%                 M((n+1)*points+m,n*points+m,i) = Dt*phe_diff(i)/(h*h) ;
%             end
%             if n>0
%                 M((n-1)*points+m,n*points+m,i) = Dt*phe_diff(i)/(h*h) ; 
%             end
%         end
%     end
% end
% 
% %Mutation matrix
% Mutation = zeros(phenotypes_number, phenotypes_number) ;
% Mutation(1,1) = -2;
% Mutation(1,2) = 0.5 ;
% Mutation(1,3) = 0;
% Mutation(2,1) = 0.5;
% Mutation(2,2) = -1.5 ;
% Mutation(2,3) = 0 ;
% Mutation(3,1) = 0.5 ;
% Mutation(3,2) = 0.5 ;
% Mutation(3,3) = -1 ;
% 
% 
% Mutation = zeros(phenotypes_number, phenotypes_number) ;
% epsilon = 0.01 ;
% 
% 
% 
% 
% % Evolution 
% t = 0 ;
% while t < tfinal 
%     for n=1:points^2
%     A(n) = a(n) - sum(E(:,n)) ;
%     end
%     
%     for i=1:phenotypes_number
%         for n=1:points^2 
%             MM=0;
%             for k=1:phenotypes_number
%                 MM = MM+Mutation(i,k)*E(k,n) ;
%             end 
%             M(n,n,i) = 1 + Dt*(phe_diff(i)*(-4)/(h*h) + A(n) +  epsilon*MM);
%         end
%     end 
%     for i=1:phenotypes_number
%         E(i,:) = E(i,:)*M(:,:,i);
%     end,'FaceAlpha',0.1
%     if mod(t,50)==0 
%         figure(1)
%         for i=1:phenotypes_number
%             x=1:points;
%             y=1:points;
%             C0=0;
%             CO(:,:,3) = ones(points)*(i)/(phenotypes_number);
%             CO(:,:,1) = ones(points)*(phenotypes_number-i)/(phenotypes_number);
%             CO(:,:,2) = ones(points)*(phenotypes_number-i*0.5)/(phenotypes_number);
%             surf(x,y,reshape(E(i,:),points,points),CO) ;
%             hold on
%         end
%         hold off
%     end
% 
%     t=t+1;
% end 
%% Equation with mutations rewriten (completly corect)


% parameters
x_length = 10 ; % space dimension
points = 100 ; % iteration number space
Dt = 0.001 ; % time step
tfinal = 400 ; %final time
nt = tfinal/Dt ; % iteration number time
h = x_length/points ; % space step
afficher=true;

%Diffusion constants
phenotypes_number = 3 ; %number of phenotypes
phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
phe_diff(1) = 0.1 ;
phe_diff(2) = 0.4 ;
phe_diff(3) = 0.9 ;

%Environment matrix
E = zeros(1,phenotypes_number*points) ;
E(20) = 1 ;
E(points+50) = 1 ;
E(2*points+80) = 1 ;

%Sources vectors
a = -pi/2:2*pi/(points-1):3*pi/2;% fonction de répartition de la nouriture
a= cos(a);
a=a*diag(a)*0.8+0.2;
A=zeros(1,3*phenotypes_number);
for i=0:phenotypes_number-1
A(i*points+1:i*points+points)=a;%répartition de la nouriture par espèce
end


%Mutation matrix
Mutation = zeros(phenotypes_number, phenotypes_number) ;
Mutation(1,1) = -1;
Mutation(1,2) = 0.5 ;
Mutation(1,3) = 0.5;
Mutation(2,1) = 0.5;
Mutation(2,2) = -1;
Mutation(2,3) = 0.5 ;
Mutation(3,1) = 0.5 ;
Mutation(3,2) = 0.5 ;
Mutation(3,3) = -1 ;

epsilon=0.01;

% linear Discretization matrix
M = zeros(points* phenotypes_number, points* phenotypes_number) ;

for  i=0:phenotypes_number-1
    for n=2:points-1
        k=i*points+n;
        % termes liés à la difusion
        M(k,k) = phe_diff(i+1)*(-2)/(h*h) ;
        if n<points
        M(k+1,k) = phe_diff(i+1)/(h*h) ;
        end
        if n>1
        M(k-1,k) = phe_diff(i+1)/(h*h) ;
        end
        for j=0:phenotypes_number-1
            p=j*points+n;
            % termes liés aux mutations
            M(p,k)= M(p,k)+(epsilon*Mutation(i+1,j+1));
        end
    end
    
    % termes aux limites
    M(i*points+1,i*points+1)=-phe_diff(i+1)/(h*h);
    M(i*points+2,i*points+1)= phe_diff(i+1)/(h*h);
    M((i+1)*points,(i+1)*points)=-phe_diff(i+1)/(h*h);
    M((i+1)*points-1,(i+1)*points)= +phe_diff(i+1)/(h*h);
    
    
end

% non linear Discretization matrix
M2 = zeros(points* phenotypes_number, points* phenotypes_number) ;

for i=0:phenotypes_number-1
    for n=1:points
        k=i*points+n;
        for j=0:phenotypes_number-1
            M2(j*points+n,i*points+n)=-1;  % termes liés à la consomation de ressources
        end
    end
    
end







% Evolution 
t = 0 ;
while t < nt
     E=E+ Dt*(E*M+E*(diag(A+E*M2)));
    
    if mod(t,1000)==0 && afficher 
        figure(1)
        Es=0;
        for i=1:phenotypes_number
            plot((h:h:x_length), E((i-1)*points+1:(i-1)*points+points),'Color',[i/phenotypes_number,(phenotypes_number-i)/phenotypes_number,1]) ;
            hold on
            Es=Es+E((i-1)*points+1:(i-1)*points+points);
        end
        plot((h:h:x_length),Es,'r--')
        hold on
        plot((h:h:x_length),a,'g--')
        hold off
    end
    t=t+1;
end

%% Equation with Phenotype aparition


% parameters
x_length = 10 ; % space dimension
points = 100 ; % iteration number space
tfinal = 2500 ; %final time
Dt = tfinal/200000 ; % time step
nt = tfinal/Dt ; % iteration number time
h = x_length/points ; % space step
afficher=false;


phenotypes_number = 10 ; %number max of phenotypes
apparition= 0.2*nt/phenotypes_number;% temps d'aparition des nouveau phénotypes 
phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
phe_diff(1) = 0.2 ;
for i=2:phenotypes_number
    phe_diff(i)=phe_diff(i-1)*0.5; % chaquenouveau phénotype aura une vitesse de difusion plus faible que le précédent
end

%Environment matrix
E = zeros(1,phenotypes_number*points) ;
E(1:points) = 1 ;

%Sources vectors
a = -pi/2:2*pi/(points-1):3*pi/2;% fonction de répartition de la nouriture
a= cos(a);
a=a*diag(a)*0.8+0.2;
A=zeros(1,3*phenotypes_number);
for i=0:phenotypes_number-1
A(i*points+1:i*points+points)=a;%répartition de la nouriture par espèce
end


%Mutation matrix
Mutation = zeros(phenotypes_number, phenotypes_number) ;

epsilon=0.001;


% linear Discretization matrix
M = zeros(points* phenotypes_number, points* phenotypes_number) ;

for  i=0:phenotypes_number-1
    for n=2:points-1
        k=i*points+n;
        % termes liés à la difusion
        M(k,k) = phe_diff(i+1)*(-2)/(h*h) ;
        if n<points
        M(k+1,k) = phe_diff(i+1)/(h*h) ;
        end
        if n>1
        M(k-1,k) = phe_diff(i+1)/(h*h) ;
        end
        for j=0:phenotypes_number-1
            p=j*points+n;
            % termes liés aux mutations
            M(p,k)= M(p,k)+(epsilon*Mutation(i+1,j+1));
        end
    end
    
    % termes aux limites
    M(i*points+1,i*points+1)=-phe_diff(i+1)/(h*h);
    M(i*points+2,i*points+1)= phe_diff(i+1)/(h*h);
    M((i+1)*points,(i+1)*points)=-phe_diff(i+1)/(h*h);
    M((i+1)*points-1,(i+1)*points)= phe_diff(i+1)/(h*h);
    
    
end

% non linear Discretization matrix
M2 = zeros(points* phenotypes_number, points* phenotypes_number) ;

for i=0:phenotypes_number-1
    for n=1:points
        k=i*points+n;
        for j=0:phenotypes_number-1
            M2(j*points+n,i*points+n)=-1;  % termes liés à la consomation de ressources
        end
    end
    
end


C=zeros(200,1000,3);
C(:,1:100-1,1) = 0.55;
C(:,101:200-1,2) = 0.60;
C(:,201:300-1,3) = 0.65;
C(:,301:400-1,1) = 0.70;
C(:,401:500-1,2) = 0.75;
C(:,501:600-1,3) = 0.80;
C(:,601:700-1,1) = 0.85;
C(:,701:800-1,2) = 0.9;
C(:,801:900-1,3) = 0.95;
C(:,901:1000,1) = 1;



% Evolution 

Mem=zeros(nt/1000,phenotypes_number*points);
t = 0 ;

nbphen=1; % nombre de phénotype dans l'experimentation en cours
while t < nt-1
    t=t+1;
    
    if mod(t,1000)==0
        Mem(201-(t/1000 +1),:)=E;
        figure(3)
        surf(Mem,C,'LineStyle',':', 'FaceAlpha',1,'FaceLighting','gouraud') ;

        
    end
    
    if t>apparition
        %E
        %pause(1)
    end
    if mod(t,apparition)==10 && nbphen~=phenotypes_number % aparition d'un nouveau phenotype
        nbphen=nbphen+1;
        nbphen
        Mold=Mutation;%old mutation Matrix
        for i=1:nbphen %nouvelle Matrice de mutation avec les hypothèses suivantes: equiprobabilité de mutaton, conservation de la proportion de mutants
            for j=1:nbphen
                if i == j
                    Mutation(i,j)=-1;
                else
                    Mutation(i,j)=1/(nbphen-1);
                end
                
            end
            
        end
        for i=0:phenotypes_number-1
            for n=1:points
                k=i*points+n;
                for j=0:phenotypes_number-1
                p=j*points+n;
                M(p,k)= M(p,k)+epsilon*(Mutation(i+1,j+1)-Mold(i+1,j+1));% on recalcule la nouvelle matrice M pour la nouvelle matrice des mutations
                end
            end
    
        end
        
        
        
        
    end
      
    

    E(:) =E+ Dt*(E*M+E*(diag(A+E*M2)));
    
    if mod(t,1000)==0 && afficher 
        figure(1)
        Es=0;
        for i=1:phenotypes_number
            plot((h:h:x_length), E((i-1)*points+1:(i-1)*points+points),'Color',[i/phenotypes_number,(phenotypes_number-i)/phenotypes_number,1]) ;
            hold on
            Es=Es+E((i-1)*points+1:(i-1)*points+points);
        end
        plot((h:h:x_length),Es,'r--')
        hold on
        plot((h:h:x_length),a,'g--')
        hold off
        
    end
    
    
    
end

figure(1)
Es=0;
for i=1:phenotypes_number
    plot((h:h:x_length), E((i-1)*points+1:(i-1)*points+points),'Color',[(mod(i,3)==1)*i/phenotypes_number,(mod(i,3)==2)*i/phenotypes_number,(mod(i,3)==0)*i/phenotypes_number]) ;
    hold on
    Es=Es+E((i-1)*points+1:(i-1)*points+points);
end
plot((h:h:x_length),Es,'r--')
hold on
plot((h:h:x_length),a,'g--')
hold off

figure(3)
    C(:,1:100,1) = ones(200,100)*0.55;
    C(:,101:200,2) = ones(200,100)*0.60;
    C(:,201:300,3) = ones(200,100)*0.65;
    C(:,301:400,1) = ones(200,100)*0.70;
    C(:,401:500,2) = ones(200,100)*0.75;
    C(:,501:600,3) = ones(200,100)*0.80;
    C(:,601:700,1) = ones(200,100)*0.85;
    C(:,701:800,2) = ones(200,100)*0.9;
    C(:,801:900,3) = ones(200,100)*0.95;
    C(:,901:1000,2) = ones(200,100)
    surf(Mem,C,'LineStyle',':', 'FaceAlpha',1) ;
    

%% Equation with Memory 


% parameters
x_length = 10 ; % space dimension
points = 100 ; % iteration number space
tfinal = 400 ; %final time
Dt = tfinal/200000 ; % time step
nt = tfinal/Dt ; % iteration number time
h = x_length/points ; % space step
afficher=true;

%Diffusion constants
phenotypes_number = 5 ; %number of phenotypes
phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
phe_diff(1) = 0.01 ;
for i=2:phenotypes_number
    phe_diff(i) = 3*phe_diff(i-1) ;
end

%Environment matrix
E = zeros(1,phenotypes_number*points) ;
for i=1:phenotypes_number
E((i-1)*points+1:(i-1)*points+20) = 0.1:-0.005:0.005 ; 
end


%Sources vectors
a = -pi/2:2*pi/(points-1):3*pi/2;% fonction de répartition de la nouriture
a= cos(a);
a=a*diag(a)*0.8+0.2;
A=zeros(1,3*phenotypes_number);
for i=0:phenotypes_number-1
A(i*points+1:i*points+points)=a;%répartition de la nouriture par espèce
end


%Mutation matrix, conservation and equipropability of mutation
Mutation=zeros(phenotypes_number);
for i=1:phenotypes_number 
    for j=1:phenotypes_number
        if i == j
            Mutation(i,j)=-1;
        else
            Mutation(i,j)=1/(phenotypes_number-1);
        end

    end

end

epsilon=0.00;

% linear Discretization matrix
M = zeros(points* phenotypes_number, points* phenotypes_number) ;

for  i=0:phenotypes_number-1
    for n=2:points-1
        k=i*points+n;
        % termes liés à la difusion
        M(k,k) = phe_diff(i+1)*(-2)/(h*h) ;
        if n<points
        M(k+1,k) = phe_diff(i+1)/(h*h) ;
        end
        if n>1
        M(k-1,k) = phe_diff(i+1)/(h*h) ;
        end
        for j=0:phenotypes_number-1
            p=j*points+n;
            % termes liés aux mutations
            M(p,k)= M(p,k)+(epsilon*Mutation(i+1,j+1));
        end
    end
    
    % termes aux limites
    M(i*points+1,i*points+1)=-phe_diff(i+1)/(h*h);
    M(i*points+2,i*points+1)= phe_diff(i+1)/(h*h);
    M((i+1)*points,(i+1)*points)=-phe_diff(i+1)/(h*h);
    M((i+1)*points-1,(i+1)*points)= +phe_diff(i+1)/(h*h);
    
    
end

% non linear Discretization matrix
M2 = zeros(points* phenotypes_number, points* phenotypes_number) ;

for i=0:phenotypes_number-1
    for n=1:points
        k=i*points+n;
        for j=0:phenotypes_number-1
            M2(j*points+n,i*points+n)=-1;  % termes liés à la consomation de ressources
            M2(k,k)=-1;  % termes liés à la consomation de ressources
        end
    end
    
end







% Evolution
Mem=zeros(100,phenotypes_number*points);

t = 0 ;
while t < nt-1
     E=E+ Dt*(E*M+E*(diag(A+E*M2)));
    if mod(t,1000)==0
        Mem(t/1000 +1,:)=E;
    end
    if mod(t,1000)==0 && afficher 
        figure(1)
        Es=0;
        for i=1:phenotypes_number
            plot((h:h:x_length), E((i-1)*points+1:(i-1)*points+points),'Color',[i/phenotypes_number,(phenotypes_number-i)/phenotypes_number,1]) ;
            hold on
            Es=Es+E((i-1)*points+1:(i-1)*points+points);
        end
        plot((h:h:x_length),Es,'r--')
        hold on
        plot((h:h:x_length),a,'g--')
        hold off
    end
    t=t+1;
end

figure(1)
Es=0;
for i=1:phenotypes_number
    plot((h:h:x_length), E((i-1)*points+1:(i-1)*points+points),'Color',[(mod(i,3)==1)*i/phenotypes_number,(mod(i,3)==2)*i/phenotypes_number,(mod(i,3)==0)*i/phenotypes_number]) ;
    hold on
    Es=Es+E((i-1)*points+1:(i-1)*points+points);
end
plot((h:h:x_length),Es,'r--')
hold on
plot((h:h:x_length),a,'g--')
hold off

figure(2)
clear C
    C(:,1:100,1) = ones(100,100)*0.5;
    C(:,101:200,2) = ones(100,100)*0.7;
    C(:,201:300,3) = ones(100,100)*0.8;
    C(:,301:400,1) = ones(100,100)*0.9;
    C(:,401:500,2) = ones(100,100);
    %C(:,501:600,3) = ones(100,100)*0.6;
    %C(:,601:700,1) = ones(100,100)*0.7;
    %C(:,701:800,2) = ones(100,100)*0.8;
    %C(:,801:900,3) = ones(100,100)*0.9;
    %C(:,901:1000,2) = ones(100,100)
    surf(Mem,C,'LineStyle','-', 'FaceAlpha',1) ;
    
%% Equation with Penalisation of slow dif


% parameters
x_length = 10 ; % space dimension
points = 100 ; % iteration number space
tfinal = 400 ; %final time
Dt = tfinal/500000 ; % time step
nt = tfinal/Dt ; % iteration number time
h = x_length/points ; % space step
afficher=true;

%Diffusion constants
phenotypes_number = 2 ; %number of phenotypes
phe_diff = zeros(phenotypes_number, 1) ; %vector of phenotypes diffusion
phe_diff(1) = 0.1 ;
phe_diff(2) = 0.9 ;

%Environment matrix
E = zeros(1,phenotypes_number*points) ;
for i=1:phenotypes_number
E((i-1)*points+1:(i-1)*points+20) = 0.1:-0.005:0.005 ; 
end


%Sources vectors
a = -pi/2:2*pi/(points-1):3*pi/2;% fonction de répartition de la nouriture
a= cos(a);
a=a*diag(a)*0.8+0.2;
A=zeros(1,points*phenotypes_number);
for i=0:phenotypes_number-1
A(i*points+1:i*points+points)=a;%répartition de la nouriture par espèce
end
A(1:points)=A(1:points)-0.05; %pénalisation de A

%Mutation matrix
Mutation = zeros(phenotypes_number, phenotypes_number) ;
Mutation(1,1) = -1;
Mutation(1,2) = 0.5 ;
Mutation(1,3) = 0.5;
Mutation(2,1) = 0.5;
Mutation(2,2) = -1;
Mutation(2,3) = 0.5 ;
Mutation(3,1) = 0.5 ;
Mutation(3,2) = 0.5 ;
Mutation(3,3) = -1 ;

epsilon=0.0;

% linear Discretization matrix
M = zeros(points* phenotypes_number, points* phenotypes_number) ;

for  i=0:phenotypes_number-1
    for n=2:points-1
        k=i*points+n;
        % termes liés à la difusion
        M(k,k) = phe_diff(i+1)*(-2)/(h*h) ;
        if n<points
        M(k+1,k) = phe_diff(i+1)/(h*h) ;
        end
        if n>1
        M(k-1,k) = phe_diff(i+1)/(h*h) ;
        end
        for j=0:phenotypes_number-1
            p=j*points+n;
            % termes liés aux mutations
            M(p,k)= M(p,k)+(epsilon*Mutation(i+1,j+1));
        end
    end
    
    % termes aux limites
    M(i*points+1,i*points+1)=-phe_diff(i+1)/(h*h);
    M(i*points+2,i*points+1)= phe_diff(i+1)/(h*h);
    M((i+1)*points,(i+1)*points)=-phe_diff(i+1)/(h*h);
    M((i+1)*points-1,(i+1)*points)= +phe_diff(i+1)/(h*h);
    
    
end

% non linear Discretization matrix
M2 = zeros(points* phenotypes_number, points* phenotypes_number) ;

for i=0:phenotypes_number-1
    for n=1:points
        k=i*points+n;
        for j=0:phenotypes_number-1
            M2(j*points+n,i*points+n)=-1;  % termes liés à la consomation de ressources
        end
    end
    
end







% Evolution
Mem=zeros(100,phenotypes_number*points);

t = 0 ;
while t < nt-1
     E=E+ Dt*(E*M+E*(diag(A+E*M2)));
    if mod(t,5000)==0
        Mem(t/5000 +1,:)=E;
    end
    if mod(t,1000)==0 && afficher 
        figure(1)
        Es=0;
        for i=1:phenotypes_number
            plot((h:h:x_length), E((i-1)*points+1:(i-1)*points+points),'Color',[i/phenotypes_number,(phenotypes_number-i)/phenotypes_number,1]) ;
            hold on
            Es=Es+E((i-1)*points+1:(i-1)*points+points);
        end
        plot((h:h:x_length),Es,'r--')
        hold on
        plot((h:h:x_length),a,'g--')
        hold off
    end
    t=t+1;
end
figure(1)
Es=0;
for i=1:phenotypes_number
    plot((h:h:x_length), E((i-1)*points+1:(i-1)*points+points),'Color',[(mod(i,3)==1)*i/phenotypes_number,(mod(i,3)==2)*i/phenotypes_number,(mod(i,3)==0)*i/phenotypes_number]) ;
    hold on
    Es=Es+E((i-1)*points+1:(i-1)*points+points);
end
plot((h:h:x_length),Es,'r--')
hold on
plot((h:h:x_length),a,'g--')
hold off

figure(2)
    clear C
    C(:,1:100,1) = ones(100,100)*0.5;
    C(:,101:200,2) = ones(100,100);
    C(:,101:200,3) = zeros(100,100);
    
    surf(Mem,C,'LineStyle','-', 'FaceAlpha',1) ;