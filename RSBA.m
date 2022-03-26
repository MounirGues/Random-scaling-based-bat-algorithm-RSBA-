clear


disp('********* optimization has started ************')

%Greenhouse climate datasets:::
load('1-  TinM days23 SIMUL Chi=5.2 [5.8 5.3 1.7 2.2].mat')
load('DayN-2-3.mat')
% load('DayN-3-4.mat')
% load('DayN-4-5.mat')
load('DB5Days.mat')
DayNX=DayN23;

%Bat algorithm control parameters:::::
n=100;         % Population size, typically 10 to 40
% A=0.25;      % Loudness  (constant or decreasing)
% r=0.5;      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
Qmin=0;         % Frequency minimum %frequancy tuning
Qmax=1.5;         % Frequency maximum
% Max_iterations=2000;
% Iteration parameters
tol=10^-8;
N_iter=1;       % Total number of function evaluations
d=4;            % Number of dimensions
% initial array
Q=zeros(n,1);% frequency
v=zeros(n,d);% velocities
S=zeros(n,d);% positions
Sol=zeros(n,d);%best positions
Fitness=zeros(n,1);%fitness
Convergence_curve=zeros(n,1);

% initialize loudness and pulse rate
A=ones(n,1);            %loudness for each BAT
r=rand(n,1);            %pulse emission rate for each BAT
alpha=0.9;              %constant for loudness update
gamma=0.9;              %constant for emission rate update
r0=0.3;

%Online plotting parameters::::
p=20;
po=-1;
pp=1;
Particles=cell(1,10);
RE=cell(1,10);

%Upper and lower search boundaries - Search space for each unknown parameter::::
z=100;
Lb=[  0     0    0   0];
Ub=[  z     z    z   z];


% initialize the population/solution
for i=1:n
    Sol(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
    %%%%%%% Test%%%%%%%%
    Fitness(i)=FunDB(Sol(i,:),DayNX,DB5Days,TinM);
    %%%%%%%%%%%%%%%%%%%%
end

% find the current best
[fmin,I]=min(Fitness);
best =Sol(I,:);

% start the iterations
while (N_iter<501)

% loop over all bats solutions
for i=1:n    
    Q(i)=Qmin+(Qmax-Qmin)*rand;
    v(i,:)=v(i,:)+(Sol(i,:)- best)*Q(i); % Vi the new velocities

    for j=1:length(Lb)
         if v(i,j)>Ub(j)
           v(i,j)=Lb(j)+(Ub(j)-Lb(j)).*rand;
         elseif v(i,j)<Lb(j)
           v(i,j)=Lb(j)+(Ub(j)-Lb(j)).*rand;
         end 
    end
    
        %%%%%%%%%%%%%%%%%%%%%%
        S(i,:)=Sol(i,:)+v(i,:);  % Xi the new positions (solutions)
        %%%%%%%%%%%%%%%%%%%%%%
        
% Keep the new values inside the search space::::
     for j=1:length(Lb) 
         if S(i,j)>Ub(j)
            S(i,j)=Lb(j)+(Ub(j)-Lb(j)).*rand.*0.1; 
         elseif S(i,j)<Lb(j)
            S(i,j)=Lb(j)+(Ub(j)-Lb(j)).*rand.*0.1;
         end
     end

% If bat's pulse emission increased, a solution exists  
    if rand>r
        eps=0.0000001+(0.01-(0.0000001))*rand; %Parameter Randomization using a very small low and upper boundaries (it depends on the nature of the unknown parameters) for the scaling parameter (eps)
        S(i,:)=best+eps*randn(1,d);
    end
%     evaluate new solutions
          %%%%%%%%%%%%%%%%
          Fnew=FunDB(S(i,:),DayNX,DB5Days,TinM);  %FunDB is the cost or the objective function
          %%%%%%%%%%%%%%%% 
    % Lodness decreases = Solution existance confirmation 
     if (Fnew<=Fitness(i)) && (rand<A(i))
                    Sol(i,:)=S(i,:);
                    Fitness(i)=Fnew;
                    A(i)=alpha*A(i);
                    r(i)=r0*(1-exp(-gamma*N_iter));
     end
  % Update the current best solution
          if Fnew<=fmin
                best=S(i,:);
                fmin=Fnew;
          end
  end
      
%% Online plotting section::::::
uuu=[5.8 5.3 1.7 2.2];
batsPP(1,:)=best;
Particles{1}(pp,:)=best;
RE{1}(pp,:)=(abs(batsPP-uuu)./uuu)*100;
% Online_Plot_BA;
pp=pp+1;

Convergence_curve(N_iter)=fmin;
disp(['Number of evaluations: ',num2str(N_iter),' best= ',num2str(best)]);
disp(['fmin=',num2str(fmin)]);

%Increamental iterations:
N_iter=N_iter+1;

end

disp('*************Convergance is done*************');
plot(Convergence_curve)
title('Convergence Curve')
xlabel('Iterations')
ylabel('Fitness Value')
