%% Particle Swarm Optimization
% Jonas M�ller Gon�alves
clear all;clc;
%% Inicializa��o das Vari�veis
% f = @(x,y) -20.*exp(-0.2.*sqrt(0.5.*(x.^2 + y.^2))) - exp(0.5.*(cos(2.*pi.*x) + cos(2.*pi.*y))) + exp(1) + 20;  %Fun��o de Ackley
% f = @(x,y) 0.5 + ((sin(x.^2-y.^2)).^2-0.5)./((1+0.001.*(x.^2+y.^2)).^2);
f = @(x) 100*(x(2)-x(1).^2).^2+(1-x(1)).^2;
% f = @(x,y) (1-x).^2+100.*(y-x.^2).^2;
% f = @(x,y) 0.1.*x.*y;
% f = @(x,y) (x-2.*y-7).^2+(2.*x+y-5).^2;
% f = @(x,y) (1.5 -x +x.*y).^2 + (2.25 -x + x.*y.^2).^2 + (2.625 -x + x.*y.^3).^2;
% restr = @(x,y) x.^2+y.^2-2
% restr = @(x,y) (1+0.2.*cos(8.*atan(x./y))).^2

rng(5);
population = 20; npart=population;   % Popula��o total de Part�culas 
variables = 2; nvar=variables;     % Vari�veis de Projeto
lambda1 = 2;    % Lambda 1
lambda2 = 2;    % Lambda 2
omega = 0.5;    % Propriedade de In�rcia
dominio=[-4.5 4.5;-4.5 4.5]; % Dom�nio das nvar variaveis
niter = 2000;    % N�mero de Itera��es

%% Inicializa��o dos Vetores de Projeto
% Distribuindo as Part�culas
for icol = 1:nvar
    for ivar = 1:npart
        position(ivar,icol) = dominio(icol,1)+rand*(dominio(icol,1)-dominio(icol,2));
    end
end

% Definindo as Velocidades Iniciais
velocity = zeros(npart,nvar);   % inicializa as velocidades

%% C�lculo Inicial da Melhor Part�cula

for ires=1:size(position,1)
    resultados(ires,1) = f(position(ires,:));    % Define os valores de f(x,y)
end

[resg,ord]=sort(resultados);    % Ordena os valores do menor para o maior

xlbest = zeros(npart,nvar); % Inicia Vetor de melhor local

for ixgbest=1:nvar
    xgbest(ixgbest) = [position(ord(1,1),ixgbest)];   % Define o melhor global como o primeiro do vetor ordenado
end

melhorglobal=xgbest;    % Cria um valor armazenado do melhor global
ordglobal = ord(1,1);   % Replica qual a linha do melhor global

% Vari�veis Auxiliares
parada=0; iter=0; pos_old=position; v_old=velocity; cont1=1; cont2=1;


%% Otimiza��o
while parada == 0
    iter = iter+1; % Atualiza o Contador de Itera��es
    r1=rand; r2=rand; % R's aleat�rios
    
    % C�lculo das Novas Velocidades e Posi��es de cada part�cula 
    for i=1:npart
        for j=1:nvar
            v(i,j) = omega*v_old(i)+lambda1*r1*(xlbest(i,j) - pos_old(i,j)) + lambda2*r2*(xgbest(j) - pos_old(i,j));
            pos(i,j) = pos_old(i,j) + v(i,j);
        end
    end
    
    % Verifica��o do Dom�nio das Novas Part�culas
    for ia=1:npart      
        if pos(ia,1) > dominiox(2) || pos(ia,1) >  dominiox(2)
            passoux(cont1,1) = pos(ia,1);   % Part�culas que passaram do dom�nio x
            pos(ia,1) = pos(ia,1) - (pos(ia,1)-1);  % Induz o retorno da part�cula que passou
            cont1=cont1+1;  % Atualiza o contador x      
        elseif pos(ia,2) <  dominioy(1) || pos(ia,2) <  dominioy(1)
            passouy(cont2,1) = pos(ia,2);   % Part�culas que passaram do dom�nio y
            pos(ia,2) = pos(ia,2) - (pos(ia,2)-1);  % Induz o retorno da part�cula que passou
            cont2=cont2+1;  % Atualiza o contador y
        end 
    end
    
    % Verifica��o da Restri��o 
%     for ia=1:npart
%         if pos(ia,1)^2+pos(ia,2)^2  > (pos(ia,1)^2-pos(ia,2)^2)% (1+0.2.*cos(8.*atan(pos(ia,1)./pos(ia,2)))).^2
%             pos(ia,1) = exp(pos(ia,1));
%             pos(ia,2) = exp(pos(ia,2));
%         end
%     end
    
    % C�lculo da Nova Melhor Part�cula Local
    resultados = f(pos(:,1),pos(:,2));
    [resl,ord]=sort(resultados);
    xlbest = [pos(ord,1) pos(ord,2)];
    melhorlocal = [xlbest(1,1) xlbest(1,2)];
    ordlocal = ord;
    
    v_old = v;  % Atribui a velocidade antiga como a nova velocidade calculada
    pos_old = pos;  % Atribui a posi��o antiga como a nova posi��o calculada (j� com o dom�nio verificado)
    
    % Atribui�ao da Nova Melhor Part�cula Global (caso exista)
    if resg(1,1) > resl(1,1)
        resg=resl;
        melhorglobal = resl;
        xgbest=[xlbest(1,1) xlbest(1,2)];
    end
    melhores(iter,:) = xgbest;
    
    % Condi��o de Parada
    if iter==niter
        parada=1;
    end
end
figure(1)
fsurf(f);hold on;
scatter3(melhores(:,1),melhores(:,2),f(melhores(:,1),melhores(:,2)));hold on;
% fsurf(restr,[-1.5 1.5]);

figure(2);
scatter(1:niter,f(melhores(:,1),melhores(:,2)));hold on;