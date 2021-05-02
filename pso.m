%%

clc; clear; close all; 
%%
% while 1
%Bu sat�rlar� kullanarak komut penceresinden input al�nabilir.   
prompt = 'Target=? '; 
target = input(prompt);

% figure(1)
% xlabel('x coordinate ');
% ylabel('y coordinate ');
% %zlabel('z coordinate '); % bunu ben ekledim �al��mad� de�il bu fonksiyon x
% %ve y al�r :) get ile denedim olmad� :) 
% grid on
% axis([-50 50 -50 50 ]);
% title('Click on a coordinate as the target')
% target = ginput(1); % Grafik �zerinden bir kordinata t�klanarak input al�nacak.

%%
%general format for anonymous function: 
%variableName = @(inputVariable) statment(function)
% sqr = @(x) x.^2;
% a = sqr(5)
% a =
%    25
distanceFunction = @(x,target) distance(x,target);  % Fonksiyon
%tan�mland�. asl�

 %%
nVar = 3; %3 boyutlu
VarSize = [1 nVar]; % Her par�ac���n x, y ve z  kordinatlar�n� tan�mlamak i�in
VarMin = -10; % Her par�ac���n ilk konumu rastegele se�ilirken s�n�rland�rma i�in min. ve max. de�erleri
VarMax = 10;
MaxIt = 100; % Her par�ac���n hedefe gitmesi i�in tekrarlanma say�s�
nPop = 200; % par�ac�k say�s�
w = 1; % inertia coefficient
c1 = 2; % cognitive coefficient 
c2 = 2; % social coefficient
Positions = zeros(nPop,nVar); % (par�ac�k say�s�,pozisyon de�er say�s�)
BestPos = zeros(MaxIt, 3); 
% Her par�ac��a pozisyon, h�z, mesafe, en iyi pozisyon ve en iyi mesafe
% �eklinde ba�ta bo� olmak �zere tan�mlamalar yap�ld�.
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.distance = [];
empty_particle.Best.Position = [];
empty_particle.Best.distance = [];
particle = repmat(empty_particle,nPop,1);
GlobalBest.distance = Inf;

 for i=1:nPop % Her par�ac�k i�in tekrar et 
     
 particle(i).Position = unifrnd(VarMin, VarMax, VarSize); % �lk konumu rastgele belirle
 particle(i).Velocity = zeros(VarSize); % �lk h�z� s�f�r olarak belirle
 particle(i).distance = distanceFunction(particle(i).Position,target); % �lk mesafe de�erini belirle
 particle(i).Best.Position = particle(i).Position; % Her poziyonu en iyi pozisyon olarak belirle
 particle(i).Best.distance = particle(i).distance; % Her mesafe de�erini en iyi mesafe olarak belirle
 
    if particle(i).Best.distance < GlobalBest.distance % �lk pozisyonlar i�in global en iyi de�erleri belirle (global en iyi pozisyon ve global en iyi mesafe)
       GlobalBest = particle(i).Best; 
    end   
 end
 
 for it=1:MaxIt % Tekrarlanma say�s� kadar tekrarla
     
    for i=1:nPop % Her par�ac�k i�in tekrar et 
        
        for a=1:nPop  
        Positions(a, :) = [particle(a).Position]; % Her par�ac���n pozisyon de�erlerini 'Positions' matrisine s�rayla ata.
        end
    
  % PSO sonucu belirlenen bir sonraki h�z ve pozisyon bulma form�lleri  
  %g�ncellemeler yap�l�yor
    %bir sonraki h�z� bulma 
    %eski h�z vekt�r�+ ki�isel en iyi ve global en iyi de�erleri kullan�l�yor ve yeni h�z vekt�r�
    %olu�turur
    particle(i).Velocity = w*particle(i).Velocity...
    + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position)...
    + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position); 

    %bir sonraki pozisyonu bulma - o anki konumu+ h�z vekt�r�
    particle(i).Position = particle(i).Position + particle(i).Velocity;
    
    %bir sonraki mesafeyi bulma
    particle(i).distance = distanceFunction(particle(i).Position,target);
    
         %tecr�belerim �ncekilerden iyi mi?
         if particle(i).distance < particle(i).Best.distance % E�er yeni mesafe, bir �nceki en iyi mesafesinden daha k���k ise;
         particle(i).Best.Position = particle(i).Position;   % Yeni pozisyon, en iyi pozisyonu olsun. 
         particle(i).Best.distance = particle(i).distance;   % Yeni mesafe, en iyi mesafesi olsun.
            %ondan sonra
            if particle(i).Best.distance < GlobalBest.distance % E�er yeni en iyi mesafe, global en iyi mesafeden k���k ise;
            GlobalBest = particle(i).Best;                     % O par�ac���n en iyi de�erleri global olarak atans�n.
            end
         end
    end
    BestPos(it,:) = GlobalBest.Position; % Her tekrarlanma a�amas�n�n global en iyi pozisyon de�erini 'BestPos' matrisine kaydet
    
    w = w - w*it/MaxIt; % inertia coefficient de�erini her tekrarlanmada azalt 
    
    
%%% HER TEKRARLANMADA PAR�ACIKLARIN KONUMLARINI GRAF�KTE G�STER %%%
    figure(1);
    plot3(Positions(:, 1), Positions(:, 2),Positions(:, 3), 'x');
    xlabel('x coordinate ');
    ylabel('y coordinate ');
    zlabel('z coordinate ');
    grid on
    axis([-50 50 -50 50 -50 50]);
%     surf(); %denedim olmad�    :)
    title(['Iteration = ' num2str(it) ',    Target = ' num2str(target) ',    Best Position = ' num2str(BestPos(it,:))])
    
    
    
    
    %%% BU GRAF�KLE DE TEK B�R PAR�ACI�IN HEDEFE VARANA KADAR HAREKET� G�ZLEMLENEB�L�R
%     figure(2);
%     plot(Positions(1, 1), Positions(1, 2), 'x')
%     xlabel('x coordinate ');
%     ylabel('y coordinate ');
%     grid on
%     axis([-50 50 -50 50])
    
    pause(.1) % Her tekrarda par�ac�klar�n hareketini g�zlemleyebilmek i�in, her tekrarda .1 saniye bekle
    
 end
 

 
% end







