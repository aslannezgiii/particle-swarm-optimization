%%

clc; clear; close all; 
%%
% while 1
%Bu satýrlarý kullanarak komut penceresinden input alýnabilir.   
prompt = 'Target=? '; 
target = input(prompt);

% figure(1)
% xlabel('x coordinate ');
% ylabel('y coordinate ');
% %zlabel('z coordinate '); % bunu ben ekledim çalýþmadý deðil bu fonksiyon x
% %ve y alýr :) get ile denedim olmadý :) 
% grid on
% axis([-50 50 -50 50 ]);
% title('Click on a coordinate as the target')
% target = ginput(1); % Grafik üzerinden bir kordinata týklanarak input alýnacak.

%%
%general format for anonymous function: 
%variableName = @(inputVariable) statment(function)
% sqr = @(x) x.^2;
% a = sqr(5)
% a =
%    25
distanceFunction = @(x,target) distance(x,target);  % Fonksiyon
%tanýmlandý. aslý

 %%
nVar = 3; %3 boyutlu
VarSize = [1 nVar]; % Her parçacýðýn x, y ve z  kordinatlarýný tanýmlamak için
VarMin = -10; % Her parçacýðýn ilk konumu rastegele seçilirken sýnýrlandýrma için min. ve max. deðerleri
VarMax = 10;
MaxIt = 100; % Her parçacýðýn hedefe gitmesi için tekrarlanma sayýsý
nPop = 200; % parçacýk sayýsý
w = 1; % inertia coefficient
c1 = 2; % cognitive coefficient 
c2 = 2; % social coefficient
Positions = zeros(nPop,nVar); % (parçacýk sayýsý,pozisyon deðer sayýsý)
BestPos = zeros(MaxIt, 3); 
% Her parçacýða pozisyon, hýz, mesafe, en iyi pozisyon ve en iyi mesafe
% þeklinde baþta boþ olmak üzere tanýmlamalar yapýldý.
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.distance = [];
empty_particle.Best.Position = [];
empty_particle.Best.distance = [];
particle = repmat(empty_particle,nPop,1);
GlobalBest.distance = Inf;

 for i=1:nPop % Her parçacýk için tekrar et 
     
 particle(i).Position = unifrnd(VarMin, VarMax, VarSize); % Ýlk konumu rastgele belirle
 particle(i).Velocity = zeros(VarSize); % Ýlk hýzý sýfýr olarak belirle
 particle(i).distance = distanceFunction(particle(i).Position,target); % Ýlk mesafe deðerini belirle
 particle(i).Best.Position = particle(i).Position; % Her poziyonu en iyi pozisyon olarak belirle
 particle(i).Best.distance = particle(i).distance; % Her mesafe deðerini en iyi mesafe olarak belirle
 
    if particle(i).Best.distance < GlobalBest.distance % Ýlk pozisyonlar için global en iyi deðerleri belirle (global en iyi pozisyon ve global en iyi mesafe)
       GlobalBest = particle(i).Best; 
    end   
 end
 
 for it=1:MaxIt % Tekrarlanma sayýsý kadar tekrarla
     
    for i=1:nPop % Her parçacýk için tekrar et 
        
        for a=1:nPop  
        Positions(a, :) = [particle(a).Position]; % Her parçacýðýn pozisyon deðerlerini 'Positions' matrisine sýrayla ata.
        end
    
  % PSO sonucu belirlenen bir sonraki hýz ve pozisyon bulma formülleri  
  %güncellemeler yapýlýyor
    %bir sonraki hýzý bulma 
    %eski hýz vektörü+ kiþisel en iyi ve global en iyi deðerleri kullanýlýyor ve yeni hýz vektörü
    %oluþturur
    particle(i).Velocity = w*particle(i).Velocity...
    + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position)...
    + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position); 

    %bir sonraki pozisyonu bulma - o anki konumu+ hýz vektörü
    particle(i).Position = particle(i).Position + particle(i).Velocity;
    
    %bir sonraki mesafeyi bulma
    particle(i).distance = distanceFunction(particle(i).Position,target);
    
         %tecrübelerim öncekilerden iyi mi?
         if particle(i).distance < particle(i).Best.distance % Eðer yeni mesafe, bir önceki en iyi mesafesinden daha küçük ise;
         particle(i).Best.Position = particle(i).Position;   % Yeni pozisyon, en iyi pozisyonu olsun. 
         particle(i).Best.distance = particle(i).distance;   % Yeni mesafe, en iyi mesafesi olsun.
            %ondan sonra
            if particle(i).Best.distance < GlobalBest.distance % Eðer yeni en iyi mesafe, global en iyi mesafeden küçük ise;
            GlobalBest = particle(i).Best;                     % O parçacýðýn en iyi deðerleri global olarak atansýn.
            end
         end
    end
    BestPos(it,:) = GlobalBest.Position; % Her tekrarlanma aþamasýnýn global en iyi pozisyon deðerini 'BestPos' matrisine kaydet
    
    w = w - w*it/MaxIt; % inertia coefficient deðerini her tekrarlanmada azalt 
    
    
%%% HER TEKRARLANMADA PARÇACIKLARIN KONUMLARINI GRAFÝKTE GÖSTER %%%
    figure(1);
    plot3(Positions(:, 1), Positions(:, 2),Positions(:, 3), 'x');
    xlabel('x coordinate ');
    ylabel('y coordinate ');
    zlabel('z coordinate ');
    grid on
    axis([-50 50 -50 50 -50 50]);
%     surf(); %denedim olmadý    :)
    title(['Iteration = ' num2str(it) ',    Target = ' num2str(target) ',    Best Position = ' num2str(BestPos(it,:))])
    
    
    
    
    %%% BU GRAFÝKLE DE TEK BÝR PARÇACIÐIN HEDEFE VARANA KADAR HAREKETÝ GÖZLEMLENEBÝLÝR
%     figure(2);
%     plot(Positions(1, 1), Positions(1, 2), 'x')
%     xlabel('x coordinate ');
%     ylabel('y coordinate ');
%     grid on
%     axis([-50 50 -50 50])
    
    pause(.1) % Her tekrarda parçacýklarýn hareketini gözlemleyebilmek için, her tekrarda .1 saniye bekle
    
 end
 

 
% end







