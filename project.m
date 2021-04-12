clc
clear
load('iddata-11.mat');
%% generare puteri
counter=0;
for m=1:3   % gradul polinomului
    for na=1:3
 counter=counter+1;
 nk=0; %intarzierea
 nb=na;
 n=na+nb; % numarul de parametrii           
a=zeros(1,n); % matricea cu ajutorul careia cream puterile polinomului
matricea_puterilor=zeros(1,n); % matricea puterilor
while (a(n) ~= m) % conditia de oprire a algoritmului de generare a puterilor expresiei x1^p1*x2^p2*...xn^pn, deoarece cand ultimul element din sir ajunge la gradul maxim, nu mai are rost sa generam numere in sirul 00000m ->10000m deoarece depaseste gradul polinomului
    a(1)=a(1)+1; % crestem cu o unitate primul element din sir a(1)a(2)a(3)...a(n)
    j=1; % cu ajutorul contorului j parcurgem sirul de la a(1) la a(n).
    while(a(j)>=m) % conditia de repetare a buclei in cazul in care in sir sunt elemente de valoare maxima(gradul polinomului) dupa ce am parcurs o data sirul        
        if(a(j)>=m) % conditie de verificare daca un element din sir a depasit valoarea maxima, in caz afirmativ 
            for i=1:n % luam un for care sa parcurga sirul de la inceput pana la final si verificam urmatoarele:
                if(a(i)>m && a(n)<m) % daca avem un element care depaseste gradul maxim il facem 0 si il crestem pe urmatorul cu o unitate
                    a(i)=0;          % acest lucru se poate face doar daca ultimul element din sir nu are valoarea maxima. In cazul in care ultimul element din sir are valoarea maxima am depasi lungimea sirului.
                    a(i+1)=a(i+1)+1;
                end
            end
        end
         j=j+1;
         if(j==n+1) % conditia de resetare a contorului dupa ce am terminat de verificat sirul
             j=1;
         end
    end
    if(sum(a)<=m) % in cazul in care suma elementelor din sir nu depaseste gradul polinomului, sirul generat este considerat sirul de puteri pentru expresia x1^p1*x2^p2*...xn^pn
        matricea_puterilor=cat(1,matricea_puterilor,a);
    end
end

%% generare regresori
N=length(id.u);
phi=0; % initializam matricele phi, phival, x, xval = 0 
x=0; % matricele x si xval sunt folosite pentru a ajuta la calculul elementelor matricei phi si phival

xval=0;
phival=0;

for k=1:N % parcurgem elementele din setul de date
    if k<=length(id.u) % conditie folosita pentru a optimiza codul, folosim un singur for pentru calulul ambelor valori phi si phival
        for i=1:na % cu ajutorul unui for punem in matricea x elemente din datele de identificare
            if(k-i-nk>0) % verificam daca indexul elementelor exista in setul de date. In cazul in care exista:
                x(i)=id.y(k-i); % atribuim elementelor pana la pozitia na valorile de pe iesire
                x(na+i)=id.u(k-nk-i); % de la pozitia na pana la pozitia nb=na+na atribuim elmentelor valoarea de pe intrare
            else
                x(i)=0; % in cazul in care avem indici care nu se regasec in setul de date punem valoarea 0.
                x(na+i)=0;
            end
        end    
        for i=1:length(matricea_puterilor)
            phi(k,i)=prod(x.^matricea_puterilor(i,1:end)); % calculam matricea linie cu linie si pe linie element cu element astfel: aplicam pentru fiecare element expresia x1^p1*x2^p2*...*xn^pn, unde x este vectorul format anterior si p1..pn sunt elementele din matricea puterilor alcatuita in algoritmul de generare de mai sus. Fiecarui element de pe linie din matricea phi ii corespunde o linie din matricea puterilor
        end
    end
    
    if k<=length(val.u) % urmatoarea secventa de cod se repeta pentru date de validare
        for i=1:na
            if(k-i-nk>0)
                xval(i)=val.y(k-i);
                xval(na+i)=val.u(k-i-nk);
            else
                xval(i)=0;
                xval(na+i)=0;
            end
        end
   
        for i=1:length(matricea_puterilor)            
            phival(k,i)=prod(xval.^matricea_puterilor(i,1:end));                    
        end
    end
end

%% aflare parametrii

theta =phi\id.y; % aflam parametrii polinomului prin impartire la stanga

%% validare

ypredid=phi*theta; % calculam y predictie pentru setul de identificare
ypredval=phival*theta; % calculam y predictie pentru setul de validare

if(m==3 && na==2) % salvam valoarea optima a aproximatorului
    ypredvaloptim=ypredval;
end

if(m==3 && na==3) % salvam valoarea optima a aproximatorului
    ypredidoptim=ypredid;
end

   
MSE1valpred(counter)=1/length(val.y)*sum((ypredval-val.y).^2); % calculam eroarea medie patratica pentru setul de validare pentru y predictie
MSE1idpred(counter)=1/length(id.y)*sum((ypredid-id.y).^2); % calculam eroarea medie patratica pentru setul de identificare pentru y predictie

       xvalsim=0;
       phivalsim=0;
       xidsim=0;
       phiidsim=0;
       
          
      for k2=1:length(id.u) % aplicam algoritmul de mai sus pentru sedul de validare pentru a calcula y simulat
          if(k2<=length(val.u))
            for i=1:na % diferenta intre algoritmul descris mai sus este ca la y simulat ne folosim de valorile anterioare simulate si nu de valorile de la iesirea sistemului real
                if(k2-i-nk>0)
                    xvalsim(i)=ysimval(k2-i);
                    xvalsim(na+i)=val.u(k2-nk-i);
                else
                    xvalsim(i)=0;
                    xvalsim(na+i)=0;
                end
            end
            for i=1:length(matricea_puterilor)
                phivalsim(k2,i)=prod(xvalsim.^matricea_puterilor(i,1:end));
            end
            ysimval(k2)=phivalsim(k2,1:end)*theta;
          end
          
          
            for i=1:na 
                if(k2-i-nk>0)
                    xidsim(i)=ysimid(k2-i);
                    xidsim(na+i)=id.u(k2-nk-i);
                else
                    xidsim(i)=0;
                    xidsim(na+i)=0;
                end
            end
            for i=1:length(matricea_puterilor)
                phiidsim(k2,i)=prod(xidsim.^matricea_puterilor(i,1:end));
                
            end
            ysimid(k2)=phiidsim(k2,1:end)*theta;
          
          
      end
      
      MSE2valsim(counter)=1/length(val.y)*sum((ysimval'-val.y).^2); % calculam eroarea medie patratica pentru y simulat pe datele de validare
      MSE2idsim(counter)=1/length(id.y)*sum((ysimid'-id.y).^2); % calculam eroarea medie patratica pentru y simulat pe datele de identificare
      
if(m==2 && na==1) % salvam valoarea optima pentru y simulat
    ysimoptim=ysimval;
end
      
if(m==3 && na==3) % salvam valoarea optima a aproximatorului
    ysimidoptim=ysimid;
end   
   
   
   [m na] % trebuie sters
      
    end
end
%%
clc
grad_polinom=0;
numar_parametrii=0;
for i=1:m
    grad_polinom=cat(1,grad_polinom,i*ones(na,1));
    numar_parametrii=cat(1,numar_parametrii,(1:na)');
    
%     grad_polinom=[1*ones(na,1);2*ones(na,1);3*ones(na,1);4*ones(na,1);5*ones(na,1)]; % algoritmul de mai jos este pentru a crea un tabel cu valorile MSE1 si MSE2 in functie de m si na=nb
%     numar_parametrii=[(1:na)'; (1:na)'; (1:na)'; (1:na)'; (1:na)']; % algoritmul functioneaza doar in cazul in care gradul polinomului nu este mai mare decat 3
end

grad_polinom=grad_polinom(2:end,1);
numar_parametrii=numar_parametrii(2:end,1);

MSE1valpred=MSE1valpred';
table_val_pred=table(grad_polinom,numar_parametrii,MSE1valpred)

MSE2valsim=MSE2valsim';
table_val_sim=table(grad_polinom,numar_parametrii,MSE2valsim)

MSE1idpred=MSE1idpred';
table_id_pred=table(grad_polinom,numar_parametrii,MSE1idpred)

MSE2idsim=MSE2idsim';
table_id_sim=table(grad_polinom,numar_parametrii,MSE2idsim)

figure   % afisarea sub forma de grafic: y predictie optim, y simulat optim, y real
   plot(ypredvaloptim,'b')
   hold on
   plot(val.y,'k')
   plot(ysimoptim,'r')
   legend('Y predictie','Y real', 'Y simulat')
   title('Y predictie optim/Y real/Y simulat optim setul de validare')
   
figure   % afisarea sub forma de grafic: y predictie optim, y simulat optim, y real
   plot(ypredidoptim,'b')
   hold on
   plot(id.y,'k')
   plot(ysimidoptim,'r')
   legend('Y predictie','Y real', 'Ysimulat')
   title('Y predictie optim/Y real/Y simulat optim setul de identificare')