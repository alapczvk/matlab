clear all; close all;
%zadanie 41 i 42
A=[1,2,1,3;4,1,5,6;3,2,1,2;3,4,5,3]; %macierz 4x4
B=[1,2,3,4,5;2,3,4,5,1;2,4,5,6,1;2,3,4,2,1;3,5,2,1,4]; %macierz 5x5
AC=[1,4,3;4,2,1;3,1,3]; %macierz ymetryczna do dekompozycji choleskiego
b=1; % zmiana choice 
myLu(A,b);
[L1,U1]=lu(A),
spr1=L1*U1,

myLu(AC,b); %dekompozycja choleskiego
function [L,U] = myLu(A,b) %funkcja wypluwa macierze L i U, a przyjmuje macierz A do rozkladu i b
A=[1,2,1,3;4,1,5,6;3,2,1,2;3,4,5,3]; %macierz 4x4
%A=[1,2,3,4,5;2,3,4,5,1;2,4,5,6,1;2,3,4,2,1;3,5,2,1,4]; %macierz 5x5
%A=[1,2,3;4,5,6;3,2,1];
%b=0;
[N,N] = size(A); %3x3- A=macierz kwadratowa

if (b==0) % prosciej, wolniej, algorytm doolittle- przyrost o nowa zmienna   w kazdym elemencie kolumny oblizanych macierzy
  L = eye(N); U = zeros(N,N); %L= macierz diagonalna o wymiarze N, U=macierz  wypelniona zerami o wymiarze N
  for i = 1:N
      for j=i:N %iterujemy po wszystkich elementach
          U(i,j) = A(i,j) - L(i,1:i-1)*U(1:i-1,j);
      end
      for j=i+1:N
          L(j,i) = 1/U(i,i) * ( A(j,i) - L(j,1:i-1)*U(1:i-1,i) );
      end
  end
  %dodanie mozliwosci obliczen
elseif(b==1)
    L=eye(N); U=zeros(N,N);
    for i=1:N
        for j=1:i-1
            L(i,j)=1/U(j,j)*(A(i,j)-L(i,1:j-1)*U(1:j-1,j));
        end
        for j=i:N
            U(i,j)=A(i,j)-L(i,1:i-1)*U(1:i-1,j);
        end    
    end
 elseif(b==2) %dekompozycja Choleskiego
   L = eye(N);
   for j=1:N
       value=0;
       for k=1:j-1
           value = value + L(j,k)*L(j,k);
       end
       L(j,j) = sqrt(A(j,j) - value);
       for i=j+1:N
           value = 0;
           for k=1:j-1
               value = value + L(i,k)*L(j,k);
           end
           L(i,j) = (1/L(j,j)) * (A(i,j) - value);
       end
       
   end
   U = L.'; %U to macierz trojkatna transponowana
else % trudniej, szybciej ----------------------------------------
U=A; L=eye(N);
  for i=1:N-1
     for j=i+1:N
        L(j,i) = U(j,i) / U(i,i);
        U(j,i:N) = U(j,i:N) - L(j,i)*U(i,i:N);
     end
  end
  
end
disp("U:");
disp(U);
disp("L:");
disp(L);
spr=L*U,
%[L1,U1]=lu(A),
%chol(AC);
end

