function [A]=OMPerrn(D,X,errorGoal,T0m)
% function [A]=OMPerrn(D,X,errorGoal,T0m)
% 
% Modified by Ariel Habshush, September 2014.
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  errorGoal - the maximal allowed representation error for
%                  each signal component.
%                  T0m - maximum sparsity (number of non-zeros) allowed.
% output arguments: A - sparse coefficient matrix.
%=============================================
[~,P]=size(X);
[~,K]=size(D);
% Etot = errorGoal;
E2 = errorGoal;
% E2 = errorGoal^2;
maxNumCoef = T0m;
% A = sparse(K,P);
% A = spalloc(K,P,maxNumCoef*P);
A = zeros(K,P);
for k=1:1:P,
    %a=[];
    x=X(:,k);
    residual=x;
	indx = [];
	a = [];
	currResNorm2 = sum((abs(residual)).^2);
%     currResSum = sum( abs(residual) );
	j = 0;
    while currResNorm2>E2 && j < maxNumCoef,
%     while currResSum > Etot && j < maxNumCoef,
		j = j+1;
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-(D(:,indx(1:j))*a);
		currResNorm2 = sum((abs(residual)).^2);
    end;
   if (~isempty(indx))
       A(indx,k)=a;
   end
end;
A = sparse(A);
return;
