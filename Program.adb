Program
With Text IO,Ada.Numerics.Generic Elementary Functions;
With Ada.Numerics;
use Text IO,Ada.Numerics;
PROCEDURE eutecada IS
package real io is new FLOAT IO(float);
use real io;
package integer io is new INTEGER IO(integer);
use integer io;
package elem fct is new Ada.Numerics.Generic Elementary Functions (Float);
use elem fct;
Relax: CONSTANT :=0.1; — Under - relaxation
Rg: CONSTANT :=8.314510; — R
MaxConst : CONSTANT :=20; — Max number of components
Precision : CONSTANT := 0.00001; — Accuracy
To : CONSTANT := 273.15; — Temperature reference
N const: integer;
CaractH : ARRAY (1..MaxConst ) OF float;
CaractT : ARRAY (1..MaxConst ) OF float;
i,j:integer;
Jc : ARRAY (1..MaxConst,1..MaxConst ) OF Float ;
Dx,X,F : ARRAY (1..MaxConst ) OF Float ;
coef,g,residu,Tp : Float ;
– Jc jacobian matrix inversion then multiply by F vector
PROCEDURE SolvJc
IS
i,j,k,l,n:integer;
x: Float ;
BEGIN
FOR I IN 1.. N Const LOOP
BEGIN
X:=Jc(i,i);
Jc(i,i):=1.0;
FOR J IN 1.. N Const LOOP
Jc(I,J):=Jc(I,J)/X;
END LOOP;
FOR K IN 1.. N Const LOOP
IF K/=I THEN
BEGIN
X:=Jc(K,I);
Jc(K,I):=0.0;
684 L. Brunet, J. Caillard & P. Andr´e
FOR L IN 1.. N Const LOOP
Jc(K,L):=Jc(K,L)-X*Jc(I,L);
END LOOP;
END;
END IF;
END LOOP;
END;
END LOOP;
FOR i IN 1.. N Const LOOP
BEGIN
Dx(i):=0.0;
FOR j IN 1.. N const LOOP
Dx(i):=Dx(i)+Jc(j,i)*F(j);
END LOOP;
END;
END LOOP;
END;
BEGIN
Put(”.....................................................................”);New Line;
Put(”.....................................................................”);New Line;
Put(”.....................................................................”);New Line;
Put(” ”);New Line;
Put(” N Components Eutectic ”);New Line;
Put(” ”);New Line;
Put(” L.Brunet 2002 ”);New Line;
Put(”.....................................................................”);New Line;
Put(”.....................................................................”);New Line;
Put(”.....................................................................”);New Line;
Put(”N components= ”);
Get(N const);
FOR i IN 1.. N const LOOP
BEGIN
Put(”Enter Hfus and Tfus for ”);Put(i); New Line;
Put(” H(kJ/mol)=”);
Get(CaractH(i));
Put(” T(K)=”);
Get(CaractT(i));
END;
END LOOP;
Put(”Running...”);New Line;
— Create jacobian matrix
FOR i IN 1.. N Const LOOP
BEGIN
Thermodynamic Calculation of n-Component Eutectic Mixtures 685
FOR j IN 1.. N const LOOP
Jc(i,j):=0.0;
END LOOP;
Dx(i):=0.0;
F(i):=0.0;
X(i):=1.0/float(N Const);
END;
END LOOP;
X(N Const):=CaractT(1);
LOOP
FOR i IN 1.. N const LOOP
BEGIN
IF i¡N Const THEN
Jc(i,i):=1.0/X(i)
;
ELSE
BEGIN
coef:=1.0;
FOR j IN 1.. N Const-1 LOOP
coef:=coef-X(j);
END LOOP;
FOR j IN 1.. N Const-1 LOOP
Jc(j,i):=-1.0/coef;
END LOOP;
END;
END IF;
Jc(N const,i):=-CaractH(i)/Rg/(X(N Const)**2);
END;
END LOOP;
— Create vector F
FOR i IN 1.. N Const LOOP
BEGIN
IF i¡N Const THEN
g:=LOG(X(i))+
CaractH(i)/Rg*(1.0/X(N const)-1.0/CaractT(i));
ELSE
g:=LOG(coef)+
CaractH(i)/Rg*(1.0/X(N const)-1.0/CaractT(i));
END IF;
F(i):=g;
END;
END LOOP;
— Compute DeltaX
686 L. Brunet, J. Caillard & P. Andr´e
SolvJc;
— Newton-Raphson method
FOR i IN 1.. N Const LOOP
X(i):=X(i)-relax*Dx(i);
END LOOP;
— Calcul du residu
residu:=0.0;
FOR i IN 1.. N Const-1 LOOP
residu:=residu+abs(Dx(i));
END LOOP;
EXIT WHEN residu¡Precision;
END LOOP;
— Display results
Put(”——————————–RESULT————————”);
New Line;
coef:=1.0;
FOR j IN 1.. N Const-1 LOOP
coef:=coef-X(j);
Put(”X”); Put(j);Put(”=”);Put(X(j));New Line;
END LOOP;
Put(”X”);Put(N Const);Put(”=”);Put(coef);New Line;
Put(”T(K)=”);Put(X(N Const));Put(” = ”);Put(X(N Const)-To,5,2);Put(” øC”);
END;