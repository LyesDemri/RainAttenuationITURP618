function hR = ituRP839(lon,lat)

load h0.mat;

%on passe de la latitude et la longitude à la ligne et colonne de la
%matrice
colonne = round(lon/1.5)+1;
ligne = round(-lat/1.5)+61;

hR = h0(ligne,colonne) + 0.36;