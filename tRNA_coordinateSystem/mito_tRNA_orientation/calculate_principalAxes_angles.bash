#!/bin/bash


vector1=$1
vector2=$2

if [ "$vector1" = "AtRNA_acceptor" ] && [ "$vector2" = "PtRNA_acceptor" ]
then
	referencevector=AtRNA_ASL
elif [ "$vector1" = "AtRNA_ASL" ] && [ "$vector2" = "PtRNA_ASL" ]
then
	referencevector=AtRNA_acceptor
fi

#paxis3 contains the minor inertia axis.

#calculate cross product between the 2 vectors: 
paste $vector1\_paxis3.xvg $vector2\_paxis3.xvg | awk 'function acos(x) {return atan2(sqrt(1-x*x), x)}{magA=sqrt($2^2+$3^2+$4^2); magB=sqrt($6^2+$7^2+$8^2); Cx=($3*$8-$4*$7); Cy=($4*$6-$2*$8); Cz=($2*$7-$3*$6); magC=sqrt(Cx^2+Cy^2+Cz^2); AdotB=($2*$6+$3*$7+$4*$8); angle=180/3.141592*acos(AdotB/(magA*magB)); print Cx,Cy,Cz,angle}' > $vector1\_X_$vector2    #format of this file: Cx,Cy,Cz,angle in degrees

#calculate  the overlap of the cross product with the reference vector to determine the angle
paste $vector1\_X_$vector2 $referencevector\_paxis3.xvg | awk 'function acos(x) { return atan2(sqrt(1-x*x), x) }{CdotREF=($1*$6+$2*$7+$3*$8); magC=sqrt($1^2+$2^2+$3^2); magREF=($6^2+$7^2+$8^2); angle2=180/3.141592*acos(CdotREF/(magC*magREF)); if(CdotREF>0){dir=+1}; if(CdotREF<0){dir=-1}; print $5,angle2,dir*$4}' > $vector1\_$vector2\_angles   #format: time, angle between cross product vector and reference vector(should be near 0 or near 180), angle between $vector1 and $vector2 ... both in degrees


