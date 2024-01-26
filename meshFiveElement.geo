// Gmsh project created on Sun Jan 21 16:21:53 2024

width = 1;
height = 1;
patchoffset = (width/2)*0.6;
patchheightset = (height/2)*0.5;

// Zuerst die Punkte
Point(1) = {0,0,0};
Point(2) = {width,0,0};
Point(3) = {width,height*1.5,0};
Point(4) = {0,height,0};
Point(5) = {patchoffset,patchheightset,0};
Point(6) = {width-patchheightset,patchoffset,0};
Point(7) = {width-patchoffset,height-patchheightset,0};
Point(8) = {patchheightset,height-patchoffset,0};


// Dann die äußeren Linien
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// diagonalen
Line(5) = {1,5};
Line(6) = {2,6};
Line(7) = {3,7};
Line(8) = {4,8};
// innere
Line(9) = {5,6};
Line(10) = {6,7};
Line(11) = {7,8};
Line(12) = {8,5};

// unteres Element
Curve Loop(100) = {1,2,3,4};

// Element der Mitte
Curve Loop(105) = {9,10,11,12};

Surface(100) = {100};

Physical Curve("Left",4)   = {4};
Physical Curve("Right",2)  = {2};
Physical Curve("Bottom",1) = {1};
Physical Curve("Top",3)    = {3};

Physical Surface("Domain",1) = {100};

// // Viele Linien, und vorallem Knoten
Transfinite Line {1,2,3,4}  = 2 Using Progression 1;


Coherence;
// SetOrder 1;
Recombine Surface {100};
Mesh 2;
RenumberMeshNodes;
Save "patchC.msh4";