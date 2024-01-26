// Gmsh project created on Sun Jan 21 16:21:53 2024

width = 1;
height = 1;
patchoffset = (width/2)*0.2;
patchheightset = (height/2)*0.2;
// Zuerst die Punkte
Point(1) = {0,0,0};
Point(2) = {width,0,0};
Point(3) = {width,height,0};
Point(4) = {0,height,0};
Point(5) = {width/2+patchoffset,0,0};
Point(6) = {width,height/2+patchheightset,0};
Point(7) = {width/2-patchoffset,height,0};
Point(8) = {0,height/2-patchheightset,0};
Point(9) = {width/2,height/2,0};


// Dann die äußeren Linien
Line(1) = {1,5};
Line(2) = {5,2};
Line(3) = {2,6};
Line(4) = {6,3};
Line(5) = {3,7};
Line(6) = {7,4};
Line(7) = {4,8};
Line(8) = {8,1};
// Jetzt die inneren
Line(9) = {5,9};
Line(10) = {6,9};
Line(11) = {7,9};
Line(12) = {8,9};

Curve Loop(100) = {1,9,-12,8};
Curve Loop(101) = {2,3,10,-9};
Curve Loop(102) = {-10,4,5,11};
Curve Loop(103) = {12,-11,6,7};

Surface(100) = {100};
Surface(101) = {101};
Surface(102) = {102};
Surface(103) = {103};

Physical Curve("Left",4)   = {7,8};
Physical Curve("Right",2)  = {3,4};
Physical Curve("Bottom",1) = {1,2};
Physical Curve("Top",3)    = {5,6};

Physical Surface("Domain",1) = {100,101,102,103};

// // Viele Linien, und vorallem Knoten
Transfinite Line {1,9,-12,8}  = 2 Using Progression 1;
Transfinite Line {2,3,10,-9}  = 2 Using Progression 1;
Transfinite Line {-10,4,5,11} = 2 Using Progression 1;
Transfinite Line {12,-11,6,7} = 2 Using Progression 1;


Coherence;
// SetOrder 1;
Recombine Surface {100,101,102,103};
Mesh 2;
RenumberMeshNodes;
Save "patchE.msh4";