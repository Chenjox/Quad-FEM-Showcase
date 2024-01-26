// Gmsh project created on Sun Jan 21 16:21:53 2024

width = 1;
height = 1;
patchoffset = (width/2)*0.8;
patchheightset = (height/2);

// Zuerst die Punkte
Point(1) = {0,0,0};
Point(2) = {width,0,0};
Point(3) = {width,height+patchheightset,0};
Point(4) = {0,height,0};


// Dann die äußeren Linien
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(100) = {1,2,3,4};

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
Save "patchD.msh4";