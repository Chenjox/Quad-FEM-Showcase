# Quad-FEM-Showcase
A small rust program to show facettes of FEM with quadrilaterial elements in linear isotropic elasticity



# Quad Element Indexing

Numbers on Verteces denote Node Numbers
Numbers on top of Edges denote Edge Numbers

Edge evaluations will be conforming to the reference directions

```
       2
 3 ->- - ->- 2
 
 |     2     |
 ^     ^     ^
 |     |     |
       |
3|     o-->1 |1
 
 |           |
 ^           ^
 |           |
 
 0 ->- - ->- 1
       0
```