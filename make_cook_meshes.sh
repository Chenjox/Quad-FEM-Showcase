#!/bin/bash

alias gmsh="C:/Users/Chenjox/Desktop/Daten/Programme/gmsh/gmsh-4.12.1-Windows64/gmsh-4.12.1-Windows64/gmsh.exe"

gmsh - -setnumber num_nodes 6 CooksMembrane.geo

gmsh - -setnumber num_nodes 12 CooksMembrane.geo

gmsh - -setnumber num_nodes 24 CooksMembrane.geo

gmsh - -setnumber num_nodes 48 CooksMembrane.geo

gmsh - -setnumber num_nodes 96 CooksMembrane.geo

gmsh - -setnumber num_nodes 192 CooksMembrane.geo

gmsh - -setnumber num_nodes 384 CooksMembrane.geo