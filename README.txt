This is a program that searches for reduced spherical pictures over certain
1-generator, 1-relator relative group presentations.  Usage instructions can be
found by running the program with no input.

Here is the input required to find pictures for three interesting presentations.

<G,x | x^2*g*x^-1*g^3> where G = <g | g^9>:
  -o 9 -r xxX -e 0 1 3 -w BBB -n 6 -h 27

<G,x | x^2*g*x^-1*g^-3> where G = <g | g^9>:
  -o 9 -r xxX -e 0 1 -3 -w B -n 7 -h 15

<G,x | x^3*g*x^-1*g^-3> where G = <g | g^6>:
  -o 6 -r xxxX -e 0 0 2 1 -w DD -n 3 -h 10

-Matthias Merzenich
