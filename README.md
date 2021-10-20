# TrafficScienceProject
A set of FORTRAN tools to evaluate traffic behavior and its conection to Evolutionary Game Theory

I do not garantee the behavior of any routine in the program, but if you find a bug or want to share something you can reach me in thisisajobrelatedmail@gmail.com

I pray you have as much fun running it as I had making it and feel free to make any modification you please, but I ask you to make it publically available as I did. 

I used the GNU fortran compiler (gfortran). After you uncomment the subroutine you wish to run in the MAIN section on the final lines, I recomend you go to the folder where the program is saved and run: 

$ gfortran  -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace Teste.f95 -o TrafficProject

$ ./TrafficProject

All data is compatible with the graphic interface I use, the gnuplot. You can change in the program if this cause problems for you


Best regards and have fun!


