#include <iostream>
#include <cmath>
#include "general.hpp"

double** SetPoints(int n,double a,double b)
{
    double** pointsweights = Allocate_matrix(2, 2*n +1); //Initialises 2 arrays, for the points and weights
    double h = (b-a)/(2*n);

    for (int j=0; j<=2*n; j++)
    {
        
        pointsweights[0][j] = a + h*(double)(j); // Creates array of points in the interval
        if (j == 0 || j == 2*n) // creates weights depending on the value of j
        {
            pointsweights[1][j] = h/3.0;
        }
        else if (j % 2 == 1) // Tests if j is odd
        {
            pointsweights[1][j] = 4.0*h/3.0;
        }
        else
        {
            pointsweights[1][j] = 2.0*h/3.0;
        }
    }
    return pointsweights;
}

double CompSimpson(double (*pFunction)(double x), int n, double a, double b)
{
    double** pointsweights = SetPoints(n ,a ,b ); // Uses the previous function
    double sum = 0.0; //Initiates sum to calculate integral.
    for (int i=0; i<=2*n; i++)
    {
        sum +=  (*pFunction)(pointsweights[0][i]) * pointsweights[1][i];
    }
    return sum;
    Deallocate_matrix(2,pointsweights); //Deletes points and weights function.
}
