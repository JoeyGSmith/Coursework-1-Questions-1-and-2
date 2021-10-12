#include <iostream>
#include <cmath>
#include "general.hpp"


double* A_vec_multiplication(double* A0, double* A1, double* A2, int n, double* vec);
double Dot_Product(double* x, double* y, int n);
double vec_Norm(double*x, int n);


// The function takes the matrix A as three vectors, A0 being the lower diagonal,
// A1 is the leading diagonal and A2 is the upper diagonal.

double* BiCGstab(double* A0, double* A1, double* A2, int n, double* b, int maxIter, double Tol)
{
    int count = 0; //This variable wants to suck your blood! Also counts the number
    //of iterations required at the end of the program.

    //Creates an ungodly number of vectors, for the various steps
    //in the BiCGstab method.
    double* rk = Allocate_vector(n);
    double* pPrev = Allocate_vector(n);
    double* pk = Allocate_vector(n);
    double* vk = Allocate_vector(n);
    double* rhat = Allocate_vector(n);
    double* xhat = Allocate_vector(n);
    double* Ax = Allocate_vector(n);
    double* h = Allocate_vector(n);
    double* Ah = Allocate_vector(n);
    double* test = Allocate_vector(n);
    double* s = Allocate_vector(n);
    double* t = Allocate_vector(n);

    double rhoprev = 1.0, alpha = 1.0, omega = 1.0, rhok, beta;

    for (int i=0;i<n;i++) // Sets up zero vector for x0, p0 and v0
    {
        xhat[i] = 0.0;
        pPrev[i] = 0.0;
        vk[i] = 0.0;
    }

    Ax = A_vec_multiplication(A0,A1,A2,n,xhat); // Calculates Ax0

    for (int i=0;i<n;i++) //Creates r0 and rhat
    {
        rk[i] = b[i] - Ax[i];
        rhat[i] = rk[i];
    }

    for (int k=1; k <= maxIter; k++)
    {
        std::cout << "Norm of residual = " << vec_Norm(rk,n) << std::endl;
        rhok = Dot_Product(rhat, rk,n); // step (i)
        beta = (rhok*alpha)/(rhoprev*omega); // step (ii)
        rhoprev = rhok; // stores previous rho for step k+1
        for (int i=0; i<n;i++) // step (iii)
        {
            pk[i] = rk[i] + beta*(pPrev[i] - omega*vk[i]);
            pPrev[i] = pk[i]; //stores previous pk vector for next iteration
        }
        vk = A_vec_multiplication(A0,A1,A2,n,pk); // step (iv)
        alpha = rhok/(Dot_Product(rhat,vk,n)); // step (v)
        for (int i=0; i<n;i++) // step (vi)
        {
            h[i] = xhat[i] + alpha*pk[i];
        }

        Ah = A_vec_multiplication(A0,A1,A2,n,h);
        for (int i=0; i<n;i++) // creates error vector
        {
            test[i] = b[i] - Ah[i];
        }
        if (vec_Norm(test,n) < Tol) // step (vii)
        {
            for (int i=0; i<n;i++)
            {
                xhat[i] = h[i];
            }
            break;
        }
        for (int i=0; i<n;i++) // step (viii)
        {
            s[i] = rk[i] - alpha*vk[i];
        }
        t = A_vec_multiplication(A0,A1,A2,n,s); // step (ix)
        omega = Dot_Product(t,s,n)/Dot_Product(t,t,n); // step (x)
        for (int i=0; i<n;i++) // steps (xi) and (xii)
        {
            xhat[i] = h[i] + omega*s[i];
            rk[i] = s[i] - omega*t[i];
        }

        if (vec_Norm(rk,n) < Tol) // step (xiii)
        {
            break;
        }
        count += 1; // step (xiv)
    }
    std::cout << "Number of iterations required = " << count << std::endl;


    //Deletes all the vectors except xhat, the output.
    Deallocate_vector(rk);
    Deallocate_vector(pPrev);
    Deallocate_vector(pk);
    Deallocate_vector(vk);
    Deallocate_vector(rhat);
    Deallocate_vector(Ax);
    Deallocate_vector(h);
    Deallocate_vector(Ah);
    Deallocate_vector(test);
    Deallocate_vector(s);
    Deallocate_vector(t);
    return xhat;
}


//Uses the three vectors that represent A to compute A times a given vector.
double* A_vec_multiplication(double* A0, double* A1, double* A2, int n, double* vec)
{
    double* AvecProd = Allocate_vector(n);
    AvecProd[0] = A1[0]*vec[0] + A2[0]*vec[1];
    AvecProd[n-1] = A0[n-2]*vec[n-2] + A1[n-1]*vec[n-1];
    for (int i=1; i < n-1; i++)
    {
        AvecProd[i] = A0[i-1]*vec[i-1] + A1[i]*vec[i] + A2[i]*vec[i+1];
    }
    return AvecProd;
}


//Finds the dot product of two n-dimensional vectors
double Dot_Product(double* x, double* y, int n)
{
    double dotproduct=0;

    for (int i=0; i<n; i++)
    {
        dotproduct += x[i]*y[i];
    }
    return dotproduct;
}


//Finds the norm of a given vector
double vec_Norm(double*x, int n)
{
    double norm = 0;

    for (int i=0; i<n;i++)
    {
        norm += pow(x[i],2);
    }
    return pow(norm,0.5);
}
