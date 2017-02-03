#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
//#include <R_ext/Linpack.
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#define MAT2(mat, i, j, nI) (mat[i+j*nI])
#define MAT3(mat, i, j, a, nI, nJ) (mat[i+j*nI+a*nI*nJ])


void modificationc5(
                    int *p,
                    double *InfoInv,
                    double *nu_stu,
                    double *nu_s_tu,
                    double *mod
                    )
{
    double k1r,k2r,k3r,nu_ab,cons;
    int r,s,t,u;
    
    for (r=0; r<p[0]; r++)
    {
        k2r = 1/MAT2(InfoInv,r,r,p[0]);
        k1r = 0.0;
        k3r = 0.0;
        
        for (s=0; s<p[0]; s++)
        {
            for (t=0; t<p[0]; t++)
            {
                for (u=t; u<p[0]; u++)
                {
                    if(t!=u)cons=2.0;
                    else cons=1.0;
                    
                    k3r = (k3r+cons*MAT2(InfoInv,r,s,p[0])*MAT2(InfoInv,r,t,p[0])*MAT2(InfoInv,r,u,p[0])*
                           MAT3(nu_stu,s,t,u,p[0],p[0])*R_pow_di(k2r,3));
                    nu_ab = (MAT2(InfoInv,t,u,p[0])-MAT2(InfoInv,t,r,p[0])*MAT2(InfoInv,r,u,p[0])
                             /MAT2(InfoInv,r,r,p[0]));
                    k1r = (k1r -0.5*cons*MAT2(InfoInv,r,s,p[0])*k2r*nu_ab*(MAT3(nu_stu,s,t,u,p[0],p[0])
                                                                           +MAT3(nu_s_tu,s,t,u,p[0],p[0])));
                    
                }
            }
        }
        mod[r] = -k1r + k3r/(6*k2r);
    }
}	


void modificationc7(
                    int *p,
                    double *InfoInv,
                    double *nu_stu,
                    double *nu_s_tu,
                    double *mod
                    )
{
    double cons;
    int r,s,t,u,P=p[0];
    
    for (r=0; r<P; r++)
    {
        for (s=0; s<P; s++)
        {
            for (t=0; t<P; t++)
            {
                for (u=t; u<P; u++)
                {
                    if(t==u)cons=1.0;
                    else cons=2.0;
                    mod[r] = (mod[r] +cons*0.5*MAT2(InfoInv,r,s,P)*(MAT2(InfoInv,t,u,P)*(MAT3(nu_stu,s,t,u,P,P)+MAT3(nu_s_tu,s,t,u,P,P))
							-MAT2(InfoInv,r,t,P)*MAT2(InfoInv,r,u,P)*(2.0*MAT3(nu_stu,s,t,u,P,P)/3.0+MAT3(nu_s_tu,s,t,u,P,P))/MAT2(InfoInv,r,r,P)) );
                }
            }
        }

    }
}	
