#pragma once

#include <cmath>

inline int svd_hestenes(int m_m, int n_n, double * a, double * u, double * v, double * sigma)
{
double thr=1.E-4f, nul=1.E-16f;
int n,m,i,j,l,k,lort,iter, in,ll,kk;
double alfa,betta,hamma,eta,t,cos0,sin0,buf,s ;
  n = n_n;
  m = m_m;
	for(i=0;i<n;i++)
		  { in=i*n;
		for(j=0;j<n;j++)
			if(i==j) v[in+j]=1.;
			else v[in+j]=0.;
			}
	for(i=0;i<m;i++)
		  {       in=i*n;
		for(j=0;j<n;j++)
					 {
			u[in+j]=a[in+j];
					  }
			 }

	iter=0;
	while(1)
		{
		 lort=0;
		iter++;
		for(l=0;l<n-1;l++)
			for(k=l+1;k<n;k++)
				{
				alfa=0.; betta=0.; hamma=0.;
				for(i=0;i<m;i++)
					{
						 in=i*n;
						 ll=in+l;
						 kk=in+k;
					alfa += u[ll]*u[ll];
					betta+= u[kk]*u[kk];
					hamma+= u[ll]*u[kk];
					}

				if( sqrt(alfa*betta) < nul )	continue;
				if(fabs(hamma)/sqrt(alfa*betta)<thr) continue;

				lort=1;
				 eta=(betta-alfa)/(2.f*hamma);
				 t=(double)((eta/fabs(eta))/(fabs(eta)+sqrt(1.+eta*eta)));
				 cos0=(double)(1./sqrt(1.+t*t));
				 sin0=t*cos0;

				for(i=0;i<m;i++)
					{
													 in=i*n;
					buf=u[in+l]*cos0-u[in+k]*sin0;
					u[in+k]=u[in+l]*sin0+u[in+k]*cos0;
					u[in+l]=buf;

					if(i>=n) continue;
					buf=v[in+l]*cos0-v[in+k]*sin0;
					v[in+k]=v[in+l]*sin0+v[in+k]*cos0;
					v[in+l]=buf;
					}
				}

		if(!lort) break;
		}

	for(i=0;i<n;i++)
		{
			  s=0.;
		for(j=0;j<m;j++)	s+=u[j*n+i]*u[j*n+i];
		s=(double)sqrt(s);
		sigma[i]=s;
			if( s < nul )	continue;
		for(j=0;j<m;j++)	u[j*n+i]/=s;
		}
//======= Sortirovka ==============
	for(i=0;i<n-1;i++)
		for(j=i;j<n;j++)
			if(sigma[i]<sigma[j])
				{
				  s=sigma[i]; sigma[i]=sigma[j]; sigma[j]=s;
				  for(k=0;k<m;k++)
				  { s= u[i+k*n]; u[i+k*n]=u[j+k*n]; u[j+k*n]=s;}
				  for(k=0;k<n;k++)
				  { s= v[i+k*n]; v[i+k*n]=v[j+k*n]; v[j+k*n]=s;}
				}

	return iter;
}


/**
 * To speedup the original algorithm...
 * 
 * Assumes a is a matrix all elements of which are determined by only
 * a single line: a[i][j] = a[abs(i - j)];
 */
inline int svd_toepliz(int m_m, int n_n, double * a, double * u, double * v, double * sigma)
{
    double thr = 1.E-4f, nul = 1.E-16f; int iter_max = 100000;
    int n, m, i, j, l, k, lort, iter, in, ll, kk;
    double alfa, betta, hamma, eta, t, cos0, sin0, buf, s;
    n = n_n;
    m = m_m;
    for (i = 0; i < n; i++)
    {
        in = i*n;
        for (j = 0; j < n; j++)
            if (i == j) v[in + j] = 1.;
            else v[in + j] = 0.;
    }
    for (i = 0; i < m; i++)
    {
        in = i*n;
        for (j = 0; j < n; j++)
        {
            u[in + j] = a[abs(i - j)];
        }
    }

    iter = 0;
    while (1)
    {
        lort = 0;
        iter++;
        if (iter > iter_max)
        {
            return -1;
        }
        for (l = 0; l < n - 1; l++)
            for (k = l + 1; k < n; k++)
            {
            alfa = 0.; betta = 0.; hamma = 0.;
            for (i = 0; i < m; i++)
            {
                in = i*n;
                ll = in + l;
                kk = in + k;
                alfa += u[ll] * u[ll];
                betta += u[kk] * u[kk];
                hamma += u[ll] * u[kk];
            }

            if (sqrt(alfa*betta) < nul)	continue;
            if (fabs(hamma) / sqrt(alfa*betta) < thr) continue;

            lort = 1;
            eta = (betta - alfa) / (2.f*hamma);
            t = (double) ((eta / fabs(eta)) / (fabs(eta) + sqrt(1. + eta*eta)));
            cos0 = (double) (1. / sqrt(1. + t*t));
            sin0 = t*cos0;

            for (i = 0; i < m; i++)
            {
                in = i*n;
                buf = u[in + l] * cos0 - u[in + k] * sin0;
                u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
                u[in + l] = buf;

                if (i >= n) continue;
                buf = v[in + l] * cos0 - v[in + k] * sin0;
                v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
                v[in + l] = buf;
            }
            }

        if (!lort) break;
    }

    for (i = 0; i < n; i++)
    {
        s = 0.;
        for (j = 0; j < m; j++)	s += u[j*n + i] * u[j*n + i];
        s = (double) sqrt(s);
        sigma[i] = s;
        if (s < nul)	continue;
        for (j = 0; j < m; j++)	u[j*n + i] /= s;
    }
    //======= Sortirovka ==============
    for (i = 0; i < n - 1; i++)
        for (j = i; j < n; j++)
            if (sigma[i] < sigma[j])
            {
        s = sigma[i]; sigma[i] = sigma[j]; sigma[j] = s;
        for (k = 0; k < m; k++)
        {
            s = u[i + k*n]; u[i + k*n] = u[j + k*n]; u[j + k*n] = s;
        }
        for (k = 0; k < n; k++)
        {
            s = v[i + k*n]; v[i + k*n] = v[j + k*n]; v[j + k*n] = s;
        }
            }

    return iter;
}
