/*****************************************************************************
 * intrapot     Function for implementation of intra-molecular potentials    *
 *                                                                           *
 * BendForce()  Nimmt Position dreier Teilchen entgegen und berechnet        *
 *              aufgrund des Bending-Potentials des eingeschlossenen Winkels *
 *              die jeweilige Kraft auf die Teilchen.                        *
 * TorForce()   Nimmt die Positionen von vier aufeinanderfolgenden Teilchen  *
 *              entgegen (1-2-3-4) und berechnet Kraft und Potential auf-    *
 *              grund der Drehung um die 2-3 Achse.                          *
 * Chedin2Force(),							     *
 * Chedin34Force(),							     *
 * Chedin5Force() and							     *   
 * Chedin6Force() take the coordinates of 3 particles (C-O-O) and, thinking  *
 *  		  about them to be carbon dioxide, calculate forces and      *
 *		  potential in correspondence with Chedin formula on         *
 *		  2th, 3-4th, 5th and 6th order terms respectively.	     *
 * BenzeneForce() Berechnet die Kraefte aufgrund von internen Koordinaten    *
 *                im Benzolring.                                             *
 *****************************************************************************
 * $Log:$
 *
 */

#ifndef lint
static char *RCSid = "$Header:$";
#endif

/*****************************************************************************
 * Program include files                                                     *
 *****************************************************************************/
#include        "defs.h"
#include        "structs.h"
#include        "messages.h"

/*****************************************************************************
 * Library include files                                                     *
 *****************************************************************************/
#include        <math.h>

/*****************************************************************************
 * External function declarations                                            *
 *****************************************************************************/
#if defined(ANSI) || defined(__STDC__)
 void		message(int *, ...);	/* Write a warning or error message  */
#else
 void		message();		/* Write a warning or error message  */
#endif
void	KreuzProd();

/*****************************************************************************
 * Macros                                                                    *
 *****************************************************************************/
#define C(a,b) SCALPROD(a,b)
#define D(a,b) (C(a,a)*C(b,b)-SQR(C(a,b)))

static void VecProd(a, b, c)    /* ACHTUNG: Fall a/b = c */
 vec_mp         a, b, c;
{
  (*c)[0] = ( (*a)[1] * (*b)[2] - (*a)[2] * (*b)[1] );
  (*c)[1] = ( (*a)[2] * (*b)[0] - (*a)[0] * (*b)[2] );
  (*c)[2] = ( (*a)[0] * (*b)[1] - (*a)[1] * (*b)[0] );
}

static void da(r1, r0, d1)
 vec_mp         r1, r0, d1;
{
  (*d1)[0] = (*r1)[0] - (*r0)[0];
  (*d1)[1] = (*r1)[1] - (*r0)[1];
  (*d1)[2] = (*r1)[2] - (*r0)[2];
}

/*****************************************************************************
 * AngleForce                                                                *
 *  Potential v(phi) = 1/2 k_phi ( phi - phi_0 )^2                           *
 *  phi = PI - theta mit theta nach Allen/Tildesley p329                     *
 *****************************************************************************/
void AngleForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte (out)             */
 double         *pot,                           /* Energie (out)             */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = k_phi            */
                                                /* par[1] = phi_0            */
{
  vec_mt         d0, d1;
  double         theta,                /* Winkel nach Allen/Tildesley s329   */
                 c1, c2, c3,           /* temporaere Variable                */
                 dv;                   /* d v(phi)/d cos(theta)              */
  double         k_phi = par[0],
                 phi_0 = par[1];
  
  da(r+2, r+1, (vec_mp)d1);
  da(r+1, r, (vec_mp)d0);  

  theta = acos(SCALPROD(d1,d0)/(ABS3(d1)*ABS3(d0))); /* Nenner duerfte un-   */
  						     /* problematisch sein,  */
  						     /* d's bei vernuenftigen*/
  						     /* Molekuelen endlich   */
  *pot += 0.5 * k_phi * SQR(PI-theta-phi_0);

  dv = k_phi/sin(theta)*(PI-theta-phi_0); /* 1.0/sin(x) produziert bei HP's  */
  					  /* ab x < 1.0e-75 UNDERFLOW !      */
  c1 = -dv * ( -1.0/sqrt((C(d1,d1)*C(d0,d0))));
  c2 = C(d1,d0)/C(d1,d1);
  
  f[2][0] = c1 * ( c2 * d1[0] - d0[0] );
  f[2][1] = c1 * ( c2 * d1[1] - d0[1] );
  f[2][2] = c1 * ( c2 * d1[2] - d0[2] );
  
  c3 = C(d1,d0)/C(d0,d0);
  
  f[1][0] = -c1 * ( (c2 + 1) * d1[0] - (c3 + 1) * d0[0] );
  f[1][1] = -c1 * ( (c2 + 1) * d1[1] - (c3 + 1) * d0[1] );
  f[1][2] = -c1 * ( (c2 + 1) * d1[2] - (c3 + 1) * d0[2] );
  
  f[0][0] = -c1 * ( c3 * d0[0] - d1[0] );
  f[0][1] = -c1 * ( c3 * d0[1] - d1[1] );
  f[0][2] = -c1 * ( c3 * d0[2] - d1[2] );
}

/*****************************************************************************
 * TorForce                                                                  *
 *                                                                           *
 * Torsionspotential unter Verwendung von 4 Sites.                           *
 *                                                                           *
 *  V_tor = 1/2 V_2 (1-cos 2(phi-phi0))                                      *
 *  Berechnung der Kraefte nach Allen/Tildesley p330                         *
 *****************************************************************************/
void TorForce(r, f, par, pot)
 vec_mp  r,                  /* die 4 Ortskoordinaten                  (in)  */
         f;                  /* resultierende Kraefte auf die Teilchen (out) */
 double  *pot,               /* Potential                              (out) */
         par[];              /* Potentialparameter                           */
                             /* par[0] = V_2                                 */
                             /* par[1] = phi0                                */
{
  vec_mt         d0, d1, d2,
                 t, u, tu;
  double         phi, c1, c2, c3, c4,   /* Zur Aufbewahrung von Zwischen-    */
                 c5, c6, c00, c01, c11, /* ergebnissen.                      */
                 c02, c12, c22, c5d1,
                 c5d0, d01, d12,
                 dv;                    /* d v(phi)/d cos(phi)               */
  /*
  vec_mt static	 old_f[4];
   */
  
  
  da(r+1, r, (vec_mp)d0);
  da(r+2, r+1, (vec_mp)d1);
  da(r+3, r+2, (vec_mp)d2);
  
  c00 = C(d0,d0);
  c01 = C(d0,d1);
  c02 = C(d0,d2);
  c11 = C(d1,d1);
  c12 = C(d1,d2);
  c22 = C(d2,d2);
  d01 = D(d1,d0);
  d12 = D(d1,d2);
  
  c6   =  -1.0 / sqrt(d12 * d01);               /* kein underflow bei sqrt ! */
  c5   =  ( c12 * c01 - c02 * c11);
  
  c4 = c5 * c6;
  
  if ((c4 > -1.0) && (c4 < 1.0))              /* manchmal ist c4 = 1.0 + eps */
   phi  =  acos(c4);
  else
   if (c4 >= 1.0)
      phi = 0.0;
   else
      phi = PI;
  
  KreuzProd(d0, d1, t);
  KreuzProd(d1, d2, u);
  KreuzProd(t, u, tu);
  
  if ( SCALPROD(d1, tu) > 0.0 )
   phi = -phi;   
  
  *pot +=  0.5*par[0]*(1.0-cos(2.0*(phi-par[1])));
  
  /*  if ( (phi >= 0.005) && (phi <= PI - 0.005) ) */
  
  
  if ( (phi != 0.0) && (phi != PI) )  
  {
    dv  = -2.0*par[0]*cos(phi-par[1])*(cos(par[1])-sin(par[1])/tan(phi));       
                                               /* d v(phi)/d cos(phi)        */

    /* Behandlung der Singularitaet fuer phi = n * pi: */
    /* Falls tan(phi) == 0 wird, verwende Kraefte aus  */
    /* dem letzten Zeitschritt. Resultat:              */
    /* Ausser fuer T == 0 werden die 4 Sites niemals   */
    /* in einer Ebene zu liegen kommen, sondern immer  */
    /* schwingen.                                      */
    /* Dies werden sie jedoch bei endlicher Tempe-     */
    /* ratur ohnehin; bleibt die Frage nach der        */
    /* Unstetigkeit der Kraefte.                       */
    /*                                                 */
    /* GEAENDERT, 10/95: Obiges Vorgehen fuehrt zu un- */
    /* physikalischem Verhalten, insb. bei kleinen     */
    /* Temperaturen. Jetzt werden die Kraefte einfach  */
    /* auf 0 gesetzt.                                  */
                                              
    c6 *= -dv;
  
    c1 =  -c12 * c6;
    c2 =   c11 * c6;
  
    c5d0 = c5/d01;
    c5d1 = c5/d12;
    
    c3 =  c5d0 * c11 * c6;
    c4 = -c5d0 * c01 * c6;
    
    /*
    old_f[0][0] = f[0][0] = ( c1 + c4 ) * d1[0] + c2 * d2[0] + c3 * d0[0];
    old_f[0][1] = f[0][1] = ( c1 + c4 ) * d1[1] + c2 * d2[1] + c3 * d0[1];
    old_f[0][2] = f[0][2] = ( c1 + c4 ) * d1[2] + c2 * d2[2] + c3 * d0[2];
     */
    
    
    f[0][0] = ( c1 + c4 ) * d1[0] + c2 * d2[0] + c3 * d0[0];
    f[0][1] = ( c1 + c4 ) * d1[1] + c2 * d2[1] + c3 * d0[1];
    f[0][2] = ( c1 + c4 ) * d1[2] + c2 * d2[2] + c3 * d0[2];
    
    c1 = c6 * ( -c01 - c11 - c5d1 * c12 );
    c2 = c6 * ( c12 + 2 * c02 + ( c00 + c01) * c5d0 + c5d1 * c22 );
    c3 = c6 * ( -c12 - c5d0 * ( c11 + c01));
    
    /*
    old_f[1][0] = f[1][0] = c1 * d2[0] + c2 * d1[0] + c3 * d0[0];
    old_f[1][1] = f[1][1] = c1 * d2[1] + c2 * d1[1] + c3 * d0[1];
    old_f[1][2] = f[1][2] = c1 * d2[2] + c2 * d1[2] + c3 * d0[2];
     */
    
    
    f[1][0] = c1 * d2[0] + c2 * d1[0] + c3 * d0[0];
    f[1][1] = c1 * d2[1] + c2 * d1[1] + c3 * d0[1];
    f[1][2] = c1 * d2[2] + c2 * d1[2] + c3 * d0[2];
    
    c1 = c6 * ( c01 + c5d1 * ( c11 + c12 ));
    c2 = c6 * ( -c01 - 2 * c02 - c5d0 * c00 - c5d1 * (c22 + c12) );
    c3 = c6 * ( c12 + c11 + c5d0 * c01 );
    
    /*
    old_f[2][0] = f[2][0] = c1 * d2[0] + c2 * d1[0] + c3 * d0[0];
    old_f[2][1] = f[2][1] = c1 * d2[1] + c2 * d1[1] + c3 * d0[1];
    old_f[2][2] = f[2][2] = c1 * d2[2] + c2 * d1[2] + c3 * d0[2];
     */
    
    
    f[2][0] = c1 * d2[0] + c2 * d1[0] + c3 * d0[0];
    f[2][1] = c1 * d2[1] + c2 * d1[1] + c3 * d0[1];
    f[2][2] = c1 * d2[2] + c2 * d1[2] + c3 * d0[2];
    
    c1 =  c6 * ( c01 + c5d1 * c12 );
    c2 = -c6 * c11;
    c3 = -c6 * c5d1 * c11;
    
    /*
    old_f[3][0] = f[3][0] = c1 * d1[0] + c2 * d0[0] + c3 * d2[0];
    old_f[3][1] = f[3][1] = c1 * d1[1] + c2 * d0[1] + c3 * d2[1];
    old_f[3][2] = f[3][2] = c1 * d1[2] + c2 * d0[2] + c3 * d2[2];
     */
    
    
    f[3][0] = c1 * d1[0] + c2 * d0[0] + c3 * d2[0];
    f[3][1] = c1 * d1[1] + c2 * d0[1] + c3 * d2[1];
    f[3][2] = c1 * d1[2] + c2 * d0[2] + c3 * d2[2];
  }
  else
  {
     /*
     f[0][0] = old_f[0][0]; f[0][1] = old_f[0][1]; f[0][2] = old_f[0][2];
     f[1][0] = old_f[1][0]; f[1][1] = old_f[1][1]; f[1][2] = old_f[1][2];
     f[2][0] = old_f[2][0]; f[2][1] = old_f[2][1]; f[2][2] = old_f[2][2];
     f[3][0] = old_f[3][0]; f[3][1] = old_f[3][1]; f[3][2] = old_f[3][2];     
      */
     
     
     f[0][0] = f[1][0] = f[2][0] = f[3][0] = 0.0;
     f[0][1] = f[1][1] = f[2][1] = f[3][1] = 0.0;
     f[0][2] = f[1][2] = f[2][2] = f[3][2] = 0.0;
/*     f[0][3] = f[1][3] = f[2][3] = f[3][3] = 0.0;*/
  }

}

/*****************************************************************************
 * BellForce                                                                 *
 *                                                                           *
 * Torsionspotential nach Bell (1945). Unter Einbeziehung von 6 Sites.       *
 *                                                                           *
 * Reihenfolge der Sites: i, j, k, l, m, n                                   *
 *                                                                           *
 *         m             n                                                   *
 *          \           /                                                    *
 *           j---------k                                                     *
 *          /           \                                                    *
 *         i             l                                                   *
 *                                                                           *
 *                                                                           *
 *  V = 1/2 V_2 (1-cos (Phi_H + Phi_C))                                      *
 *                                                                           *
 *  Berechnung der Kraefte nach Allen/Tildesley p330                         *
 *****************************************************************************/
void BellForce(r, f, par, pot)
   vec_mp  r,                  /* die 6 Ortskoordinaten                  (in)  */
           f;                  /* resultierende Kraefte auf die Teilchen (out) */
   double  *pot,               /* Potential                              (out) */
           par[];              /* Potentialparameter                           */
                               /* par[0] = V_2                                 */
{
   vec_mt         d0, d1, d2, t, u, tu;
   double         phi[2], c1, c2, c3, c4,  /* Zur Aufbewahrung von Zwischen-  */
      c5, c6, c00, c01, c11,               /* ergebnissen.                    */
      c02, c12, c22, c5d1,
      c5d0, d01, d12,
      dv;                                /* d v(phi)/d cos(phi)               */
   int  i,j;
   
   
      
   f[0][0] = f[1][0] = f[2][0] = f[3][0] = f[4][0] = f[5][0] = 0.0;
   f[0][1] = f[1][1] = f[2][1] = f[3][1] = f[4][1] = f[5][1] = 0.0;
   f[0][2] = f[1][2] = f[2][2] = f[3][2] = f[4][2] = f[5][2] = 0.0;
   
   for ( i = 0; i < 2; i++)
   {
      
      if ( i == 0 ) {
       da(r+1, r,   (vec_mp)d0);
       da(r+2, r+1, (vec_mp)d1);
       da(r+3, r+2, (vec_mp)d2);
      }
      else {
       da(r+1, r+4, (vec_mp)d0);
       da(r+2, r+1, (vec_mp)d1);
       da(r+5, r+2, (vec_mp)d2);	 	 
      }
      
      c00 = C(d0,d0);
      c01 = C(d0,d1);
      c02 = C(d0,d2);
      c11 = C(d1,d1);
      c12 = C(d1,d2);
      c22 = C(d2,d2);
      d01 = D(d1,d0);
      d12 = D(d1,d2);
      
      c6   =  -1.0 / sqrt(d12 * d01);            /* kein underflow bei sqrt !   */
      c5   =  ( c12 * c01 - c02 * c11);
      
      c4 = c5 * c6;
      
      if ((c4 > -1.0) && (c4 < 1.0))             /* manchmal ist c4 = 1.0 + eps */
	 phi[i]  =  acos(c4);
      else
	 if (c4 >= 1.0)
	    phi[i] = 0.0;
	 else
	    phi[i] = PI;
      
      KreuzProd(d0, d1, t);
      KreuzProd(d1, d2, u);
      KreuzProd(t, u, tu);
  
      if ( SCALPROD(d1, tu) > 0.0 )
       phi[i]= -phi[i];         
      
      if ( (phi[i] != 0.0) && (phi[i] != PI) )
      {	 
	 c6 /= sin(phi[i]);
	 
	 c1 =  -c12 * c6;
	 c2 =   c11 * c6;
	 
	 c5d0 = c5/d01;
	 c5d1 = c5/d12;
	 
	 c3 =  c5d0 * c11 * c6;
	 c4 = -c5d0 * c01 * c6;
	 
	 f[0 + 4 * i][0] += ( c1 + c4 ) * d1[0] + c2 * d2[0] + c3 * d0[0];
	 f[0 + 4 * i][1] += ( c1 + c4 ) * d1[1] + c2 * d2[1] + c3 * d0[1];
	 f[0 + 4 * i][2] += ( c1 + c4 ) * d1[2] + c2 * d2[2] + c3 * d0[2];
	 
	 c1 = c6 * ( -c01 - c11 - c5d1 * c12 );
	 c2 = c6 * ( c12 + 2 * c02 + ( c00 + c01) * c5d0 + c5d1 * c22 );
	 c3 = c6 * ( -c12 - c5d0 * ( c11 + c01));
	 
	 f[1][0] += c1 * d2[0] + c2 * d1[0] + c3 * d0[0];
	 f[1][1] += c1 * d2[1] + c2 * d1[1] + c3 * d0[1];
	 f[1][2] += c1 * d2[2] + c2 * d1[2] + c3 * d0[2];
	 
	 c1 = c6 * ( c01 + c5d1 * ( c11 + c12 ));
	 c2 = c6 * ( -c01 - 2 * c02 - c5d0 * c00 - c5d1 * (c22 + c12) );
	 c3 = c6 * ( c12 + c11 + c5d0 * c01 );
	 
	 f[2][0] += c1 * d2[0] + c2 * d1[0] + c3 * d0[0];
	 f[2][1] += c1 * d2[1] + c2 * d1[1] + c3 * d0[1];
	 f[2][2] += c1 * d2[2] + c2 * d1[2] + c3 * d0[2];
	 
	 c1 =  c6 * ( c01 + c5d1 * c12 );
	 c2 = -c6 * c11;
	 c3 = -c6 * c5d1 * c11;
	 
	 f[3 + 2 * i][0] += c1 * d1[0] + c2 * d0[0] + c3 * d2[0];
	 f[3 + 2 * i][1] += c1 * d1[1] + c2 * d0[1] + c3 * d2[1];
	 f[3 + 2 * i][2] += c1 * d1[2] + c2 * d0[2] + c3 * d2[2];
      }
   }
   
   *pot +=  0.5 * par[0] * ( 1.0 - cos( phi[0] + phi[1] ) );
  
   dv  = 0.5 * par[0] * sin(phi[0] + phi[1]);       
   
   for ( i = 0; i < 6; i++)
    for ( j = 0; j < 3; j++)
     f[i][j] *= dv;
}

/*****************************************************************************
 * WagForce                                                                  *
 *  Potential nach B. Mannfors 1991                                          *
 *            v(r) = 1/2 K ( theta - theta_0 )^2                             *
 *                                                                           *
 *              2                                                            *
 *     ^       /                                                             *
 *  theta 1 - 0  } phi                                                       *
 *     |       \                                                             *
 *              3                                                            *
 *                                                                           *
 *****************************************************************************/
void WagForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = K         (in)   */
                                                /* par[1] = theta_0   (in)   */
{
   int         	j;
   double 	theta, phi, tan_theta, sin_phi,
        	cos_phi, cos_theta, tan_theta_sin_sqr,
		cos_theta_sin_phi,
		r_ab, r_ac, r_ad;
   vec_mt       d0, d1, d2, k;

   da(r+1, r, (vec_mp)d0);
   da(r+2, r, (vec_mp)d1);
   da(r+3, r, (vec_mp)d2);

   r_ab = ABS3(d0);
   r_ac = ABS3(d1);
   r_ad = ABS3(d2);

   for ( j = 0; j < 3; j++)
   {
    d0[j] /=  r_ab;
    d1[j] /=  r_ac;
    d2[j] /=  r_ad;
   }

   cos_phi = SCALPROD(d1, d2);
   phi     = acos( cos_phi );
   sin_phi = sin(phi);

   KreuzProd( d2, d1, k);

   theta     = SCALPROD(d0, k) / sin_phi;
   theta     = asin(theta);
   tan_theta = tan(theta);
   cos_theta = cos(theta);
   cos_theta_sin_phi = cos_theta * sin_phi;
   tan_theta_sin_sqr = tan_theta / SQR(sin_phi);

   r_ab = - par[0] * (theta - par[1]) / r_ab;
   r_ac = - par[0] * (theta - par[1]) / r_ac;
   r_ad = - par[0] * (theta - par[1]) / r_ad;

   /*
   r_ab = r_ac = r_ad = 1.0;
    */

   f[1][0] = ((d2[1] * d1[2] - d2[2] * d1[1])/cos_theta_sin_phi
	      - d0[0] * tan_theta) * r_ab;

   f[1][1] = ((d2[2] * d1[0] - d2[0] * d1[2])/cos_theta_sin_phi
	      - d0[1] * tan_theta) * r_ab;

   f[1][2] = ((d2[0] * d1[1] - d2[1] * d1[0])/cos_theta_sin_phi
	      - d0[2] * tan_theta) * r_ab;

   f[2][0] = ((d0[1] * d2[2] - d0[2] * d2[1])/cos_theta_sin_phi
	      - tan_theta_sin_sqr * ( d1[0] - d2[0] * cos_phi )) * r_ac;

   f[2][1] = ((d0[2] * d2[0] - d0[0] * d2[2])/cos_theta_sin_phi
	      - tan_theta_sin_sqr * ( d1[1] - d2[1] * cos_phi )) * r_ac;

   f[2][2] = ((d0[0] * d2[1] - d0[1] * d2[0])/cos_theta_sin_phi
	      - tan_theta_sin_sqr * ( d1[2] - d2[2] * cos_phi )) * r_ac;

   f[3][0] = ((d1[1] * d0[2] - d1[2] * d0[1])/cos_theta_sin_phi
	      - tan_theta_sin_sqr * ( d2[0] - d1[0] * cos_phi )) * r_ad;

   f[3][1] = ((d1[2] * d0[0] - d1[0] * d0[2])/cos_theta_sin_phi
	      - tan_theta_sin_sqr * ( d2[1] - d1[1] * cos_phi )) * r_ad;

   f[3][2] = ((d1[0] * d0[1] - d1[1] * d0[0])/cos_theta_sin_phi
	      - tan_theta_sin_sqr * ( d2[2] - d1[2] * cos_phi )) * r_ad;

   f[0][0] =  (- f[1][0] - f[2][0] - f[3][0]);
   f[0][1] =  (- f[1][1] - f[2][1] - f[3][1]);
   f[0][2] =  (- f[1][2] - f[2][2] - f[3][2]);

   (*pot) += 0.5 * par[0] * SQR(theta - par[1]);
}

/*****************************************************************************
 * StretchForce                                                              *
 *  Potential v(r) = 1/2 K ( |r_1 - r_2| - R_0 )^2                           *
 *****************************************************************************/
void StretchForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = K         (in)   */
                                                /* par[1] = phi_0     (in)   */
{
  vec_mt		r_diff;
  double		r_abs, fak1;

  r_diff[0] = r[1][0] - r[0][0];
  r_diff[1] = r[1][1] - r[0][1];
  r_diff[2] = r[1][2] - r[0][2];

  r_abs = ABS3(r_diff);

  fak1 = par[0] * (r_abs - par[1])/r_abs;
  
  f[0][0] = r_diff[0] * fak1;
  f[0][1] = r_diff[1] * fak1;
  f[0][2] = r_diff[2] * fak1;
  f[1][0] = -r_diff[0] * fak1;
  f[1][1] = -r_diff[1] * fak1;
  f[1][2] = -r_diff[2] * fak1;
  
  (*pot) += 0.5 * par[0] * SQR(r_abs - par[1]);
 
}

/*****************************************************************************
 * Stretch-Stretch Force                                                              *
 *  Potential v(r) = K ( |r_1 - r_2| - R_012 )*( |r_1 - r_3| - R_013 )                          *
 *****************************************************************************/
void StretchCrossForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = K         (in)   */
                                                /* par[1] = R_012     (in)   */
                                                /* par[2] = R_013     (in)   */
{
  vec_mt		r_diff1, r_diff2;
  double		r_abs1, r_abs2, fak1, fak2;

  r_diff1[0] = r[1][0] - r[0][0];
  r_diff1[1] = r[1][1] - r[0][1];
  r_diff1[2] = r[1][2] - r[0][2];
  r_diff2[0] = r[2][0] - r[0][0];
  r_diff2[1] = r[2][1] - r[0][1];
  r_diff2[2] = r[2][2] - r[0][2];

  r_abs1 = ABS3(r_diff1);
  r_abs2 = ABS3(r_diff2);

  fak1 = - par[0] * (r_abs2 - par[2])/r_abs1;
  fak2 = - par[0] * (r_abs1 - par[1])/r_abs2;
  
  f[1][0] = r_diff1[0] * fak1;
  f[1][1] = r_diff1[1] * fak1;
  f[1][2] = r_diff1[2] * fak1;
  f[2][0] = r_diff2[0] * fak2;
  f[2][1] = r_diff2[1] * fak2;
  f[2][2] = r_diff2[2] * fak2;
  f[0][0] = - f[1][0] - f[2][0];
  f[0][1] = - f[1][1] - f[2][1];
  f[0][2] = - f[1][2] - f[2][2];
  
  (*pot) += par[0] * (r_abs1 - par[1]) * (r_abs2 - par[2]);
}

/*****************************************************************************
 * XY2Z2SymStrForce                                                          *
 *  Potential v(r) = 1/2 Sum (K_ij S_i S_j), i,j = 2, 4, 6, 9                *
 *    (Wang et al. J.Phys.Chem. A101(1997)5919)				     *
 *****************************************************************************/
void XY2Z2SymStrForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0..4] = K_ij (in)     */
                                                /* par[5..6] = r_eq (in)     */
{
  vec_mt	n_12, n_13, n_14, n_15;
  double	r_abs12, r_abs13, r_abs14, r_abs15, 
		sym[4], fak[4], temp;

					/* Definition of the local coords   */
  da(r+1, r, (vec_mp)n_12);  		/* r_X(1)Y(2) 			    */
  da(r+2, r, (vec_mp)n_13);  		/* r_X(1)Y(3)			    */
  da(r+3, r, (vec_mp)n_14); 		/* r_X(1)Z(4)			    */
  da(r+4, r, (vec_mp)n_15);  		/* r_X(1)Z(5)			    */

  r_abs12 = ABS3(n_12);
  r_abs13 = ABS3(n_13);
  r_abs14 = ABS3(n_14);
  r_abs15 = ABS3(n_15);
					 /* Symmetry coordinates 	    */
  sym[0] = r_abs14 + r_abs15 - 2*par[6]; /* SQRT(2) is transported from     */
  sym[1] = r_abs12 + r_abs13 - 2*par[5]; /* denominators to pot expression  */
  sym[2] = r_abs14 - r_abs15;
  sym[3] = r_abs12 - r_abs13;

  r_abs12 = 1.0/r_abs12;
  r_abs13 = 1.0/r_abs13;
  r_abs14 = 1.0/r_abs14;
  r_abs15 = 1.0/r_abs15;

  n_12[0] = n_12[0]*r_abs12; 		/* Unit-vector n_12 		    */
  n_12[1] = n_12[1]*r_abs12; 
  n_12[2] = n_12[2]*r_abs12; 
  n_13[0] = n_13[0]*r_abs13; 		/* Unit-vector n_13		    */
  n_13[1] = n_13[1]*r_abs13; 
  n_13[2] = n_13[2]*r_abs13; 
  n_14[0] = n_14[0]*r_abs14; 		/* Unit-vector n_14 		    */
  n_14[1] = n_14[1]*r_abs14; 
  n_14[2] = n_14[2]*r_abs14; 
  n_15[0] = n_15[0]*r_abs15; 		/* Unit-vector n_15		    */
  n_15[1] = n_15[1]*r_abs15; 
  n_15[2] = n_15[2]*r_abs15; 

 (*pot) += 0.25 * ( par[0]*SQR(sym[0]) + par[1]*SQR(sym[1])  
		 + par[2]*SQR(sym[2]) + par[3]*SQR(sym[3]) )
		 + 0.5*par[4]*sym[0]*sym[1];

 fak[0] = par[0]*sym[0] + par[4]*sym[1]; 
 fak[1] = par[4]*sym[0] + par[1]*sym[1]; 
 fak[2] = par[2]*sym[2]; 
 fak[3] = par[3]*sym[3]; 

 temp = (fak[1] + fak[3]) * 0.5; /* 2= SQRT(2) from sym. coordinates defs   */
 f[1][0] = - temp*n_12[0];       /* multiplied by SQRT(2) from f defs       */
 f[1][1] = - temp*n_12[1];
 f[1][2] = - temp*n_12[2];

 temp = (fak[1] - fak[3]) * 0.5;
 f[2][0] = - temp*n_13[0];
 f[2][1] = - temp*n_13[1];
 f[2][2] = - temp*n_13[2];

 temp = (fak[0] + fak[2]) * 0.5;
 f[3][0] = - temp*n_14[0];
 f[3][1] = - temp*n_14[1];
 f[3][2] = - temp*n_14[2];

 temp = (fak[0] - fak[2]) * 0.5;
 f[4][0] = - temp*n_15[0];
 f[4][1] = - temp*n_15[1];
 f[4][2] = - temp*n_15[2];

 f[0][0] = - f[1][0] - f[2][0] - f[3][0] - f[4][0];
 f[0][1] = - f[1][1] - f[2][1] - f[3][1] - f[4][1];
 f[0][2] = - f[1][2] - f[2][2] - f[3][2] - f[4][2];
}

/*****************************************************************************
 * XY2Z2SymBendForce                                                         *
 *  Potential v(r) = 1/2 Sum (K_ij S_i S_j), i,j = 1, 3, 5, 7, 8             *
 *    (Wang et al. J.Phys.Chem. A101(1997)5919)				     *
 *****************************************************************************/
void XY2Z2SymBendForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0..5] = K_ij (in)     */
                                                /* par[6..8] = phi_0(in)     */
{
  vec_mt	n_12, n_13, n_14, n_15;
  double	r_abs12, r_abs13, r_abs14, r_abs15, 
		cos213, cos214, cos215, cos314, cos315, cos415,
		theta213, theta214, theta215, theta314, theta315, theta415,
		sym[5], fak[5], b[6], a[7], temp;

					/* Definition of the local coords   */
  da(r+1, r, (vec_mp)n_12);  		/* r_X(1)Y(2) 			    */
  da(r+2, r, (vec_mp)n_13);  		/* r_X(1)Y(3)			    */
  da(r+3, r, (vec_mp)n_14); 		/* r_X(1)Z(4)			    */
  da(r+4, r, (vec_mp)n_15);  		/* r_X(1)Z(5)			    */

  r_abs12 = ABS3(n_12);
  r_abs13 = ABS3(n_13);
  r_abs14 = ABS3(n_14);
  r_abs15 = ABS3(n_15);

  r_abs12 = 1.0/r_abs12;
  r_abs13 = 1.0/r_abs13;
  r_abs14 = 1.0/r_abs14;
  r_abs15 = 1.0/r_abs15;

  n_12[0] = n_12[0]*r_abs12; 		/* Unit-vector n_12 		    */
  n_12[1] = n_12[1]*r_abs12; 
  n_12[2] = n_12[2]*r_abs12; 
  n_13[0] = n_13[0]*r_abs13; 		/* Unit-vector n_13		    */
  n_13[1] = n_13[1]*r_abs13; 
  n_13[2] = n_13[2]*r_abs13; 
  n_14[0] = n_14[0]*r_abs14; 		/* Unit-vector n_14 		    */
  n_14[1] = n_14[1]*r_abs14; 
  n_14[2] = n_14[2]*r_abs14; 
  n_15[0] = n_15[0]*r_abs15; 		/* Unit-vector n_15		    */
  n_15[1] = n_15[1]*r_abs15; 
  n_15[2] = n_15[2]*r_abs15; 

 cos213 = SCALPROD(n_12,n_13);
 cos214 = SCALPROD(n_12,n_14);
 cos215 = SCALPROD(n_12,n_15);
 cos314 = SCALPROD(n_13,n_14);
 cos315 = SCALPROD(n_13,n_15);
 cos415 = SCALPROD(n_14,n_15);

 theta213 = acos(cos213); /* Y(2)-X(1)-Y(3) */ 
 theta214 = acos(cos214); /* Y(2)-X(1)-Z(4) */ 
 theta215 = acos(cos215); /* Y(2)-X(1)-Z(5) */ 
 theta314 = acos(cos314); /* Y(3)-X(1)-Z(4) */ 
 theta315 = acos(cos315); /* Y(3)-X(1)-Z(5) */ 
 theta415 = acos(cos415); /* Z(4)-X(1)-Z(5) */ 

/* Wang */
 a[0] = 0.408248;	a[1] = 3*0.408248;	a[2] = 0.408248;
 a[3] = 3*0.204124;	a[4] = 0.204124;	a[5] = 0.204124;
 a[6] = 0.5;
/* Davis */ /*
 a[0] = 0.0;	a[1] = 0.8940;	a[2] = 0.2240;
 a[3] = 0.9002;	a[4] = 0.0;	a[5] = 0.2177;
 a[6] = - 0.5;*/
					/* Symmetry coordinates 	    */
 sym[0] = a[0]*(theta213 -  par[6]) + a[1]*(theta415 - par[7]) 
        - a[2]*(theta214 + theta215 + theta314 + theta315 - 4*par[8]);
 sym[1] = a[3]*(theta213 -  par[6]) + a[4]*(theta415 - par[7]) 
        - a[5]*(theta214 + theta215 + theta314 + theta315 - 4*par[8]);
 sym[2] = a[6]*(theta214 - theta215 - theta314 + theta315); 
 sym[3] = a[6]*(theta214 - theta215 + theta314 - theta315); 
 sym[4] = a[6]*(theta214 + theta215 - theta314 - theta315);

 (*pot) += 0.5 * ( par[0]*SQR(sym[0]) + par[1]*SQR(sym[1])
		 + par[2]*SQR(sym[2]) + par[3]*SQR(sym[3]) 
		 + par[4]*SQR(sym[4]) ) + par[5]*sym[0]*sym[1];

 fak[0] =  par[0]*sym[0] + par[5]*sym[1]; 
 fak[1] =  par[5]*sym[0] + par[1]*sym[1]; 
 fak[2] =  par[2]*sym[2]; 
 fak[3] =  par[3]*sym[3]; 
 fak[4] =  par[4]*sym[4]; 

					  /* 1.0/sin(x) produziert bei HP's  */
  					  /* ab x < 1.0e-75 UNDERFLOW !      */
 temp = a[2]*fak[0] + a[5]*fak[1];
 b[0] = 1.0/sin(theta213)*(a[0]*fak[0] + a[3]*fak[1]); 
 b[1] = 1.0/sin(theta214)*((fak[2] + fak[3] + fak[4])*a[6] - temp);
 b[2] = 1.0/sin(theta215)*((fak[4] - fak[3] - fak[2])*a[6] - temp);
 b[3] = 1.0/sin(theta314)*((fak[3] - fak[2] - fak[4])*a[6] - temp);
 b[4] = 1.0/sin(theta315)*((fak[2] - fak[3] - fak[4])*a[6] - temp);
 b[5] = 1.0/sin(theta415)*(a[1]*fak[0] + a[4]*fak[1]);

 temp = b[0]*cos213 + b[1]*cos214 + b[2]*cos215;
 f[1][0] = (b[0]*n_13[0] + b[1]*n_14[0] + b[2]*n_15[0] - temp*n_12[0])*r_abs12;
 f[1][1] = (b[0]*n_13[1] + b[1]*n_14[1] + b[2]*n_15[1] - temp*n_12[1])*r_abs12;
 f[1][2] = (b[0]*n_13[2] + b[1]*n_14[2] + b[2]*n_15[2] - temp*n_12[2])*r_abs12;

 temp = b[0]*cos213 + b[3]*cos314 + b[4]*cos315;
 f[2][0] = (b[0]*n_12[0] + b[3]*n_14[0] + b[4]*n_15[0] - temp*n_13[0])*r_abs13;
 f[2][1] = (b[0]*n_12[1] + b[3]*n_14[1] + b[4]*n_15[1] - temp*n_13[1])*r_abs13;
 f[2][2] = (b[0]*n_12[2] + b[3]*n_14[2] + b[4]*n_15[2] - temp*n_13[2])*r_abs13;

 temp = b[1]*cos214 + b[3]*cos314 + b[5]*cos415;
 f[3][0] = (b[1]*n_12[0] + b[3]*n_13[0] + b[5]*n_15[0] - temp*n_14[0])*r_abs14;
 f[3][1] = (b[1]*n_12[1] + b[3]*n_13[1] + b[5]*n_15[1] - temp*n_14[1])*r_abs14;
 f[3][2] = (b[1]*n_12[2] + b[3]*n_13[2] + b[5]*n_15[2] - temp*n_14[2])*r_abs14;

 temp = b[2]*cos215 + b[4]*cos315 + b[5]*cos415;
 f[4][0] = (b[2]*n_12[0] + b[4]*n_13[0] + b[5]*n_14[0] - temp*n_15[0])*r_abs15;
 f[4][1] = (b[2]*n_12[1] + b[4]*n_13[1] + b[5]*n_14[1] - temp*n_15[1])*r_abs15;
 f[4][2] = (b[2]*n_12[2] + b[4]*n_13[2] + b[5]*n_14[2] - temp*n_15[2])*r_abs15;

 f[0][0] = - f[1][0] - f[2][0] - f[3][0] - f[4][0];
 f[0][1] = - f[1][1] - f[2][1] - f[3][1] - f[4][1];
 f[0][2] = - f[1][2] - f[2][2] - f[3][2] - f[4][2];

}

/*****************************************************************************
 * XY2Z2SymIntForce                                                          *
 *  Potential v(r) = 1/2 Sum (K_ij S_i S_j), ij = 12, 14, 23, 34, 67, 89     *
 *    (Wang et al. J.Phys.Chem. A101(1997)5919)				     *
 *****************************************************************************/
void XY2Z2SymIntForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0..5] = K_ij (in)     */
                                                /* par[6..7] = r_eq (in)     */
                                                /* par[8..10] = phi_0(in)    */
{
  vec_mt	n_12, n_13, n_14, n_15, d0, d1;
  double	r_abs12, r_abs13, r_abs14, r_abs15, 
		cos213, cos214, cos215, cos314, cos315, cos415,
		theta213, theta214, theta215, theta314, theta315, theta415,
		sym[8], fak[8], b[6], a[7], temp;

					/* Definition of the local coords   */
  da(r+1, r, (vec_mp)n_12);  		/* r_X(1)Y(2) 			    */
  da(r+2, r, (vec_mp)n_13);  		/* r_X(1)Y(3)			    */
  da(r+3, r, (vec_mp)n_14); 		/* r_X(1)Z(4)			    */
  da(r+4, r, (vec_mp)n_15);  		/* r_X(1)Z(5)			    */

  r_abs12 = ABS3(n_12);
  r_abs13 = ABS3(n_13);
  r_abs14 = ABS3(n_14);
  r_abs15 = ABS3(n_15);
					 /* Symmetry coords for stretches   */
 sym[1] = (r_abs14 + r_abs15 - 2*par[7])*0.707107;
 sym[3] = (r_abs12 + r_abs13 - 2*par[6])*0.707107;
 sym[4] = (r_abs14 - r_abs15)*0.707107;
 sym[7] = (r_abs12 - r_abs13)*0.707107;

  r_abs12 = 1.0/r_abs12;
  r_abs13 = 1.0/r_abs13;
  r_abs14 = 1.0/r_abs14;
  r_abs15 = 1.0/r_abs15;

  n_12[0] = n_12[0]*r_abs12; 		/* Unit-vector n_12 		    */
  n_12[1] = n_12[1]*r_abs12; 
  n_12[2] = n_12[2]*r_abs12; 
  n_13[0] = n_13[0]*r_abs13; 		/* Unit-vector n_13		    */
  n_13[1] = n_13[1]*r_abs13; 
  n_13[2] = n_13[2]*r_abs13; 
  n_14[0] = n_14[0]*r_abs14; 		/* Unit-vector n_14 		    */
  n_14[1] = n_14[1]*r_abs14; 
  n_14[2] = n_14[2]*r_abs14; 
  n_15[0] = n_15[0]*r_abs15; 		/* Unit-vector n_15		    */
  n_15[1] = n_15[1]*r_abs15; 
  n_15[2] = n_15[2]*r_abs15; 

 cos213 = SCALPROD(n_12,n_13);
 cos214 = SCALPROD(n_12,n_14);
 cos215 = SCALPROD(n_12,n_15);
 cos314 = SCALPROD(n_13,n_14);
 cos315 = SCALPROD(n_13,n_15);
 cos415 = SCALPROD(n_14,n_15);

 theta213 = acos(cos213); /* Y(2)-X(1)-Y(3) */ 
 theta214 = acos(cos214); /* Y(2)-X(1)-Z(4) */ 
 theta215 = acos(cos215); /* Y(2)-X(1)-Z(5) */ 
 theta314 = acos(cos314); /* Y(3)-X(1)-Z(4) */ 
 theta315 = acos(cos315); /* Y(3)-X(1)-Z(5) */ 
 theta415 = acos(cos415); /* Z(4)-X(1)-Z(5) */ 

/* Wang */
 a[0] = 0.408248;	a[1] = 3*0.408248;	a[2] = 0.408248;
 a[3] = 3*0.204124;	a[4] = 0.204124;	a[5] = 0.204124;
 a[6] = 0.5;
/* Davis */ /*
 a[0] = 0.0;	a[1] = 0.8940;	a[2] = 0.2240; 
 a[3] = 0.9002;	a[4] = 0.0;	a[5] = 0.2177;
 a[6] = - 0.5;*/
 					/* Symmetry coordinates for bends   */
 sym[0] = a[0]*(theta213 -  par[8]) + a[1]*(theta415 - par[9]) 
        - a[2]*(theta214 + theta215 + theta314 + theta315 - 4*par[10]);
 sym[2] = a[3]*(theta213 -  par[8]) + a[4]*(theta415 - par[9]) 
        - a[5]*(theta214 + theta215 + theta314 + theta315 - 4*par[10]);
 sym[5] = a[6]*(theta214 - theta215 + theta314 - theta315);
 sym[6] = a[6]*(theta214 + theta215 - theta314 - theta315);

 (*pot) += sym[0]*(par[0]*sym[1] + par[1]*sym[3]) 
         + sym[2]*(par[2]*sym[1] + par[3]*sym[3]) 
	 + sym[4]*par[4]*sym[5] + sym[6]*par[5]*sym[7];
 
 fak[0] = par[0]*sym[1] + par[1]*sym[3]; 
 fak[1] = par[0]*sym[0] + par[2]*sym[2]; 
 fak[2] = par[2]*sym[1] + par[3]*sym[3]; 
 fak[3] = par[1]*sym[0] + par[3]*sym[2]; 
 fak[4] = par[4]*sym[5]; 
 fak[5] = par[4]*sym[4]; 
 fak[6] = par[5]*sym[7]; 
 fak[7] = par[5]*sym[6]; 

					  /* 1.0/sin(x) produziert bei HP's  */
  					  /* ab x < 1.0e-75 UNDERFLOW !      */
 temp = a[2]*fak[0] + a[5]*fak[2];
 b[0] =   1.0/sin(theta213)*(a[0]*fak[0] + a[3]*fak[2]);
 b[1] =   1.0/sin(theta214)*((fak[5] + fak[6])*a[6] - temp);
 b[2] =   1.0/sin(theta215)*((fak[6] - fak[5])*a[6] - temp);
 b[3] =   1.0/sin(theta314)*((fak[5] - fak[6])*a[6] - temp);
 b[4] = - 1.0/sin(theta315)*((fak[5] + fak[6])*a[6] + temp);
 b[5] =   1.0/sin(theta415)*(a[1]*fak[0] + a[4]*fak[2]);

 temp = b[0]*cos213 + b[1]*cos214 + b[2]*cos215;
 f[1][0] = (b[0]*n_13[0] + b[1]*n_14[0] + b[2]*n_15[0] - temp*n_12[0])*r_abs12
	   - (fak[3] + fak[7])*n_12[0]*0.707107;
 f[1][1] = (b[0]*n_13[1] + b[1]*n_14[1] + b[2]*n_15[1] - temp*n_12[1])*r_abs12
	   - (fak[3] + fak[7])*n_12[1]*0.707107;
 f[1][2] = (b[0]*n_13[2] + b[1]*n_14[2] + b[2]*n_15[2] - temp*n_12[2])*r_abs12
	   - (fak[3] + fak[7])*n_12[2]*0.707107;

 temp = b[0]*cos213 + b[3]*cos314 + b[4]*cos315;
 f[2][0] = (b[0]*n_12[0] + b[3]*n_14[0] + b[4]*n_15[0] - temp*n_13[0])*r_abs13
	   - (fak[3] - fak[7])*n_13[0]*0.707107;
 f[2][1] = (b[0]*n_12[1] + b[3]*n_14[1] + b[4]*n_15[1] - temp*n_13[1])*r_abs13
	   - (fak[3] - fak[7])*n_13[1]*0.707107;
 f[2][2] = (b[0]*n_12[2] + b[3]*n_14[2] + b[4]*n_15[2] - temp*n_13[2])*r_abs13
	   - (fak[3] - fak[7])*n_13[2]*0.707107;

 temp = b[1]*cos214 + b[3]*cos314 + b[5]*cos415;
 f[3][0] = (b[1]*n_12[0] + b[3]*n_13[0] + b[5]*n_15[0] - temp*n_14[0])*r_abs14
	   - (fak[1] + fak[4])*n_14[0]*0.707107;
 f[3][1] = (b[1]*n_12[1] + b[3]*n_13[1] + b[5]*n_15[1] - temp*n_14[1])*r_abs14
	   - (fak[1] + fak[4])*n_14[1]*0.707107;
 f[3][2] = (b[1]*n_12[2] + b[3]*n_13[2] + b[5]*n_15[2] - temp*n_14[2])*r_abs14
	   - (fak[1] + fak[4])*n_14[2]*0.707107;

 temp = b[2]*cos215 + b[4]*cos315 + b[5]*cos415;
 f[4][0] = (b[2]*n_12[0] + b[4]*n_13[0] + b[5]*n_14[0] - temp*n_15[0])*r_abs15
	   - (fak[1] - fak[4])*n_15[0]*0.707107;
 f[4][1] = (b[2]*n_12[1] + b[4]*n_13[1] + b[5]*n_14[1] - temp*n_15[1])*r_abs15
	   - (fak[1] - fak[4])*n_15[1]*0.707107;
 f[4][2] = (b[2]*n_12[2] + b[4]*n_13[2] + b[5]*n_14[2] - temp*n_15[2])*r_abs15
	   - (fak[1] - fak[4])*n_15[2]*0.707107;

 f[0][0] = - f[1][0] - f[2][0] - f[3][0] - f[4][0];
 f[0][1] = - f[1][1] - f[2][1] - f[3][1] - f[4][1];
 f[0][2] = - f[1][2] - f[2][2] - f[3][2] - f[4][2];
}


/*****************************************************************************
 * XY2Z2Force                                                         *
 *  Potential v(r) = 1/2 Sum (K_ij S_i S_j), i,j = 1, 3, 5, 7, 8             *
 *    (Wang et al. J.Phys.Chem. A101(1997)5919)				     *
 *****************************************************************************/
void XY2Z2Force(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
{
  vec_mt	n_12, n_13, n_14, n_15;
  double	r_abs12, r_abs13, r_abs14, r_abs15, r0_CH, r0_CCl, 
		cos213, cos214, cos215, cos314, cos315, cos415,
		theta213, theta214, theta215, theta314, theta315, theta415,
		thetaHCH_0, thetaHCCl_0, thetaClCCl_0,
		sym[9], fak[9], b[10], a[7], temp, para[22];

 if (par[0] > 0)
  {
   para[0] = 632.927*100; para[1] = 2297.45*100; para[2] = 351.091*100;
   para[3] = 3064.07*100; para[4] = 357.113*100; para[5] = 1764.49*100;
   para[6] = 397.461*100; para[7] = 495.622*100; para[8] = 3034.56*100;

    para[9] =   148.145*100;
   para[10] = - 137.907*100;
   para[11] =         0*100;
   para[12] = - 133.089*100;
   para[13] =         0*100;
   para[14] =         0*100;
   para[15] = - 333.024*100;
   para[16] =         0*100;
  }
 else
  {
   para[0] = 632.927*100*169; para[1] = 2297.45*100*64; para[2] = 351.091*100*36;
   para[3] = 3064.07*100*4; para[4] = 357.113*100*16; para[5] = 1764.49*100*32;
   para[6] = 397.461*100*16; para[7] = 495.622*100*16; para[8] = 3034.56*100;

    para[9] = 0;
   para[10] = 0;
   para[11] = 0;
   para[12] = 0;
   para[13] = 0;
   para[14] = 0;
   para[15] = 0;
   para[16] = 0;
  }
 para[17] = 1.075; para[18] = 1.770;
 para[19] = 1.92231; para[20] = 1.89665; para[21] = 1.95512;

 a[0] = 0; a[1] = 0.8940; a[2] = 0.2240;
 a[3] = 0.9002; a[4] = 0; a[5] = 0.2177;
 a[6] = -1.0;

					/* Definition of the local coords   */
  da(r+1, r, (vec_mp)n_12);  		/* r_X(1)Y(2) 			    */
  da(r+2, r, (vec_mp)n_13);  		/* r_X(1)Y(3)			    */
  da(r+3, r, (vec_mp)n_14); 		/* r_X(1)Z(4)			    */
  da(r+4, r, (vec_mp)n_15);  		/* r_X(1)Z(5)			    */

  r_abs12 = ABS3(n_12);
  r_abs13 = ABS3(n_13);
  r_abs14 = ABS3(n_14);
  r_abs15 = ABS3(n_15);

  sym[1] = (r_abs14 + r_abs15 - 2*para[18])/1.41421;
  sym[3] = (r_abs12 + r_abs13 - 2*para[17])/1.41421;
  sym[5] = (r_abs14 - r_abs15)/1.41421;
  sym[8] = (r_abs12 - r_abs13)/1.41421;

  r_abs12 = 1.0/r_abs12;
  r_abs13 = 1.0/r_abs13;
  r_abs14 = 1.0/r_abs14;
  r_abs15 = 1.0/r_abs15;

  n_12[0] = n_12[0]*r_abs12; 		/* Unit-vector n_12 		    */
  n_12[1] = n_12[1]*r_abs12; 
  n_12[2] = n_12[2]*r_abs12; 
  n_13[0] = n_13[0]*r_abs13; 		/* Unit-vector n_13		    */
  n_13[1] = n_13[1]*r_abs13; 
  n_13[2] = n_13[2]*r_abs13; 
  n_14[0] = n_14[0]*r_abs14; 		/* Unit-vector n_14 		    */
  n_14[1] = n_14[1]*r_abs14; 
  n_14[2] = n_14[2]*r_abs14; 
  n_15[0] = n_15[0]*r_abs15; 		/* Unit-vector n_15		    */
  n_15[1] = n_15[1]*r_abs15; 
  n_15[2] = n_15[2]*r_abs15; 

 cos213 = SCALPROD(n_12,n_13);
 cos214 = SCALPROD(n_12,n_14);
 cos215 = SCALPROD(n_12,n_15);
 cos314 = SCALPROD(n_13,n_14);
 cos315 = SCALPROD(n_13,n_15);
 cos415 = SCALPROD(n_14,n_15);

 theta213 = acos(cos213); /* Y(2)-X(1)-Y(3) */ 
 theta214 = acos(cos214); /* Y(2)-X(1)-Z(4) */ 
 theta215 = acos(cos215); /* Y(2)-X(1)-Z(5) */ 
 theta314 = acos(cos314); /* Y(3)-X(1)-Z(4) */ 
 theta315 = acos(cos315); /* Y(3)-X(1)-Z(5) */ 
 theta415 = acos(cos415); /* Z(4)-X(1)-Z(5) */ 

					/* Symmetry coordinates 	    */
 sym[0] = a[0]*(theta213 - para[19]) + a[1]*(theta415 - para[21]) 
        - a[2]*(theta214 + theta215 + theta314 + theta315 - 4*para[20]);
 sym[2] = a[3]*(theta213 - para[19]) + a[4]*(theta415 - para[21]) 
        - a[5]*(theta214 + theta215 + theta314 + theta315 - 4*para[20]);
 sym[4] = a[6]*(theta214 - theta215 - theta314 + theta315)/2; 
 sym[6] = a[6]*(theta214 - theta215 + theta314 - theta315)/2; 
 sym[7] = a[6]*(theta214 + theta215 - theta314 - theta315)/2;

 (*pot) += 0.5 * ( para[0]*sym[0]*sym[0] + para[1]*sym[1]*sym[1]
		 + para[2]*sym[2]*sym[2] + para[3]*sym[3]*sym[3] 
		 + para[4]*sym[4]*sym[4] + para[5]*sym[5]*sym[5] 
		 + para[6]*sym[6]*sym[6] + para[7]*sym[7]*sym[7] 
		 + para[8]*sym[8]*sym[8] ) 
		 +  para[9]*sym[0]*sym[1] + para[10]*sym[0]*sym[2] 
		 + para[11]*sym[0]*sym[3] + para[12]*sym[1]*sym[2] 
		 + para[13]*sym[1]*sym[3] + para[14]*sym[2]*sym[3]
		 + para[15]*sym[5]*sym[6] + para[16]*sym[7]*sym[8];

 fak[0] =   para[0]*sym[0] + para[9] *sym[1] + para[10]*sym[2] + para[11]*sym[3]; 
 fak[1] =   para[9]*sym[0] + para[1] *sym[1] + para[12]*sym[2] + para[13]*sym[3]; 
 fak[2] =  para[10]*sym[0] + para[12]*sym[1] +  para[2]*sym[2] + para[14]*sym[3]; 
 fak[3] =  para[11]*sym[3] + para[13]*sym[3] + para[14]*sym[3] +  para[3]*sym[3]; 
 fak[4] =  para[4]*sym[4]; 
 fak[5] =   para[5]*sym[5] + para[15]*sym[6]; 
 fak[6] =  para[15]*sym[5] +  para[6]*sym[6]; 
 fak[7] =   para[7]*sym[7] + para[16]*sym[8]; 
 fak[8] =  para[16]*sym[7] +  para[8]*sym[8]; 

					  /* 1.0/sin(x) produziert bei HP's  */
  					  /* ab x < 1.0e-75 UNDERFLOW !      */
 b[0] = (fak[3] + fak[8])/1.41421;
 b[1] = (fak[3] - fak[8])/1.41421;
 b[2] = (fak[1] + fak[5])/1.41421;
 b[3] = (fak[1] - fak[5])/1.41421;
 b[4] =   1.0/sin(theta213)*(a[0]*fak[0] + a[3]*fak[2]);
 b[5] =   1.0/sin(theta214)*(( fak[4] + fak[6] + fak[7])*a[6]/2 - a[2]*fak[0] - a[5]*fak[2]);
 b[6] =   1.0/sin(theta215)*((-fak[4] - fak[6] + fak[7])*a[6]/2 - a[2]*fak[0] - a[5]*fak[2]);
 b[7] =   1.0/sin(theta314)*((-fak[4] + fak[6] - fak[7])*a[6]/2 - a[2]*fak[0] - a[5]*fak[2]);
 b[8] =   1.0/sin(theta315)*(( fak[4] - fak[6] - fak[7])*a[6]/2 - a[2]*fak[0] - a[5]*fak[2]);
 b[9] =   1.0/sin(theta415)*(a[1]*fak[0] + a[4]*fak[2]);

 f[1][0] = -(b[0] + (b[4]*cos213 + b[5]*cos214 + b[6]*cos215)*r_abs12)*n_12[0] 
	   +(b[4]*n_13[0] + b[5]*n_14[0] + b[6]*n_15[0])*r_abs12;
 f[1][1] = -(b[0] + (b[4]*cos213 + b[5]*cos214 + b[6]*cos215)*r_abs12)*n_12[1] 
	   +(b[4]*n_13[1] + b[5]*n_14[1] + b[6]*n_15[1])*r_abs12;
 f[1][2] = -(b[0] + (b[4]*cos213 + b[5]*cos214 + b[6]*cos215)*r_abs12)*n_12[2] 
	   +(b[4]*n_13[2] + b[5]*n_14[2] + b[6]*n_15[2])*r_abs12;

 f[2][0] = -(b[1] + (b[4]*cos213 + b[7]*cos314 + b[8]*cos315)*r_abs13)*n_13[0] 
	   +(b[4]*n_12[0] + b[7]*n_14[0] + b[8]*n_15[0])*r_abs13;
 f[2][1] = -(b[1] + (b[4]*cos213 + b[7]*cos314 + b[8]*cos315)*r_abs13)*n_13[1] 
	   +(b[4]*n_12[1] + b[7]*n_14[1] + b[8]*n_15[1])*r_abs13;
 f[2][2] = -(b[1] + (b[4]*cos213 + b[7]*cos314 + b[8]*cos315)*r_abs13)*n_13[2] 
	   +(b[4]*n_12[2] + b[7]*n_14[2] + b[8]*n_15[2])*r_abs13;

 f[3][0] = -(b[2] + (b[5]*cos214 + b[7]*cos314 + b[9]*cos415)*r_abs14)*n_14[0] 
	   +(b[5]*n_12[0] + b[7]*n_13[0] + b[9]*n_15[0])*r_abs14;
 f[3][1] = -(b[2] + (b[5]*cos214 + b[7]*cos314 + b[9]*cos415)*r_abs14)*n_14[1] 
	   +(b[5]*n_12[1] + b[7]*n_13[1] + b[9]*n_15[1])*r_abs14;
 f[3][2] = -(b[2] + (b[5]*cos214 + b[7]*cos314 + b[9]*cos415)*r_abs14)*n_14[2] 
	   +(b[5]*n_12[2] + b[7]*n_13[2] + b[9]*n_15[2])*r_abs14;

 f[4][0] = -(b[3] + (b[6]*cos215 + b[8]*cos315 + b[9]*cos415)*r_abs15)*n_15[0] 
	   +(b[6]*n_12[0] + b[8]*n_13[0] + b[9]*n_14[0])*r_abs15;
 f[4][1] = -(b[3] + (b[6]*cos215 + b[8]*cos315 + b[9]*cos415)*r_abs15)*n_15[1] 
	   +(b[6]*n_12[1] + b[8]*n_13[1] + b[9]*n_14[1])*r_abs15;
 f[4][2] = -(b[3] + (b[6]*cos215 + b[8]*cos315 + b[9]*cos415)*r_abs15)*n_15[2] 
	   +(b[6]*n_12[2] + b[8]*n_13[2] + b[9]*n_14[2])*r_abs15;

 f[0][0] = - f[1][0] - f[2][0] - f[3][0] - f[4][0];
 f[0][1] = - f[1][1] - f[2][1] - f[3][1] - f[4][1];
 f[0][2] = - f[1][2] - f[2][2] - f[3][2] - f[4][2];

}

/*****************************************************************************
 * XY4Force                                                          *
 *  Potential v(r) = 1/2 Sum (K_ij S_i S_j), ij = 12, 14, 23, 34, 67, 89     *
 *    (Shimanouchi et al. J.Mol.Spectr. 19(1966)78)				     *
 *****************************************************************************/
void XY4Force(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0..4] = K_ij (in)     */
                                                /* par[5]    = r_eq (in)     */
                                                /* par[6]    = phi_0(in)     */
{
  vec_mt	n_12, n_13, n_14, n_15, d0, d1;
  double	r_abs12, r_abs13, r_abs14, r_abs15, 
		cos213, cos214, cos215, cos314, cos315, cos415,
		theta213, theta214, theta215, theta314, theta315, theta415,
		sym[9], fak[9], b[6], a[7], temp;

					/* Definition of the local coords   */
  da(r+1, r, (vec_mp)n_12);  		/* r_X(1)Y(2) 			    */
  da(r+2, r, (vec_mp)n_13);  		/* r_X(1)Y(3)			    */
  da(r+3, r, (vec_mp)n_14); 		/* r_X(1)Y(4)			    */
  da(r+4, r, (vec_mp)n_15);  		/* r_X(1)Y(5)			    */

  r_abs12 = ABS3(n_12);
  r_abs13 = ABS3(n_13);
  r_abs14 = ABS3(n_14);
  r_abs15 = ABS3(n_15);
					 /* Symmetry coords for stretches   */
 sym[0] = (r_abs12 + r_abs13 + r_abs14 + r_abs15 - 4*par[5])*0.5;
 sym[3] = (3*r_abs12 - r_abs13 - r_abs14 - r_abs15)*0.288675;
 sym[5] = (2*r_abs13 - r_abs14 - r_abs15)*0.408248;
 sym[7] = (r_abs14 - r_abs15)*0.707107;

  r_abs12 = 1.0/r_abs12;
  r_abs13 = 1.0/r_abs13;
  r_abs14 = 1.0/r_abs14;
  r_abs15 = 1.0/r_abs15;

  n_12[0] = n_12[0]*r_abs12; 		/* Unit-vector n_12 		    */
  n_12[1] = n_12[1]*r_abs12; 
  n_12[2] = n_12[2]*r_abs12; 
  n_13[0] = n_13[0]*r_abs13; 		/* Unit-vector n_13		    */
  n_13[1] = n_13[1]*r_abs13; 
  n_13[2] = n_13[2]*r_abs13; 
  n_14[0] = n_14[0]*r_abs14; 		/* Unit-vector n_14 		    */
  n_14[1] = n_14[1]*r_abs14; 
  n_14[2] = n_14[2]*r_abs14; 
  n_15[0] = n_15[0]*r_abs15; 		/* Unit-vector n_15		    */
  n_15[1] = n_15[1]*r_abs15; 
  n_15[2] = n_15[2]*r_abs15; 

 cos213 = SCALPROD(n_12,n_13);
 cos214 = SCALPROD(n_12,n_14);
 cos215 = SCALPROD(n_12,n_15);
 cos314 = SCALPROD(n_13,n_14);
 cos315 = SCALPROD(n_13,n_15);
 cos415 = SCALPROD(n_14,n_15);

 theta213 = acos(cos213); /* Y(2)-X(1)-Y(3) */ 
 theta214 = acos(cos214); /* Y(2)-X(1)-Y(4) */ 
 theta215 = acos(cos215); /* Y(2)-X(1)-Y(5) */ 
 theta314 = acos(cos314); /* Y(3)-X(1)-Y(4) */ 
 theta315 = acos(cos315); /* Y(3)-X(1)-Y(5) */ 
 theta415 = acos(cos415); /* Y(4)-X(1)-Y(5) */ 

 					/* Symmetry coordinates for bends   */
 sym[1] = (2*theta213 - theta214 - theta215 
		- theta314 - theta315 + 2*theta415)*0.288675;
 sym[2] = (theta214 - theta215 - theta314 + theta315)*0.5;
 sym[4] = (theta213 + theta214 + theta215 
		- theta314 - theta315 - theta415)*0.408248;
 sym[6] = (2*theta213 - theta214 - theta215 
		+ theta314 + theta315 - 2*theta415)*0.288675;
 sym[8] = (theta214 - theta215 + theta314 - theta315)*0.5;

 (*pot) += 0.5*(par[0]*sym[0]*sym[0] + par[1]*(sym[1]*sym[1] + sym[2]*sym[2]) 
	 + par[2]*(sym[3]*sym[3] + sym[5]*sym[5] + sym[7]*sym[7])
	 + par[3]*(sym[4]*sym[4] + sym[6]*sym[6] + sym[8]*sym[8]))
	 + par[4]*(sym[3]*sym[4] + sym[5]*sym[6] + sym[7]*sym[8]);
 
 fak[0] = par[0]*sym[0]; 
 fak[1] = par[1]*sym[1]; 
 fak[2] = par[1]*sym[2]; 
 fak[3] = par[2]*sym[3] + par[4]*sym[4]; 
 fak[4] = par[4]*sym[3] + par[3]*sym[4]; 
 fak[5] = par[2]*sym[5] + par[4]*sym[6]; 
 fak[6] = par[4]*sym[5] + par[3]*sym[6]; 
 fak[7] = par[2]*sym[7] + par[4]*sym[8]; 
 fak[8] = par[4]*sym[7] + par[3]*sym[8]; 

					  /* 1.0/sin(x) produziert bei HP's  */
  					  /* ab x < 1.0e-75 UNDERFLOW !      */
 b[0] =   1.0/sin(theta213)*(fak[1] + fak[4]*0.707107 + fak[6])*0.57735;
 b[1] =   1.0/sin(theta214)*((fak[2] + fak[8])*0.5 
			      - (fak[1] - fak[4]*1.41421 + fak[6])*0.288675);
 b[2] = - 1.0/sin(theta215)*((fak[2] + fak[8])*0.5 
			      + (fak[1] - fak[4]*1.41421 + fak[6])*0.288675);
 b[3] = - 1.0/sin(theta314)*((fak[2] - fak[8])*0.5 
			      + (fak[1] + fak[4]*1.41421 - fak[6])*0.288675);
 b[4] =   1.0/sin(theta315)*((fak[2] - fak[8])*0.5 
			      - (fak[1] + fak[4]*1.41421 - fak[6])*0.288675);
 b[5] =   1.0/sin(theta415)*(fak[1] - fak[4]*0.707107 - fak[6])*0.57735;
 
 temp = (fak[0] + fak[3]*1.73205)*0.5 
			+ (b[0]*cos213 + b[1]*cos214 + b[2]*cos215)*r_abs12;
 f[1][0] = (b[0]*n_13[0] + b[1]*n_14[0] + b[2]*n_15[0])*r_abs12 - temp*n_12[0];
 f[1][1] = (b[0]*n_13[1] + b[1]*n_14[1] + b[2]*n_15[1])*r_abs12 - temp*n_12[1];
 f[1][2] = (b[0]*n_13[2] + b[1]*n_14[2] + b[2]*n_15[2])*r_abs12 - temp*n_12[2];

 temp = fak[0]*0.5 - fak[3]*0.288675 + fak[5]*0.816497
			+ (b[0]*cos213 + b[3]*cos314 + b[4]*cos315)*r_abs13;
 f[2][0] = (b[0]*n_12[0] + b[3]*n_14[0] + b[4]*n_15[0])*r_abs13 - temp*n_13[0];
 f[2][1] = (b[0]*n_12[1] + b[3]*n_14[1] + b[4]*n_15[1])*r_abs13 - temp*n_13[1];
 f[2][2] = (b[0]*n_12[2] + b[3]*n_14[2] + b[4]*n_15[2])*r_abs13 - temp*n_13[2];

 temp = fak[0]*0.5 - fak[3]*0.288675 - fak[5]*0.408248 + fak[7]*0.707107
			+ (b[1]*cos214 + b[3]*cos314 + b[5]*cos415)*r_abs14;
 f[3][0] = (b[1]*n_12[0] + b[3]*n_13[0] + b[5]*n_15[0])*r_abs14 - temp*n_14[0];
 f[3][1] = (b[1]*n_12[1] + b[3]*n_13[1] + b[5]*n_15[1])*r_abs14 - temp*n_14[1];
 f[3][2] = (b[1]*n_12[2] + b[3]*n_13[2] + b[5]*n_15[2])*r_abs14 - temp*n_14[2];

 temp = fak[0]*0.5 - fak[3]*0.288675 - fak[5]*0.408248 - fak[7]*0.707107
			+ (b[2]*cos215 + b[4]*cos315 + b[5]*cos415)*r_abs15;
 f[4][0] = (b[2]*n_12[0] + b[4]*n_13[0] + b[5]*n_14[0])*r_abs15 - temp*n_15[0];
 f[4][1] = (b[2]*n_12[1] + b[4]*n_13[1] + b[5]*n_14[1])*r_abs15 - temp*n_15[1];
 f[4][2] = (b[2]*n_12[2] + b[4]*n_13[2] + b[5]*n_14[2])*r_abs15 - temp*n_15[2];

 f[0][0] = - f[1][0] - f[2][0] - f[3][0] - f[4][0];
 f[0][1] = - f[1][1] - f[2][1] - f[3][1] - f[4][1];
 f[0][2] = - f[1][2] - f[2][2] - f[3][2] - f[4][2];
 }

/*****************************************************************************
 * NoBondForce  					                     *
 *  Potential v(r) = par[0] * ( 2 ( par[1] / r )^9 - 3 ( par[1] / r )^6 )    *
 *****************************************************************************/
void NoBondForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = epsilon   (in)   */
                                                /* par[1] = r^*       (in)   */
{
  vec_mt		r_diff;
  double		r_abs, fak1;
  double                po1, po2;

  r_diff[0] = r[1][0] - r[0][0];
  r_diff[1] = r[1][1] - r[0][1];
  r_diff[2] = r[1][2] - r[0][2];

  r_abs = ABS3(r_diff);
  
  po1 = pow(par[1]/r_abs,9.0);
  po2 = pow(par[1]/r_abs,6.0);

  fak1 = par[0] * ( - 2.0 * 9.0 * po1 + 3.0 * 6.0 * po2 ) / SQR(r_abs);
  
  f[0][0] = r_diff[0] * fak1;
  f[0][1] = r_diff[1] * fak1;
  f[0][2] = r_diff[2] * fak1;
  f[1][0] = -r_diff[0] * fak1;
  f[1][1] = -r_diff[1] * fak1;
  f[1][2] = -r_diff[2] * fak1;
  
  (*pot) += par[0] * ( 2.0 * po1 - 3.0 * po2 );

}

/*****************************************************************************
 * Chedin2Force (quadratic terms in potential, 3 parameters)  		     *
 *  Potential: see 							     *
 *    Chedin A. and Teffo J.-L. J. Mol. Spectr. 107 (1984) 333  	     *
 *    Simons G., Parr R.G. and Finlan J.M. J. Chem. Phys. 59 (1973) 3229     *
 *    Carney G.D., Sprandel L.L. and Kern C.W. Adv. Chem. Phys. 37 (1977) 305*
 *****************************************************************************/
void Chedin2Force(r, f, par, pot)
 vec_mp         r,                              /* Particles coordinates (in)*/
                f;                              /* Forces (out)              */
 double         *pot,                           /* Energy (out)              */
                par[];                          /* Potential parameters (in) */
						/* par[0..2] = 11, 22, 33    */
{
  vec_mt          force1s, force3s,	    /* Stretch forces 		     */
  		  force1t, force3t,         /* Angle (torque) forces         */
                  vr21,			    /* Vector r_21 = r(O1)-r(C)      */
                  vr23, 		    /* Vector r_23 = r(O2)-r(C)      */
                  n21,			    /* Unit vector in direction r_21 */
                  n23,			    /* Unit vector in direction r_23 */
                  ncross; 		    /* Cross product of n_21 and n_23*/
  double          r21,r23, nabs,		    /* Distances between atoms       */
                  alpha, sinr23, cosr23,/* Angle between n_21 and n_23   */	
                  xi1,xi2,xi3,		    /* Intrisic SPF coordinates      */
                  dudrq1,dudrq2,dudrq3,	    /* Derivatives of potential on ic*/
                  rq1,rq2,rq3,		    /* Quasi-normal SPF coordinates  */
                  dwdrq1,dwdrq2,dwdrq3,	    /* Derivatives of potential on qn*/
                  re, a1, a2, a3;	    /* Equilibrium C-O distance and  */
                  			    /* parameters of transition to   */ 
                  	 		    /* quasi-normal coordinates      */ 
 re = 1.161; a1 = 20.294; a2 = 7.665; a3 = 14.248;

  vr21[0] = r[1][0] - r[0][0];				
  vr21[1] = r[1][1] - r[0][1];
  vr21[2] = r[1][2] - r[0][2];

  r21 = ABS3(vr21);

  vr23[0] = r[2][0] - r[0][0];				
  vr23[1] = r[2][1] - r[0][1];
  vr23[2] = r[2][2] - r[0][2];

  r23 = ABS3(vr23);

  KreuzProd(vr21,vr23,ncross);
  nabs = ABS3(ncross);
  alpha = asin(nabs/(r21*r23));

  xi1 = 1 - re/r21;
  xi2 = alpha;
  xi3 = 1 - re/r23;

  rq1 = a1*(xi1+xi3);
  rq2 = a2*xi2;
  rq3 = a3*(xi1-xi3);

  *pot += par[0]*rq1*rq1 + par[1]*rq2*rq2 + par[2]*rq3*rq3;

  dudrq1 = 2*par[0]*rq1;
  dudrq2 = 2*par[1]*rq2;
  dudrq3 = 2*par[2]*rq3;
					     /* Forces (0 is number for C) */

  f[1][0] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[0]
  	    -dudrq2*a2/(nabs/r23)*(vr23[0]/r23-vr21[0]*cos(PI-alpha)/r21);
  f[1][1] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[1]
  	    -dudrq2*a2/(nabs/r23)*(vr23[1]/r23-vr21[1]*cos(PI-alpha)/r21);
  f[1][2] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[2]
  	    -dudrq2*a2/(nabs/r23)*(vr23[2]/r23-vr21[2]*cos(PI-alpha)/r21);
  
  f[2][0] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[0]
  	    -dudrq2*a2/(nabs/r21)*(vr21[0]/r21-vr23[0]*cos(PI-alpha)/r23);
  f[2][1] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[1]
  	    -dudrq2*a2/(nabs/r21)*(vr21[1]/r21-vr23[1]*cos(PI-alpha)/r23);
  f[2][2] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[2]
  	    -dudrq2*a2/(nabs/r21)*(vr21[2]/r21-vr23[2]*cos(PI-alpha)/r23);
    
  f[0][0] = - f[1][0] - f[2][0];
  f[0][1] = - f[1][1] - f[2][1];
  f[0][2] = - f[1][2] - f[2][2];
}

/*****************************************************************************
 * Chedin34Force (cubic & 4th order terms in potential, 3 + 6 = 9 parameters)*
 *  Potential: see Chedin A. and Teffo J.-L. J. Mol. Spectr. 107 (1984) 333  *
 *    Chedin A. and Teffo J.-L. J. Mol. Spectr. 107 (1984) 333  	     *
 *    Simons G., Parr R.G. and Finlan J.M. J. Chem. Phys. 59 (1973) 3229     *
 *    Carney G.D., Sprandel L.L. and Kern C.W. Adv. Chem. Phys. 37 (1977) 305*
 *****************************************************************************/
void Chedin34Force(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte (out)             */
 double         *pot,                           /* Energie (out)             */
                par[];                          /* Potentialparameter (in)   */
						/* par[0..2] = 111, 122, 133 */
						/* par[3..8] =1111,1122,1133 */
						/* 	      2222,2233,3333 */
{
  vec_mt          force1s, force3s,	    /* Stretch forces 		     */
  		  force1t, force3t,         /* Angle (torque) forces         */
                  vr21,			    /* Vector r_21 = r(O1)-r(C)      */
                  vr23, 		    /* Vector r_23 = r(O2)-r(C)      */
                  n21,			    /* Unit vector in direction r_21 */
                  n23,			    /* Unit vector in direction r_23 */
                  ncross; 		    /* Cross product of n_21 and n_23*/
  double          r21,r23, nabs,		    /* Distances between atoms       */
                  alpha, sinr23, cosr23,/* Angle between n_21 and n_23   */	
                  xi1,xi2,xi3,		    /* Intrisic SPF coordinates      */
                  dudrq1,dudrq2,dudrq3,	    /* Derivatives of potential on ic*/
                  rq1,rq2,rq3,		    /* Quasi-normal SPF coordinates  */
                  dwdrq1,dwdrq2,dwdrq3,	    /* Derivatives of potential on qn*/
                  re, a1, a2, a3;	    /* Equilibrium C-O distance and  */
                  			    /* parameters of transition to   */ 
                  	 		    /* quasi-normal coordinates      */ 
 re = 1.161; a1 = 20.294; a2 = 7.665; a3 = 14.248;

  vr21[0] = r[1][0] - r[0][0];				
  vr21[1] = r[1][1] - r[0][1];
  vr21[2] = r[1][2] - r[0][2];

  r21 = ABS3(vr21);

  vr23[0] = r[2][0] - r[0][0];				
  vr23[1] = r[2][1] - r[0][1];
  vr23[2] = r[2][2] - r[0][2];

  r23 = ABS3(vr23);

  KreuzProd(vr21,vr23,ncross);
  nabs = ABS3(ncross);
  alpha = asin(nabs/(r21*r23));

  xi1 = 1 - re/r21;
  xi2 = alpha;
  xi3 = 1 - re/r23;

  rq1 = a1*(xi1+xi3);
  rq2 = a2*xi2;
  rq3 = a3*(xi1-xi3);

  *pot += par[0]*rq1*rq1*rq1 + par[1]*rq1*rq2*rq2 + par[2]*rq1*rq3*rq3
        + par[3]*rq1*rq1*rq1*rq1 + par[4]*rq1*rq1*rq2*rq2 
        + par[5]*rq1*rq1*rq3*rq3 + par[6]*rq2*rq2*rq2*rq2 
        + par[7]*rq2*rq2*rq3*rq3 + par[8]*rq3*rq3*rq3*rq3;

  dudrq1 = 3*par[0]*rq1*rq1 + 4*par[3]*rq1*rq1*rq1 
           + (par[1] + 2*par[4]*rq1)*rq2*rq2 + (par[2] + 2*par[5]*rq1)*rq3*rq3;
  dudrq2 = 2*rq2*(rq1*(par[1] + par[4]*rq1) + 2*par[6]*rq2*rq2 + par[7]*rq3*rq3);
  dudrq3 = 2*rq3*(rq1*(par[2] + par[5]*rq1) + par[7]*rq2*rq2 + 2*par[8]*rq3*rq3);		
					
					     /* Forces (0 is number for C) */

  f[1][0] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[0]
  	    -dudrq2*a2/(nabs/r23)*(vr23[0]/r23-vr21[0]*cos(PI-alpha)/r21);
  f[1][1] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[1]
  	    -dudrq2*a2/(nabs/r23)*(vr23[1]/r23-vr21[1]*cos(PI-alpha)/r21);
  f[1][2] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[2]
  	    -dudrq2*a2/(nabs/r23)*(vr23[2]/r23-vr21[2]*cos(PI-alpha)/r21);
  
  f[2][0] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[0]
  	    -dudrq2*a2/(nabs/r21)*(vr21[0]/r21-vr23[0]*cos(PI-alpha)/r23);
  f[2][1] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[1]
  	    -dudrq2*a2/(nabs/r21)*(vr21[1]/r21-vr23[1]*cos(PI-alpha)/r23);
  f[2][2] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[2]
  	    -dudrq2*a2/(nabs/r21)*(vr21[2]/r21-vr23[2]*cos(PI-alpha)/r23);
    
  f[0][0] = - f[1][0] - f[2][0];
  f[0][1] = - f[1][1] - f[2][1];
  f[0][2] = - f[1][2] - f[2][2];
}

/*****************************************************************************
 * Chedin5Force (5th order terms in potential, 6 parameters)		     *
 *  Potential: see Chedin A. and Teffo J.-L. J. Mol. Spectr. 107 (1984) 333  *
 *    Chedin A. and Teffo J.-L. J. Mol. Spectr. 107 (1984) 333  	     *
 *    Simons G., Parr R.G. and Finlan J.M. J. Chem. Phys. 59 (1973) 3229     *
 *    Carney G.D., Sprandel L.L. and Kern C.W. Adv. Chem. Phys. 37 (1977) 305*
 *****************************************************************************/
void Chedin5Force(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte (out)             */
 double         *pot,                           /* Energie (out)             */
                par[];                          /* Potentialparameter (in)   */
						/* par[0..5] = 11111, 11122, */
						/* 	       11133, 12222, */
						/* 	       12233, 13333  */
{
  vec_mt          force1s, force3s,	    /* Stretch forces 		     */
  		  force1t, force3t,         /* Angle (torque) forces         */
                  vr21,			    /* Vector r_21 = r(O1)-r(C)      */
                  vr23, 		    /* Vector r_23 = r(O2)-r(C)      */
                  n21,			    /* Unit vector in direction r_21 */
                  n23,			    /* Unit vector in direction r_23 */
                  ncross; 		    /* Cross product of n_21 and n_23*/
  double          r21,r23, nabs,		    /* Distances between atoms       */
                  alpha, sinr23, cosr23,/* Angle between n_21 and n_23   */	
                  xi1,xi2,xi3,		    /* Intrisic SPF coordinates      */
                  dudrq1,dudrq2,dudrq3,	    /* Derivatives of potential on ic*/
                  rq1,rq2,rq3,		    /* Quasi-normal SPF coordinates  */
                  dwdrq1,dwdrq2,dwdrq3,	    /* Derivatives of potential on qn*/
                  re, a1, a2, a3;	    /* Equilibrium C-O distance and  */
                  			    /* parameters of transition to   */ 
                  	 		    /* quasi-normal coordinates      */ 
 re = 1.161; a1 = 20.294; a2 = 7.665; a3 = 14.248;

  vr21[0] = r[1][0] - r[0][0];				
  vr21[1] = r[1][1] - r[0][1];
  vr21[2] = r[1][2] - r[0][2];

  r21 = ABS3(vr21);

  vr23[0] = r[2][0] - r[0][0];				
  vr23[1] = r[2][1] - r[0][1];
  vr23[2] = r[2][2] - r[0][2];

  r23 = ABS3(vr23);

  KreuzProd(vr21,vr23,ncross);
  nabs = ABS3(ncross);
  alpha = asin(nabs/(r21*r23));

  xi1 = 1 - re/r21;
  xi2 = alpha;
  xi3 = 1 - re/r23;

  rq1 = a1*(xi1+xi3);
  rq2 = a2*xi2;
  rq3 = a3*(xi1-xi3);

  *pot += par[0]*rq1*rq1*rq1*rq1*rq1 + par[1]*rq1*rq1*rq1*rq2*rq2 
        + par[2]*rq1*rq1*rq1*rq3*rq3 + par[3]*rq1*rq2*rq2*rq2*rq2 
        + par[4]*rq1*rq2*rq2*rq3*rq3 + par[5]*rq1*rq3*rq3*rq3*rq3;

  dudrq1 = 5*par[0]*rq1*rq1*rq1*rq1 + 3*par[1]*rq1*rq1*rq2*rq2 + par[3]*rq2*rq2*rq2*rq2 
  	 + (3*par[2]*rq1*rq1 + par[4]*rq2*rq2)*rq3*rq3 + par[5]*rq3*rq3*rq3*rq3;
  dudrq2 = 2*rq1*rq2*(par[1]*rq1*rq1 + 2*par[3]*rq2*rq2 + par[4]*rq3*rq3);
  dudrq3 = 2*rq1*rq3*(par[2]*rq1*rq1 + par[4]*rq2*rq2 + 2*par[5]*rq3*rq3);
					
					     /* Forces (0 is number for C) */

  f[1][0] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[0]
  	    -dudrq2*a2/(nabs/r23)*(vr23[0]/r23-vr21[0]*cos(PI-alpha)/r21);
  f[1][1] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[1]
  	    -dudrq2*a2/(nabs/r23)*(vr23[1]/r23-vr21[1]*cos(PI-alpha)/r21);
  f[1][2] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[2]
  	    -dudrq2*a2/(nabs/r23)*(vr23[2]/r23-vr21[2]*cos(PI-alpha)/r21);
  
  f[2][0] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[0]
  	    -dudrq2*a2/(nabs/r21)*(vr21[0]/r21-vr23[0]*cos(PI-alpha)/r23);
  f[2][1] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[1]
  	    -dudrq2*a2/(nabs/r21)*(vr21[1]/r21-vr23[1]*cos(PI-alpha)/r23);
  f[2][2] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[2]
  	    -dudrq2*a2/(nabs/r21)*(vr21[2]/r21-vr23[2]*cos(PI-alpha)/r23);
    
  f[0][0] = - f[1][0] - f[2][0];
  f[0][1] = - f[1][1] - f[2][1];
  f[0][2] = - f[1][2] - f[2][2];
}

/*****************************************************************************
 * Chedin6Force (6th order terms in potential, 10 parameters)		     *
 *  Potential: see Chedin A. and Teffo J.-L. J. Mol. Spectr. 107 (1984) 333  *
 *    Chedin A. and Teffo J.-L. J. Mol. Spectr. 107 (1984) 333  	     *
 *    Simons G., Parr R.G. and Finlan J.M. J. Chem. Phys. 59 (1973) 3229     *
 *    Carney G.D., Sprandel L.L. and Kern C.W. Adv. Chem. Phys. 37 (1977) 305*
 *****************************************************************************/
void Chedin6Force(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte (out)             */
 double         *pot,                           /* Energie (out)             */
                par[];                          /* Potentialparameter (in)   */
						/* par[0..9]= 111111, 111122,*/
						/*    111133, 112222, 112233,*/
						/*    113333, 222222, 222233,*/
						/*   	      223333, 333333 */
{
  vec_mt          force1s, force3s,	    /* Stretch forces 		     */
  		  force1t, force3t,         /* Angle (torque) forces         */
                  vr21,			    /* Vector r_21 = r(O1)-r(C)      */
                  vr23, 		    /* Vector r_23 = r(O2)-r(C)      */
                  n21,			    /* Unit vector in direction r_21 */
                  n23,			    /* Unit vector in direction r_23 */
                  ncross; 		    /* Cross product of n_21 and n_23*/
  double          r21,r23, nabs,		    /* Distances between atoms       */
                  alpha, sinr23, cosr23,/* Angle between n_21 and n_23   */	
                  xi1,xi2,xi3,		    /* Intrisic SPF coordinates      */
                  dudrq1,dudrq2,dudrq3,	    /* Derivatives of potential on ic*/
                  rq1,rq2,rq3,		    /* Quasi-normal SPF coordinates  */
                  dwdrq1,dwdrq2,dwdrq3,	    /* Derivatives of potential on qn*/
                  re, a1, a2, a3;	    /* Equilibrium C-O distance and  */
                  			    /* parameters of transition to   */ 
                  	 		    /* quasi-normal coordinates      */ 
 re = 1.161; a1 = 20.294; a2 = 7.665; a3 = 14.248;

  vr21[0] = r[1][0] - r[0][0];				
  vr21[1] = r[1][1] - r[0][1];
  vr21[2] = r[1][2] - r[0][2];

  r21 = ABS3(vr21);

  vr23[0] = r[2][0] - r[0][0];				
  vr23[1] = r[2][1] - r[0][1];
  vr23[2] = r[2][2] - r[0][2];

  r23 = ABS3(vr23);

  KreuzProd(vr21,vr23,ncross);
  nabs = ABS3(ncross);
  alpha = asin(nabs/(r21*r23));

  xi1 = 1 - re/r21;
  xi2 = alpha;
  xi3 = 1 - re/r23;

  rq1 = a1*(xi1+xi3);
  rq2 = a2*xi2;
  rq3 = a3*(xi1-xi3);

  *pot += par[0]*rq1*rq1*rq1*rq1*rq1*rq1 + par[1]*rq1*rq1*rq1*rq1*rq2*rq2 
        + par[2]*rq1*rq1*rq1*rq1*rq3*rq3 + par[3]*rq1*rq1*rq2*rq2*rq2*rq2 
        + par[4]*rq1*rq1*rq2*rq2*rq3*rq3 + par[5]*rq1*rq1*rq3*rq3*rq3*rq3
        + par[6]*rq2*rq2*rq2*rq2*rq2*rq2 + par[7]*rq2*rq2*rq2*rq2*rq3*rq3
        + par[8]*rq2*rq2*rq3*rq3*rq3*rq3 + par[9]*rq3*rq3*rq3*rq3*rq3*rq3;

  dudrq1 = 2*rq1*(3*par[0]*rq1*rq1*rq1*rq1 + 2*par[1]*rq1*rq1*rq2*rq2 
  	 + par[3]*rq2*rq2*rq2*rq2 + (2*par[2]*rq1*rq1 + par[4]*rq2*rq2)*rq3*rq3
  	 + par[5]*rq3*rq3*rq3*rq3);
  dudrq2 = 2*rq2*(par[1]*rq1*rq1*rq1*rq1 + 2*par[3]*rq1*rq1*rq2*rq2 
  	 + 3*par[6]*rq2*rq2*rq2*rq2 + (par[4]*rq1*rq1 + 2*par[7]*rq2*rq2)*rq3*rq3
  	 + par[8]*rq3*rq3*rq3*rq3);
  dudrq3 = 2*rq3*(par[2]*rq1*rq1*rq1*rq1 + par[4]*rq1*rq1*rq2*rq2 
  	 + par[6]*rq2*rq2*rq2*rq2 + 2*(par[5]*rq1*rq1 + par[8]*rq2*rq2)*rq3*rq3
  	 + 3*par[9]*rq3*rq3*rq3*rq3);
					
					     /* Forces (0 is number for C) */

  f[1][0] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[0]
  	    -dudrq2*a2/(nabs/r23)*(vr23[0]/r23-vr21[0]*cos(PI-alpha)/r21);
  f[1][1] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[1]
  	    -dudrq2*a2/(nabs/r23)*(vr23[1]/r23-vr21[1]*cos(PI-alpha)/r21);
  f[1][2] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*vr21[2]
  	    -dudrq2*a2/(nabs/r23)*(vr23[2]/r23-vr21[2]*cos(PI-alpha)/r21);
  
  f[2][0] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[0]
  	    -dudrq2*a2/(nabs/r21)*(vr21[0]/r21-vr23[0]*cos(PI-alpha)/r23);
  f[2][1] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[1]
  	    -dudrq2*a2/(nabs/r21)*(vr21[1]/r21-vr23[1]*cos(PI-alpha)/r23);
  f[2][2] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*vr23[2]
  	    -dudrq2*a2/(nabs/r21)*(vr21[2]/r21-vr23[2]*cos(PI-alpha)/r23);
    
  f[0][0] = - f[1][0] - f[2][0];
  f[0][1] = - f[1][1] - f[2][1];
  f[0][2] = - f[1][2] - f[2][2];
}

/*****************************************************************************
 * Morse-Force                                                               *
 *  Potential v(r) = D{ Exp[-2beta(R-Req)]-2 Exp[-beta(R-Req)]}              *
 ****************************************************************************/
void MorseForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = Req         (in) */
                                                /* par[1] = D     (in)       */
                                                /* par[2] = beta     (in)    */
{
  vec_mt		r_diff;
  double		r_abs, fak1, fak2;

  r_diff[0] = r[1][0] - r[0][0];
  r_diff[1] = r[1][1] - r[0][1];
  r_diff[2] = r[1][2] - r[0][2];

  r_abs = ABS3(r_diff);

  fak1 = par[1] * par[2] * 2;
  fak2 =exp(-2*par[2]*(r_abs-par[0]))-exp(-par[2]*(r_abs-par[0]));

  
  f[0][0] = fak1 * fak2 * r_diff[0]/r_abs;
  f[0][1] = fak1 * fak2 * r_diff[1]/r_abs;
  f[0][2] = fak1 * fak2 * r_diff[2]/r_abs;
  f[1][0] = -r_diff[0] * fak1;
  f[1][1] = -r_diff[1] * fak1;
  f[1][2] = -r_diff[2] * fak1;
  
  (*pot) 
+=par[1]*(exp(-2*par[2]*(r_abs-par[0]))-2*exp(-par[2]*(r_abs-par[0])));
}


/*****************************************************************************
 * BenzeneForce                                                              *
 *                                                                           *
 *  Die Matrix lokal2S gibt an, welche der lokalen Koordinaten an der inter- *
 *  nen Koordinate beteiligt ist. J. Chem. Phys. 87(5), 2564 (1987)          *
 *                                                                           *
 ****************************************************************************/
void BenzeneForce(r, f, par, pot)
 vec_mp         r,                   /* Teilchenpositionen            (in)   */
                f;                   /* Kraefte 	              (out)  */
 double         *pot,                /* Energie 	              (out)  */
                par[];               /* Potentialparameter            (in)   */
                                     /* par[0] = r_eq(C-C)            (in)   */
                                     /* par[1] = r_eq(C-H)            (in)   */
                                     /* par[2] = theta_eq(C-C-H)      (in)   */
                                     /* par[3] = theta_eq(C-C-C)      (in)   */
                                     /* par[4] = psi_eq(C-H-C-C)      (in)   */
                                     /* par[5] = phi_eq(H-C-C-H)      (in)   */
                                             
{
 double lokal2S[30][42] = {
     {0.40825,0.40825,0.40825,0.40825,0.40825,0.40825,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0.40825,0.40825,0.40825,0.40825,0.40825,0.40825,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0.20412,-0.20412,0.20412,-0.20412,0.20412,-0.20412,0.20412,-0.20412,0.20412,-0.20412,0.20412,-0.20412,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.40825,0.40825,-0.40825,0.40825,-0.40825,0.40825,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,-0.40825,0.40825,-0.40825,0.40825,-0.40825,0.40825,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {-0.40825,0.40825,-0.40825,0.40825,-0.40825,0.40825,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,-0.20412,0.20412,0.20412,-0.20412,-0.20412,0.20412,0.20412,-0.20412,-0.20412,0.20412,0.20412,-0.20412,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.57735,0.28868,0.28868,-0.57735,0.28868,0.28868,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.5,0.5,0,-0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,-0.57735,0.28868,0.28868,-0.57735,0.28868,0.28868,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,-0.5,0.5,0,-0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {-0.28868,0.57735,-0.28868,-0.28868,0.57735,-0.28868,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {-0.5,0,0.5,-0.5,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.25,0.25,0.25,-0.25,0,0,-0.25,0.25,0.25,-0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0.28868,-0.28868,-0.14434,0.14434,-0.14434,0.14434,0.28868,-0.28868,-0.14434,0.14434,-0.14434,0.14434,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.25,-0.25,0.25,-0.25,0,0,-0.25,0.25,-0.25,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0.28868,-0.28868,0.14434,-0.14434,-0.14434,0.14434,-0.28868,0.28868,-0.14434,0.14434,0.14434,-0.14434,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {-0.5,0,0.5,0.5,0,-0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0.28868,0.57735,0.28868,-0.28868,-0.57735,-0.28868,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,-0.57735,-0.28868,0.28868,0.57735,0.28868,-0.28868,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0.5,0.5,0,-0.5,-0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.40825,0.40825,0.40825,0.40825,0.40825,0.40825,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.40825,-0.40825,0.40825,-0.40825,0.40825,-0.40825},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.40825,-0.40825,0.40825,-0.40825,0.40825,-0.40825,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0.5,0,-0.5,-0.5,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.57735,-0.28868,0.28868,0.57735,0.28868,-0.28868,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.28868,0.57735,-0.28868,-0.28868,0.57735,-0.28868},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.5,0,0.5,-0.5,0,0.5},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.5,0.5,0,-0.5,0.5,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.57735,0.28868,0.28868,-0.57735,0.28868,0.28868,0,0,0,0,0,0}
 } ;

  /* */
/* Kraftkonstanten-Matrix  */
 double force[30][30] =  { 
     {7.616,0.168,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0.168,5.554,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0.878,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0.658,-0.237,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,-0.237,5.570,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,3.956,0.300,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0.300,0.824,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0.644,0,-0.136,0,0.308,0,-0.140,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0.644,0,-0.136,0,0.308,0,-0.140,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,-0.136,0,5.511,0,0.055,0,-0.065,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,-0.136,0,5.511,0,0.055,0,-0.065,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0.308,0,0.055,0,6.690,0,-0.398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0.308,0,0.055,0,6.690,0,-0.398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,-0.140,0,-0.065,0,-0.398,0,0.895,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,-0.140,0,-0.065,0,-0.398,0,0.895,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.926,0,0.209,0,0.149,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.926,0,0.209,0,0.149,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.209,0,7.380,0,0.582,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.209,0,7.380,0,0.582,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.149,0,0.582,0,5.568,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.149,0,0.582,0,5.568,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.249,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.202,0.249,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.249,0.519,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.337,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.337,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.160,0,-0.168,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.160,0,-0.168},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.168,0,0.420,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.168,0,0.420} 
 };
  

/* Nomenklatur fuer die lokalen Koordinaten ist der Literaturstelle entnommen.
               7
          12   |   8
            \  1  /
             \/ \/
             6   2            Die Atome 1-6 sind Kohlenstoffatome,
             |   |             
             5   3            7-12 sind Wasserstoffatome.
            / \ / \
          11   4   9
               |
              10 
*/
 int lokal_bind[12][2] =
     { {2,1},{3,2},{4,3},{5,4},{6,5},{1,6},{7,1},{8,2},{9,3},{10,4},{11,5},{12,6} };
 int lokal_theta[18][3] =
     { {6,1,7},{7,1,2},{1,2,8},{8,2,3},{2,3,9},{9,3,4},{3,4,10},{10,4,5},{4,5,11},{11,5,6},{5,6,12},
       {12,6,1},{6,1,2},{1,2,3},{2,3,4},{3,4,5},{4,5,6},{5,6,1}
     };
 int lokal_psi[6][4] =
     { {1,7,2,6},{2,8,3,1},{3,9,4,2},{4,10,5,3},{5,11,6,4},{6,12,1,5} };
 int lokal_phi[6][4] =
//     { {6,1,2,3},{1,2,3,4},{2,3,4,5},{3,4,5,6},{4,5,6,1},{5,6,1,2} };
    { {7,1,2,8},{8,2,3,9},{9,3,4,10},{10,4,5,11},{11,5,6,12},{12,6,1,7} };

 int inum,                   /* Nummer des Atoms                  */
     i, j, k, l,             /* Laufzahlen                        */
     a1,a2,a3,a4;            /* Atomindizes                       */
 
 vec_mt	dS_dr[30][12], delta_R[12][12];
 vec_mt d0,d1,d2,d3,n01,n12,n012;
 vec_mt f1,f2,f3,f4,v1,v2;
 double  mdyne2prg = 60221.7364;
 double  S[30], g;
 double  cos_phi, sin_phi, tan_psi, cos_psi, cos_psi_sin_phi, tan_psi_sin_sqr;
 double  r0,r1,r2,r3,r01,r12;
 double  theta,phi,psi;
 double  c1,c2,c3,c4,c5,c6,dv,scale;        
 double  c00,c01,c02,c11,c12,c22,d01,d12,c5d0,c5d1;

// Abstandsmatrix
  for (inum=0; inum<12; inum++)
 {
     if (inum<5)
     {
	 da(r+inum+1, r+inum, (vec_mp) d0);
	 for (k=0; k<3;k++)
	 {
	     delta_R[inum][inum+1][k] = d0[k];
	     delta_R[inum+1][inum][k] = -d0[k];
	 }
     }
     if (inum==5)
     {
	 da(r, r+5, (vec_mp) d0);
	 for (k=0; k<3;k++)
	 {
	     delta_R[5][0][k] = d0[k]; 
	     delta_R[0][5][k] = -d0[k];
	 }
     }
     if (inum>5)
     {
	 da(r+inum-6, r+inum, (vec_mp) d0);
	 for (k=0; k<3;k++)
	 {
	     delta_R[inum][inum-6][k] = d0[k];
	     delta_R[inum-6][inum][k] = -d0[k];
	 }
     }
 }


  /************************
  * Schleife ueber alle  *
  * internen Koordinaten *
  ************************/
 for (i=0; i < 30; i++)
 {  
     S[i] = 0;
    
     for (inum=0; inum < 12; inum++)
     {
	 dS_dr[i][inum][0] = 0;
	 dS_dr[i][inum][1] = 0;
	 dS_dr[i][inum][2] = 0;
     }
     
     //**********************
     // lokale Bindungen    *
     //**********************
     for (j=0; j < 12; j++)  
     {
	 if (lokal2S[i][j] != 0)
	 {
	     a1 = lokal_bind[j][0]-1;
	     a2 = lokal_bind[j][1]-1;

	     r0 = ABS3(delta_R[a1][a2]);
	     c1 = lokal2S[i][j] / r0;

	     if (j<6)
	     {
		 S[i] += ( r0 - par[0] ) * lokal2S[i][j];
	     }
	     else
	     {
		 S[i] += ( r0 - par[1] ) * lokal2S[i][j];
	     }   


	     dS_dr[i][a1][0] -= delta_R[a1][a2][0] * c1;
	     dS_dr[i][a1][1] -= delta_R[a1][a2][1] * c1;
	     dS_dr[i][a1][2] -= delta_R[a1][a2][2] * c1;

	     dS_dr[i][a2][0] += delta_R[a1][a2][0] * c1;
	     dS_dr[i][a2][1] += delta_R[a1][a2][1] * c1;
	     dS_dr[i][a2][2] += delta_R[a1][a2][2] * c1;
	 }
     }

     //***************************
     // lokale Winkel            *
     //***************************
     for (j=12; j < 30; j++) 
     {
	 if (lokal2S[i][j] !=0)
	 {
	     a1 = lokal_theta[j-12][0]-1;
	     a2 = lokal_theta[j-12][1]-1;
	     a3 = lokal_theta[j-12][2]-1;
	     
	     for (k=0; k<3;k++)
	     {
		 d0[k] = delta_R[a1][a2][k];
		 d1[k] = delta_R[a2][a3][k];
	     }

	     theta = acos(SCALPROD(d1,d0)/(ABS3(d1)*ABS3(d0)));
	     dv = -1/sin(theta);
	     c1 = dv / sqrt((C(d1,d1) * C(d0,d0)));
	     c2 = C(d1,d0)/C(d1,d1);
	     c3 = C(d1,d0)/C(d0,d0);

	     if (j<25)
	     {
		 scale = lokal2S[i][j] * par[1];
		 S[i] += ( PI - theta - par[2] ) * scale;
		 c1 *= scale;
	     }
	     else
	     {
		 scale = lokal2S[i][j] * par[0];
		 S[i] += ( PI - theta - par[3] ) * scale;
		 c1 *= scale;
	     }

	     dS_dr[i][a1][0] -= c1 * ( c3 * d0[0] - d1[0] );
	     dS_dr[i][a1][1] -= c1 * ( c3 * d0[1] - d1[1] );
	     dS_dr[i][a1][2] -= c1 * ( c3 * d0[2] - d1[2] );
		 
	     dS_dr[i][a2][0] -= c1 * ( (c2+1) * d1[0] - (c3+1) * d0[0] );
	     dS_dr[i][a2][1] -= c1 * ( (c2+1) * d1[1] - (c3+1) * d0[1] );
	     dS_dr[i][a2][2] -= c1 * ( (c2+1) * d1[2] - (c3+1) * d0[2] );
		  
	     dS_dr[i][a3][0] += c1 * ( c2 * d1[0] - d0[0] );
	     dS_dr[i][a3][1] += c1 * ( c2 * d1[1] - d0[1] );
	     dS_dr[i][a3][2] += c1 * ( c2 * d1[2] - d0[2] );
	 }
     }


    //***********************
    // improper torsions    *
    //***********************   
    for (j=30; j < 36; j++) 
     {
	 if (lokal2S[i][j] !=0)
	     {
		 a1 = lokal_psi[j-30][0]-1;
		 a2 = lokal_psi[j-30][1]-1;
		 a3 = lokal_psi[j-30][2]-1;
		 a4 = lokal_psi[j-30][3]-1;

		 for (k=0; k<3;k++) 
		 {
		     d0[k] = delta_R[a1][a2][k];
		     d1[k] = delta_R[a1][a3][k];
		     d2[k] = delta_R[a1][a4][k];
		 }

		 r0 = ABS3(d0);
		 r1 = ABS3(d1);
		 r2 = ABS3(d2);
		 
		 for ( k = 0; k < 3; k++)
		 {
		     d0[k] /= r0;
		     d1[k] /= r1;
		     d2[k] /= r2;
		 }

		 cos_phi = SCALPROD(d1, d2);
		 phi     = acos(cos_phi);
		 sin_phi = sin(phi);
//printf("%f\n",sin_phi);		 
		 KreuzProd(d2, d1, d3);
		 psi     = SCALPROD(d0, d3) / sin_phi;
		 psi     = asin(psi);
		 tan_psi = tan(psi);
		 cos_psi = cos(psi);
		 cos_psi_sin_phi = cos_psi * sin_phi;
		 tan_psi_sin_sqr = tan_psi / SQR(sin_phi);

		 f1[0] = ((d2[1]*d1[2]-d2[2]*d1[1])/cos_psi_sin_phi
			  - d0[0]*tan_psi) / r0;
		 f1[1] = ((d2[2]*d1[0]-d2[0]*d1[2])/cos_psi_sin_phi
			  - d0[1]*tan_psi) / r0;
		 f1[2] = ((d2[0]*d1[1]-d2[1]*d1[0])/cos_psi_sin_phi
			  - d0[2]*tan_psi) / r0;

		 f2[0] = ((d0[1]*d2[2]-d0[2]*d2[1])/cos_psi_sin_phi
			  -tan_psi_sin_sqr*(d1[0]-d2[0]*cos_phi))/r1;
		 f2[1] = ((d0[2]*d2[0]-d0[0]*d2[2])/cos_psi_sin_phi
			  -tan_psi_sin_sqr*(d1[1]-d2[1]*cos_phi))/r1;
		 f2[2] = ((d0[0]*d2[1]-d0[1]*d2[0])/cos_psi_sin_phi
			  -tan_psi_sin_sqr*(d1[2]-d2[2]*cos_phi))/r1;

		 f3[0] = ((d1[1]*d0[2]-d1[2]*d0[1])/cos_psi_sin_phi
			  -tan_psi_sin_sqr*(d2[0]-d1[0]*cos_phi))/r2;
		 f3[1] = ((d1[2]*d0[0]-d1[0]*d0[2])/cos_psi_sin_phi
			  -tan_psi_sin_sqr*(d2[1]-d1[1]*cos_phi))/r2;
		 f3[2] = ((d1[0]*d0[1]-d1[1]*d0[0])/cos_psi_sin_phi
			  -tan_psi_sin_sqr*(d2[2]-d1[2]*cos_phi))/r2;

		 c1 = lokal2S[i][j]  ;

		 S[i] += ( psi - par[4] ) * c1;

		 dS_dr[i][a1][0] -= (f1[0] + f2[0] + f3[0])*c1;
		 dS_dr[i][a1][1] -= (f1[1] + f2[1] + f3[1])*c1;
		 dS_dr[i][a1][2] -= (f1[2] + f2[2] + f3[2])*c1;
		  
		 dS_dr[i][a2][0] += f1[0]*c1;
		 dS_dr[i][a2][1] += f1[1]*c1;
		 dS_dr[i][a2][2] += f1[2]*c1;
		  
		 dS_dr[i][a3][0] += f2[0]*c1;
		 dS_dr[i][a3][1] += f2[1]*c1;
		 dS_dr[i][a3][2] += f2[2]*c1;
		  
		 dS_dr[i][a4][0] += f3[0]*c1;
		 dS_dr[i][a4][1] += f3[1]*c1;
		 dS_dr[i][a4][2] += f3[2]*c1;
	     } 
     }

    //************************
    //  Diederwinkel         *
    //************************   
    for (j=36; j < 42; j++) 
     {
	 if (lokal2S[i][j] !=0)
	 {
	     a1 = lokal_phi[j-36][0]-1;
	     a2 = lokal_phi[j-36][1]-1;
	     a3 = lokal_phi[j-36][2]-1;
	     a4 = lokal_phi[j-36][3]-1;

	     for (k=0;k<3;k++)
	     {
		 d0[k] = delta_R[a1][a2][k];
		 d1[k] = delta_R[a2][a3][k];
		 d2[k] = delta_R[a3][a4][k];
	     }

		 r0 = ABS3(d0);
		 r1 = ABS3(d1);
		 r2 = ABS3(d2);
		 
		 for ( k = 0; k < 3; k++)
		 {
		     d0[k] /= r0;
		     d1[k] /= r1;
		     d2[k] /= r2;
		 }

	     KreuzProd(d0, d1, n01 );
	     KreuzProd(d2, d1, n12 );
	     KreuzProd(n12, n01, n012 );

	     r01 = C(n01,n01);
	     r12 = C(n12,n12);

	     cos_phi = - C(n01,n12)/(sqrt(r01)*sqrt(r12));
	     sin_phi =   C(d1,n012)/(sqrt(r01)*sqrt(r12));
	     phi = atan2(sin_phi,cos_phi);

	     scale =lokal2S[i][j] ;

	     S[i] += ( phi - par[5] ) * scale;

	     for (k=0; k<3; k++)
	     {
	       dS_dr[i][a1][k] -= scale*n01[k]/(r0*r01);
	       dS_dr[i][a2][k] += scale*((1 - C(d0,d1)*r0/r1)*n01[k]/(r0*r01)
					+ C(d1,d2)*n12[k]/(r1*r12));
	       dS_dr[i][a3][k] += scale*((1 - C(d1,d2)*r2/r1)*n12[k]/(r2*r12)
					+ C(d0,d1)*n01[k]/(r1*r01));;
	       dS_dr[i][a4][k] -= scale*n12[k]/(r2*r12);
	     }


	 } 
     } 

  } /* Schleife ueber i */

 for (inum=0; inum < 12; inum++)
 {
     for (k=0; k<3; k++) 
	 {
	     f[inum][k] = 0;
	 }

     for (i=0; i < 30; i++)
     {
	 g=0;
	 for (j=0; j < 30; j++)
	 {
	     g +=  force[i][j] * mdyne2prg  * S[j];
	    
	     if (inum==0)
	     {
		  *pot += 0.5 * force[i][j] * mdyne2prg  * S[j] * S[i];
		  //	  printf("%i %i %10.5f\n",i,j,force[i][j] * mdyne2prg  * S[j] * S[i]);
	     }
	 }

	 f[inum][0] -= g * dS_dr[i][inum][0];
	 f[inum][1] -= g * dS_dr[i][inum][1];
	 f[inum][2] -= g * dS_dr[i][inum][2];
     }

   
 }
}

/*****************************************************************************
 * Water intramolecular forces (Oxygen atom should be first)       	     *
 * for different definitions of local coordinates:			     *
 *  0: quadratic in dr01, dr02, dtheta102 (K.Toukan,A.RahmanPRB31(1985)2643) *
 *  1: quadratic in dr01, dr02, dr12	       (	-//-		   ) *
 *  2: 1 with dr01^2, dr02^2 replaced by Morse (	-//-		   ) *
 *  3: quadratic wit qubic and quartic anharm in stretches		     *
					    (Sceats and Rice,JCP71(1979)973) *
 ****************************************************************************/
void WaterForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = pot. type        */
{
  vec_mt		r_diff, n01, n02, n12, dpot;
  double		r01, r02, r12, fconst[6], req, thetaeq, theta,
			r01_2, r02_2, r01_3, r02_3, 
			sintheta, costheta, fak1, fak2;
  int			iswitch = trunc(par[0]);

  r_diff[0] = r[1][0] - r[0][0];	/* Interatomic distances */
  r_diff[1] = r[1][1] - r[0][1];
  r_diff[2] = r[1][2] - r[0][2];
  r01 = ABS3(r_diff);
  fak1 = 1.0/r01;
  n01[0] = r_diff[0]*fak1;		/* Corresponding unit vectors     */
  n01[1] = r_diff[1]*fak1;
  n01[2] = r_diff[2]*fak1;

  r_diff[0] = r[2][0] - r[0][0];
  r_diff[1] = r[2][1] - r[0][1];
  r_diff[2] = r[2][2] - r[0][2];
  r02 = ABS3(r_diff);
  fak1 = 1.0/r02;
  n02[0] = r_diff[0]*fak1;
  n02[1] = r_diff[1]*fak1;
  n02[2] = r_diff[2]*fak1;

  switch (iswitch) {
   case 0:
			/* Force constants: alpha, beta, gamma, delta */
    fconst[0] = 5091.11*100; fconst[1] = 458.285*10; /* kJ/(mol A^2) */
    fconst[2] = 137.305*10; fconst[3] =-60.8236*10;
    req = 1.0; 		/* Equilibrium OH bond length */
    thetaeq = 1.91114; 	/* Equilibrium HOH angle (109.5 Degree) */

    costheta = SCALPROD(n01, n02); /* Valence angle */
    theta = acos(costheta);
    sintheta = 1.0/sin(theta);
    fak1 = sintheta/r01;
    fak2 = sintheta/r02;
    r01 -= req;		/* Deviations from eq. values */
    r02 -= req;
    theta -= thetaeq;
			/* Derivatives on local coords */
    dpot[0] = fconst[0]*r01 + fconst[2]*req*theta + fconst[3]*r02;
    dpot[1] = fconst[0]*r02 + fconst[2]*req*theta + fconst[3]*r01;
    dpot[2] = fconst[1]*SQR(req)*theta + fconst[2]*req*(r01 + r02);

    f[1][0] = - dpot[0]*n01[0] - dpot[2]*(n01[0]*costheta - n02[0])*fak1;
    f[1][1] = - dpot[0]*n01[1] - dpot[2]*(n01[1]*costheta - n02[1])*fak1;
    f[1][2] = - dpot[0]*n01[2] - dpot[2]*(n01[2]*costheta - n02[2])*fak1;
   
    f[2][0] = - dpot[1]*n02[0] - dpot[2]*(n02[0]*costheta - n01[0])*fak2;
    f[2][1] = - dpot[1]*n02[1] - dpot[2]*(n02[1]*costheta - n01[1])*fak2;
    f[2][2] = - dpot[1]*n02[2] - dpot[2]*(n02[2]*costheta - n01[2])*fak2;
  
    f[0][0] = - f[1][0] - f[2][0];
    f[0][1] = - f[1][1] - f[2][1];
    f[0][2] = - f[1][2] - f[2][2];
  
    (*pot) += 0.5*fconst[0]*(SQR(r01) + SQR(r02)) + 0.5*fconst[1]*SQR(req*theta)
	  + fconst[2]*req*theta*(r01 + r02) + fconst[3]*r01*r02;
   break;

   case 1:
			/* Force constants: a, b, c, d */
    fconst[0] = 561926.0; fconst[1] = 137485.0; /* kJ/(mol A^2) */
    fconst[2] =-88465.2; fconst[3] = 46731.8;
    req = 1.0; 		/* Equilibrium OH bond length */
    thetaeq = 1.91114; 	/* Equilibrium HOH angle (109.5 Degree) */

    r_diff[0] = r[2][0] - r[1][0];
    r_diff[1] = r[2][1] - r[1][1];
    r_diff[2] = r[2][2] - r[1][2];
    r12 = ABS3(r_diff);
    fak1 = 1.0/r12;
    n12[0] = r_diff[0]*fak1;
    n12[1] = r_diff[1]*fak1;
    n12[2] = r_diff[2]*fak1;

    r01 -= req;		/* Deviations from eq. values */
    r02 -= req;
    r12 -= req*sqrt(2*(1 - cos(thetaeq)));
			/* Derivatives on local coords */
    dpot[0] = fconst[0]*r01 + fconst[2]*r12 + fconst[3]*r02;
    dpot[1] = fconst[0]*r02 + fconst[2]*r12 + fconst[3]*r01;
    dpot[2] = fconst[1]*r12 + fconst[2]*(r01 + r02);

    f[1][0] = - dpot[0]*n01[0] + dpot[2]*n12[0];
    f[1][1] = - dpot[0]*n01[1] + dpot[2]*n12[1];
    f[1][2] = - dpot[0]*n01[2] + dpot[2]*n12[2];
   
    f[2][0] = - dpot[1]*n02[0] - dpot[2]*n12[0];
    f[2][1] = - dpot[1]*n02[1] - dpot[2]*n12[1];
    f[2][2] = - dpot[1]*n02[2] - dpot[2]*n12[2];
    
    f[0][0] = - f[1][0] - f[2][0];
    f[0][1] = - f[1][1] - f[2][1];
    f[0][2] = - f[1][2] - f[2][2];
  
    (*pot) += 0.5*fconst[0]*(SQR(r01) + SQR(r02)) + 0.5*fconst[1]*SQR(r12)
	  + fconst[2]*r12*(r01 + r02) + fconst[3]*r01*r02;
   break;

   case 2:
			/* Force constants: a, b, c, d, D, rho */
/* Toukan-Rahman for SPC/E
    fconst[0] = 561926.0; fconst[1] = 137485.0; // kJ/(mol A^2)
    fconst[2] =-88465.2; fconst[3] = 46731.8;
    fconst[4] = 42636.7; // dissociation constant kJ/mol 
    fconst[5] = 2.566;   // A^-1 
*/
// Lawrence-Skinner for TIP4P
    fconst[0] = 561908.0; fconst[1] = 103099.0; /* kJ/(mol A^2) */
    fconst[2] =-81901.1; fconst[3] = 52994.8;
    fconst[4] = 42636.7; /* dissociation constant kJ/mol */
    fconst[5] = 2.567;   /* A^-1 */

    req = 1.0; 		/* Equilibrium OH bond length */
    thetaeq = 1.91114*0.5; 	/* Half of eq. HOH angle (109.5 Degree)/2 */
 
    r_diff[0] = r[2][0] - r[1][0];
    r_diff[1] = r[2][1] - r[1][1];
    r_diff[2] = r[2][2] - r[1][2];
    r12 = ABS3(r_diff);
    fak1 = 1.0/r12;
    n12[0] = r_diff[0]*fak1;
    n12[1] = r_diff[1]*fak1;
    n12[2] = r_diff[2]*fak1;

    r01 -= req;		/* Deviations from eq. values */
    r02 -= req;
    r12 -= 2.0*req*sin(thetaeq);
 			/* Derivatives on local coords */
    fak1 = exp(-fconst[5]*r01);
    dpot[0] = 2.0*fconst[4]*fconst[5]*fak1*(1 - fak1)
	    + fconst[2]*r12 + fconst[3]*r02;
    fak2 = exp(-fconst[5]*r02);
    dpot[1] = 2.0*fconst[4]*fconst[5]*fak2*(1 - fak2)
	    + fconst[2]*r12 + fconst[3]*r01;
    dpot[2] = fconst[1]*r12 + fconst[2]*(r01 + r02);

    f[1][0] = - dpot[0]*n01[0] + dpot[2]*n12[0];
    f[1][1] = - dpot[0]*n01[1] + dpot[2]*n12[1];
    f[1][2] = - dpot[0]*n01[2] + dpot[2]*n12[2];
   
    f[2][0] = - dpot[1]*n02[0] - dpot[2]*n12[0];
    f[2][1] = - dpot[1]*n02[1] - dpot[2]*n12[1];
    f[2][2] = - dpot[1]*n02[2] - dpot[2]*n12[2];
    
    f[0][0] = - f[1][0] - f[2][0];
    f[0][1] = - f[1][1] - f[2][1];
    f[0][2] = - f[1][2] - f[2][2];
  
    fak1 = 1 - fak1;
    fak2 = 1 - fak2;
    (*pot) += fconst[4]*(SQR(fak1) + SQR(fak2)) + 0.5*fconst[1]*SQR(r12)
	   + fconst[2]*r12*(r01 + r02) + fconst[3]*r01*r02;
   break;

   case 3:
		/* Force constants: krr, kaa, kra, krrp, krrr, krrrr */
    fconst[1] = 48225.3; 
    fconst[2] = 28328.1; 
    fconst[4] =-551489.0; fconst[5] = 833403.0;

    fconst[0] = 508581.0*(1.0 - par[1])*(1.0 - par[1]); 
    fconst[3] = 7949.22*(-1.0 + 16.0530303030303*par[1]);

    req = 0.957200004282;//1.018; 		/* Equilibrium OH bond length */
    thetaeq = 1.82421810;//1.84481; 	/* Equilibrium HOH angle (105.7 Degree) */
    fconst[4] /= req; fconst[5] /= SQR(req);

    costheta = SCALPROD(n01, n02); /* Valence angle */
    theta = acos(costheta);
    sintheta = 1.0/sin(theta);
    fak1 = sintheta/r01;
    fak2 = sintheta/r02;
    r01 -= req;		/* Deviations from eq. values */
    r02 -= req;
    theta -= thetaeq;
    r01_2 = SQR(r01);  r02_2 = SQR(r02);
    r01_3 = r01*r01_2; r02_3 = r02*r02_2;
			/* Derivatives on local coords */
    dpot[0] = fconst[0]*r01 + fconst[2]*req*theta + fconst[3]*r02
	    + 3.0*fconst[4]*r01_2 + 4.0*fconst[5]*r01_3;
    dpot[1] = fconst[0]*r02 + fconst[2]*req*theta + fconst[3]*r01
	    + 3.0*fconst[4]*r02_2 + 4.0*fconst[5]*r02_3;
    dpot[2] = fconst[1]*SQR(req)*theta + fconst[2]*req*(r01 + r02);

    f[1][0] = - dpot[0]*n01[0] - dpot[2]*(n01[0]*costheta - n02[0])*fak1;
    f[1][1] = - dpot[0]*n01[1] - dpot[2]*(n01[1]*costheta - n02[1])*fak1;
    f[1][2] = - dpot[0]*n01[2] - dpot[2]*(n01[2]*costheta - n02[2])*fak1;
   
    f[2][0] = - dpot[1]*n02[0] - dpot[2]*(n02[0]*costheta - n01[0])*fak2;
    f[2][1] = - dpot[1]*n02[1] - dpot[2]*(n02[1]*costheta - n01[1])*fak2;
    f[2][2] = - dpot[1]*n02[2] - dpot[2]*(n02[2]*costheta - n01[2])*fak2;
  
    f[0][0] = - f[1][0] - f[2][0];
    f[0][1] = - f[1][1] - f[2][1];
    f[0][2] = - f[1][2] - f[2][2];
  
    (*pot) += 0.5*fconst[0]*(r01_2 + r02_2) + 0.5*fconst[1]*SQR(req*theta)
	  + fconst[2]*req*theta*(r01 + r02) + fconst[3]*r01*r02
	  + fconst[4]*(r01_3 + r02_3) + fconst[5]*(r01*r01_3 + r02*r02_3);

   break;

   case 4:
		/* Force constants: krr, kaa, kra, krrp, krrr, krrrr */
    fconst[1] = 48225.3; 
    fconst[2] = 28328.1; 
    fconst[4] =-551489.0; fconst[5] = 833403.0;

    fconst[0] = 508581.0*(1.0 - par[1])*(1.0 - par[1]); 
    fconst[3] = 7949.22*(-1.0 + 16.0530303030303*par[1]);

    req = 0.957200004282;//1.018; 		/* Equilibrium OH bond length */
    thetaeq = 1.82421810;//1.84481; 	/* Equilibrium HOH angle (105.7 Degree) */
    fconst[4] /= req; fconst[5] /= SQR(req);

    costheta = SCALPROD(n01, n02); /* Valence angle */
    theta = acos(costheta);
    sintheta = 1.0/sin(theta);
    fak1 = sintheta/r01;
    fak2 = sintheta/r02;
    r01 -= req*1.1;		/* Deviations from eq. values */
    r02 -= req;
    theta -= thetaeq;
    r01_2 = SQR(r01);  r02_2 = SQR(r02);
    r01_3 = r01*r01_2; r02_3 = r02*r02_2;
			/* Derivatives on local coords */
    dpot[0] = fconst[0]*r01 + fconst[2]*req*theta + fconst[3]*r02
	    + 3.0*fconst[4]*r01_2 + 4.0*fconst[5]*r01_3;
    dpot[1] = fconst[0]*r02 + fconst[2]*req*theta + fconst[3]*r01
	    + 3.0*fconst[4]*r02_2 + 4.0*fconst[5]*r02_3;
    dpot[2] = fconst[1]*SQR(req)*theta + fconst[2]*req*(r01 + r02);

    f[1][0] = - dpot[0]*n01[0] - dpot[2]*(n01[0]*costheta - n02[0])*fak1;
    f[1][1] = - dpot[0]*n01[1] - dpot[2]*(n01[1]*costheta - n02[1])*fak1;
    f[1][2] = - dpot[0]*n01[2] - dpot[2]*(n01[2]*costheta - n02[2])*fak1;
   
    f[2][0] = - dpot[1]*n02[0] - dpot[2]*(n02[0]*costheta - n01[0])*fak2;
    f[2][1] = - dpot[1]*n02[1] - dpot[2]*(n02[1]*costheta - n01[1])*fak2;
    f[2][2] = - dpot[1]*n02[2] - dpot[2]*(n02[2]*costheta - n01[2])*fak2;
  
    f[0][0] = - f[1][0] - f[2][0];
    f[0][1] = - f[1][1] - f[2][1];
    f[0][2] = - f[1][2] - f[2][2];
  
    (*pot) += 0.5*fconst[0]*(r01_2 + r02_2) + 0.5*fconst[1]*SQR(req*theta)
	  + fconst[2]*req*theta*(r01 + r02) + fconst[3]*r01*r02
	  + fconst[4]*(r01_3 + r02_3) + fconst[5]*(r01*r01_3 + r02*r02_3);

   break;

   default:
    f[1][0] = f[1][1] = f[1][2] = 0;
    f[2][0] = f[2][1] = f[2][2] = 0;
    f[0][0] = f[0][1] = f[0][2] = 0;
   break;
  }
}

/*****************************************************************************
 * Ammonia intramolecular forces (Nitrogen atom should be first)       	     *
 *  0: 	    (Hoy, Mills, Strey, Molecular Physics 24 (1972), 1265) *
 *  1: 	    (Martin. Lee, Taylor,JCP97(1992)8361) *
 ****************************************************************************/
void AmmoniaForce(r, f, par, pot)
 vec_mp         r,                              /* Teilchenpositionen (in)   */
                f;                              /* Kraefte 	      (out)  */
 double         *pot,                           /* Energie 	      (out)  */
                par[];                          /* Potentialparameter (in)   */
                                                /* par[0] = pot. type        */
{
  vec_mt		r_diff;
  double		rr[3], fconst[8], req, thetaeq, theta[3], fak[3], n[3][3],
			sintheta[3], costheta[3], dpot[6], V0, sinroh, sindiff;
  int			iswitch = trunc(par[0]), i;
 
   double f2[6][6], f3[6][6][6], f4[6][6][6][6], lo[3], o[3], s[6], temp, dsdr[6][4][3], g[4][3];
   int j,k,l;
 


for (i = 0; i < 3; i++)
 {
  r_diff[0] = r[i + 1][0] - r[0][0];	/* Interatomic distances */
  r_diff[1] = r[i + 1][1] - r[0][1];
  r_diff[2] = r[i + 1][2] - r[0][2];
  rr[i] = ABS3(r_diff);
  fak[i] = 1.0/rr[i];
  n[i][0] = r_diff[0]*fak[i];		/* Corresponding unit vectors     */
  n[i][1] = r_diff[1]*fak[i];
  n[i][2] = r_diff[2]*fak[i];
 } 


  switch (iswitch) {


   case 0:
			/* Force constants: alpha, beta, gamma, delta */
    fconst[0] = 424681; fconst[1] = 38059.9; /* kJ/(mol A^2) */
    fconst[2] = 903.321; fconst[3] = -3011.07; fconst[4] = 7527.67; fconst[5] = 18186.9; fconst[6] = -462861.0; fconst[7] = 666693.0;
    req = 1.025; 		/* Equilibrium NH bond length */
    thetaeq = 1.8675023; 	/* Equilibrium HNH angle (107 Degree) */


    costheta[0] = SCALPROD(n[1], n[2]); /* Valence angle */
    costheta[1] = SCALPROD(n[0], n[2]); /* Valence angle */
    costheta[2] = SCALPROD(n[0], n[1]); /* Valence angle */
    
    for (i = 0; i < 3; i++) 
     {
      theta[i] = acos(costheta[i]);
      sintheta[i] = 1.0/(sin(theta[i]));
      rr[i] -= req;
      theta[i] -= thetaeq;
     }

			/* Derivatives on local coords */
    dpot[0] = fconst[0]*rr[0] + fconst[2]*(rr[1] + rr[2]) + fconst[4]*theta[0] + fconst[5]*(theta[1] + theta[2]);
    dpot[1] = fconst[0]*rr[1] + fconst[2]*(rr[0] + rr[2]) + fconst[4]*theta[1] + fconst[5]*(theta[0] + theta[2]);
    dpot[2] = fconst[0]*rr[2] + fconst[2]*(rr[0] + rr[1]) + fconst[4]*theta[2] + fconst[5]*(theta[0] + theta[1]);
    dpot[3] = fconst[4]*rr[0] + fconst[5]*(rr[1] + rr[2]) + fconst[1]*theta[0] + fconst[3]*(theta[1] + theta[2]);
    dpot[4] = fconst[4]*rr[1] + fconst[5]*(rr[0] + rr[2]) + fconst[1]*theta[1] + fconst[3]*(theta[0] + theta[2]);
    dpot[5] = fconst[4]*rr[2] + fconst[5]*(rr[0] + rr[1]) + fconst[1]*theta[2] + fconst[3]*(theta[0] + theta[1]);

    f[1][0] = - dpot[0]*n[0][0] + dpot[4]*(n[2][0]-n[0][0]*costheta[1])*fak[0]*sintheta[1] + dpot[5]*(n[1][0]-n[0][0]*costheta[2])*fak[0]*sintheta[2] 
    		- 3.0*fconst[6]*rr[0]*rr[0]*n[0][0] - 4.0 * fconst[7]*rr[0]*rr[0]*rr[0]*n[0][0];
    f[1][1] = - dpot[0]*n[0][1] + dpot[4]*(n[2][1]-n[0][1]*costheta[1])*fak[0]*sintheta[1] + dpot[5]*(n[1][1]-n[0][1]*costheta[2])*fak[0]*sintheta[2]
    		- 3.0*fconst[6]*rr[0]*rr[0]*n[0][1] - 4.0 * fconst[7]*rr[0]*rr[0]*rr[0]*n[0][1];
    f[1][2] = - dpot[0]*n[0][2] + dpot[4]*(n[2][2]-n[0][2]*costheta[1])*fak[0]*sintheta[1] + dpot[5]*(n[1][2]-n[0][2]*costheta[2])*fak[0]*sintheta[2]
    		- 3.0*fconst[6]*rr[0]*rr[0]*n[0][2] - 4.0 * fconst[7]*rr[0]*rr[0]*rr[0]*n[0][2];

    f[2][0] = - dpot[1]*n[1][0] + dpot[3]*(n[2][0]-n[1][0]*costheta[0])*fak[1]*sintheta[0] + dpot[5]*(n[0][0]-n[1][0]*costheta[2])*fak[1]*sintheta[2]
    		- 3.0*fconst[6]*rr[1]*rr[1]*n[1][0] - 4.0 * fconst[7]*rr[1]*rr[1]*rr[1]*n[1][0];
    f[2][1] = - dpot[1]*n[1][1] + dpot[3]*(n[2][1]-n[1][1]*costheta[0])*fak[1]*sintheta[0] + dpot[5]*(n[0][1]-n[1][1]*costheta[2])*fak[1]*sintheta[2]
    		- 3.0*fconst[6]*rr[1]*rr[1]*n[1][1] - 4.0 * fconst[7]*rr[1]*rr[1]*rr[1]*n[1][1];
    f[2][2] = - dpot[1]*n[1][2] + dpot[3]*(n[2][2]-n[1][2]*costheta[0])*fak[1]*sintheta[0] + dpot[5]*(n[0][2]-n[1][2]*costheta[2])*fak[1]*sintheta[2]
    		- 3.0*fconst[6]*rr[1]*rr[1]*n[1][2] - 4.0 * fconst[7]*rr[1]*rr[1]*rr[1]*n[1][2];
   
    f[3][0] = - dpot[2]*n[2][0] + dpot[3]*(n[1][0]-n[2][0]*costheta[0])*fak[2]*sintheta[0] + dpot[4]*(n[0][0]-n[2][0]*costheta[1])*fak[2]*sintheta[1]
    		- 3.0*fconst[6]*rr[2]*rr[2]*n[2][0] - 4.0*fconst[7]*rr[2]*rr[2]*rr[2]*n[2][0];
    f[3][1] = - dpot[2]*n[2][1] + dpot[3]*(n[1][1]-n[2][1]*costheta[0])*fak[2]*sintheta[0] + dpot[4]*(n[0][1]-n[2][1]*costheta[1])*fak[2]*sintheta[1]
    		- 3.0*fconst[6]*rr[2]*rr[2]*n[2][1] - 4.0*fconst[7]*rr[2]*rr[2]*rr[2]*n[2][1];
    f[3][2] = - dpot[2]*n[2][2] + dpot[3]*(n[1][2]-n[2][2]*costheta[0])*fak[2]*sintheta[0] + dpot[4]*(n[0][2]-n[2][2]*costheta[1])*fak[2]*sintheta[1]
    		- 3.0*fconst[6]*rr[2]*rr[2]*n[2][2] - 4.0*fconst[7]*rr[2]*rr[2]*rr[2]*n[2][2];
  
    f[0][0] = - f[1][0] - f[2][0] - f[3][0];
    f[0][1] = - f[1][1] - f[2][1] - f[3][1];
    f[0][2] = - f[1][2] - f[2][2] - f[3][2];

    
  
    (*pot) += 0.5*fconst[0]*(rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2]) + 
              0.5*fconst[1]*(theta[0]*theta[0] +  theta[1]*theta[1] + theta[2]*theta[2]) + 
	      fconst[2]*(rr[0]*rr[1] + rr[0]*rr[2] + rr[1]*rr[2]) +
              fconst[3]*(theta[0]*theta[1] +  theta[0]*theta[2] + theta[1]*theta[2]) + 
              fconst[4]*(rr[0]*theta[0] +  rr[1]*theta[1] + rr[2]*theta[2]) + 
              fconst[5]*(rr[0]*(theta[1] + theta[2]) +  rr[1]*(theta[0] + theta[2]) + rr[2]*(theta[0] + theta[1])) + 
	      fconst[6]*(rr[0]*rr[0]*rr[0] + rr[1]*rr[1]*rr[1] + rr[2]*rr[2]*rr[2]) + 
	      fconst[7]*(rr[0]*rr[0]*rr[0]*rr[0] + rr[1]*rr[1]*rr[1]*rr[1] + rr[2]*rr[2]*rr[2]*rr[2]);

   break;

   case 1:
   
   
   
    f2[0][0] = f2[1][1] = f2[2][2] = 423578.0;
    f2[0][1] = f2[0][2] = f2[1][2] = -2225.66;
    f2[0][3] = f2[1][4] = f2[2][5] = 2466.06;  
    f2[0][4] = f2[0][5] = f2[1][3] = f2[1][5] = f2[2][3] = f2[2][4] = 14537.5; 
    f2[3][3] = f2[4][4] = f2[5][5] = 38792.5;
    f2[3][4] = f2[3][5] = f2[4][5] = -3042.08;

   for (i = 0; i<6; i++)
   for (j = 0; j<i; j++) f2[i][j] = f2[j][i];
   
    f3[0][0][0] = f3[1][1][1] = f3[2][2][2] = -2715804.19;
    f3[0][0][1] = f3[0][1][1] = f3[0][0][2] = f3[0][2][2] = f3[1][1][2] = f3[1][2][2] =-736.507;
    f3[0][0][3] = f3[1][1][4] = f3[2][2][5] =-564.876;
    f3[1][1][3] = f3[2][2][3] = f3[0][0][4] = f3[2][2][4] = f3[0][0][5] = f3[1][1][5] =820.817;
    f3[0][1][2] = -4061.93;
    f3[0][1][3] = f3[0][2][3] = f3[0][1][4] = f3[1][2][4] = f3[0][2][5] = f3[1][2][5] = 5763.18;
    f3[1][2][3] = f3[0][2][4] = f3[0][1][5] = -21191.3;
    f3[0][3][3] = f3[1][4][4] = f3[2][5][5] = 3445.87;
    f3[0][3][4] = f3[1][3][4] = f3[0][3][5] = f3[2][3][5] = f3[1][4][5] = f3[2][4][5] = 5863.15;
    f3[1][3][3] = f3[2][3][3] = f3[0][4][4] = f3[2][4][4] = f3[0][5][5] = f3[1][5][5] = -12538.7;
    f3[2][3][4] = f3[1][3][5] = f3[0][4][5] = 3008.66;
    f3[3][3][3] = f3[4][4][4] = f3[5][5][5] = -40866.8;
    f3[3][4][4] = f3[3][3][4] = f3[3][3][5] = f3[3][5][5] = f3[4][4][5] = f3[4][5][5] = -1091.81;
    f3[3][4][5] = 0.0;//155.227;
    
   for (i = 0; i<6; i++)
   for (j = 0; j<i; j++)
   for (k = 0; k<j; k++)
   	{
	f3[j][k][i] = f3[k][j][i];
	f3[k][i][j] = f3[k][j][i];
	f3[j][i][k] = f3[k][j][i];
	f3[i][k][j] = f3[k][j][i];
	f3[i][j][k] = f3[k][j][i];
   	}
   for (i = 0; i<6; i++)
   for (j = 0; j<i; j++)
   	{
	f3[j][i][j] = f3[j][j][i];
	f3[i][j][j] = f3[j][j][i]; 
	f3[i][i][j] = f3[j][i][i]; 
 	f3[i][j][i] = f3[j][i][i]; 
  	}


    f4[0][0][0][0] = f4[1][1][1][1] = f4[2][2][2][2] = 15554622.05;
    f4[0][0][0][1] = f4[0][0][0][2] = f4[0][1][1][1] = f4[1][1][1][2] = f4[0][2][2][2] = f4[1][2][2][2] = -51552.5;
    f4[0][0][0][3] = f4[1][1][1][4] = f4[2][2][2][5] = 12930.1;
    f4[0][0][0][4] = f4[0][0][0][5] = f4[1][1][1][3] = f4[1][1][1][5] = f4[2][2][2][3] = f4[2][2][2][4] = -68281.4;
    f4[0][0][1][1] = f4[0][0][2][2] = f4[1][1][2][2] = -9336.72;
    f4[0][0][1][2] = f4[0][1][1][2] = f4[0][1][2][2] = -7153.7;
    f4[0][0][1][3] = f4[0][0][2][3] = f4[0][1][1][4] = f4[1][1][2][4] = f4[0][2][2][5] = f4[1][2][2][5] = -1432.67;
    f4[0][0][1][4] = f4[0][0][2][5] = f4[0][1][1][3] = f4[1][1][2][5] = f4[0][2][2][3] = f4[1][2][2][4] = -2035.48;
    f4[0][0][2][4] = f4[0][0][1][5] = f4[1][1][2][3] = f4[0][1][1][5] = f4[1][2][2][3] = f4[0][2][2][4] = 072.8679;
    f4[0][0][3][3] = f4[1][1][4][4] = f4[2][2][5][5] = -1769.91;
    f4[0][0][3][4] = f4[0][0][3][5] = f4[1][1][3][4] = f4[1][1][4][5] = f4[2][2][3][5] = f4[2][2][4][5] = 10245.5;
    f4[0][0][4][4] = f4[0][0][5][5] = f4[1][1][3][3] = f4[1][1][5][5] = f4[2][2][3][3] = f4[2][2][4][4] = -15549.2;
    f4[0][0][4][5] = f4[1][1][3][5] = f4[2][2][3][4] = 21810.4;
    f4[0][1][2][3] = f4[0][1][2][4] = f4[0][1][2][5] = -5686.1;
    f4[0][1][3][3] = f4[0][1][4][4] = f4[0][2][3][3] = f4[0][2][5][5] = f4[1][2][4][4] = f4[1][2][5][5] = 2925.55;
    f4[0][1][3][5] = f4[0][1][4][5] = f4[0][2][3][4] = f4[0][2][4][5] = f4[1][2][3][4] = f4[1][2][3][5] = 3470.56;
    f4[0][1][3][4] = f4[0][2][3][5] = f4[1][2][4][5] = 2802.7;
    f4[0][3][3][3] = f4[1][4][4][4] = f4[2][5][5][5] = -3598.23;
    f4[0][3][4][5] = f4[1][3][4][5] = f4[2][3][4][5] = 2642.51;
    f4[0][3][3][4] = f4[0][3][3][5] = f4[1][3][4][4] = f4[1][4][4][5] = f4[2][3][5][5] = f4[2][4][5][5] = 3077.91;
    f4[1][2][3][3] = f4[0][2][4][4] = f4[0][1][5][5] = 28917.1;
    f4[1][3][3][4] = f4[2][3][3][5] = f4[0][3][4][4] = f4[2][4][4][5] = f4[0][3][5][5] = f4[1][4][5][5] = 512.484;
    f4[2][3][3][4] = f4[1][3][3][5] = f4[0][4][4][5] = f4[2][3][4][4] = f4[1][3][5][5] = f4[0][4][5][5] = 8252.74;
    f4[1][3][3][3] = f4[2][3][3][3] = f4[0][4][4][4] = f4[2][4][4][4] = f4[0][5][5][5] = f4[1][5][5][5] = 23941.6;
    f4[3][3][3][3] = f4[4][4][4][4] = f4[5][5][5][5] = -48653.4;
    f4[3][3][3][4] = f4[3][3][3][5] = f4[3][4][4][4] = f4[4][4][4][5] = f4[3][5][5][5] = f4[4][5][5][5] = -3544.63;
    f4[3][3][4][4] = f4[3][3][5][5] = f4[4][4][5][5] = 19350.9;
    f4[3][3][4][5] = f4[3][4][4][5] = f4[3][4][5][5] = -6528.0;


   for (i = 0; i<6; i++)
   for (j = 0; j<i; j++)
   for (k = 0; k<j; k++)
   for (l = 0; l<k; l++)
   	{
	f4[j][k][l][i] = f4[l][k][j][i]; 
	f4[k][l][i][j] = f4[l][k][j][i]; 
	f4[l][i][j][k] = f4[l][k][j][i]; 
	f4[j][i][k][l] = f4[l][k][j][i]; 
	f4[i][k][l][j] = f4[l][k][j][i]; 
	f4[k][l][j][i] = f4[l][k][j][i]; 
	f4[l][j][i][k] = f4[l][k][j][i]; 
	f4[k][j][i][l] = f4[l][k][j][i]; 
	f4[j][i][l][k] = f4[l][k][j][i]; 
	f4[i][l][k][j] = f4[l][k][j][i]; 
	f4[l][k][i][j] = f4[l][k][j][i]; 
	f4[l][j][k][i] = f4[l][k][j][i]; 
	f4[j][k][i][l] = f4[l][k][j][i]; 
	f4[k][i][l][j] = f4[l][k][j][i]; 
	f4[i][l][j][k] = f4[l][k][j][i]; 
	f4[i][k][j][l] = f4[l][k][j][i]; 
	f4[k][j][l][i] = f4[l][k][j][i]; 
	f4[j][l][i][k] = f4[l][k][j][i]; 
	f4[l][i][k][j] = f4[l][k][j][i]; 
	f4[i][j][k][l] = f4[l][k][j][i]; 
	f4[i][j][l][k] = f4[l][k][j][i]; 
	f4[k][i][j][l] = f4[l][k][j][i]; 
	f4[j][l][k][i] = f4[l][k][j][i]; 
   	} 
	
   for (i = 0; i<6; i++)
   for (j = 0; j<i; j++)
   for (k = 0; k<j; k++)
	{
	f4[k][k][i][j] = f4[k][k][j][i];
	f4[j][k][k][i] = f4[k][k][j][i];
	f4[i][k][k][j] = f4[k][k][j][i];
	f4[j][i][k][k] = f4[k][k][j][i];
	f4[i][j][k][k] = f4[k][k][j][i];
	f4[j][k][i][k] = f4[k][k][j][i];
	f4[i][k][j][k] = f4[k][k][j][i];
	f4[k][i][k][j] = f4[k][k][j][i];
	f4[k][j][k][i] = f4[k][k][j][i];
	f4[k][i][j][k] = f4[k][k][j][i];
	f4[k][j][i][k] = f4[k][k][j][i];
	
	f4[i][j][j][k] = f4[k][j][j][i];
	f4[j][j][k][i] = f4[k][j][j][i];
	f4[j][j][i][k] = f4[k][j][j][i];
	f4[i][k][j][j] = f4[k][j][j][i];
	f4[k][i][j][j] = f4[k][j][j][i];
	f4[i][j][k][j] = f4[k][j][j][i];
	f4[k][j][i][j] = f4[k][j][j][i];
	f4[j][i][k][j] = f4[k][j][j][i];
	f4[j][k][i][j] = f4[k][j][j][i];
	f4[j][i][j][k] = f4[k][j][j][i];
	f4[j][k][j][i] = f4[k][j][j][i];
	
	f4[j][k][i][i] = f4[k][j][i][i];
	f4[i][i][k][j] = f4[k][j][i][i];
	f4[i][i][j][k] = f4[k][j][i][i];
	f4[k][i][i][j] = f4[k][j][i][i];
	f4[j][i][i][k] = f4[k][j][i][i];
	f4[i][j][k][i] = f4[k][j][i][i];
	f4[i][k][j][i] = f4[k][j][i][i];
	f4[i][j][i][k] = f4[k][j][i][i];
	f4[i][k][i][j] = f4[k][j][i][i];
	f4[j][i][k][i] = f4[k][j][i][i];
	f4[k][i][j][i] = f4[k][j][i][i];
	}
	
   for (i = 0; i<6; i++)
   for (j = 0; j<i; j++)
	{
	f4[j][i][i][j] = f4[j][j][i][i];
	f4[i][i][j][j] = f4[j][j][i][i];
	f4[i][j][j][i] = f4[j][j][i][i];
	f4[j][i][j][i] = f4[j][j][i][i];
	f4[i][j][i][j] = f4[j][j][i][i];
	}
	
   for (i = 0; i<6; i++)
   for (j = 0; j<i; j++)
	{
	f4[j][j][i][j] = f4[j][j][j][i];
	f4[j][i][j][j] = f4[j][j][j][i];
	f4[i][j][j][j] = f4[j][j][j][i];
	
	f4[i][j][i][i] = f4[j][i][i][i];
	f4[i][i][j][i] = f4[j][i][i][i];
	f4[i][i][i][j] = f4[j][i][i][i];
	}
	
  
    
    req = 1.0141; 		/* Equilibrium NH bond length */
    thetaeq = 1.84377; 	/* Equilibrium HNH angle (105.64 Degree) */

    costheta[0] = SCALPROD(n[1], n[2]); /* Valence angle */
    costheta[1] = SCALPROD(n[0], n[2]); /* Valence angle */
    costheta[2] = SCALPROD(n[0], n[1]); /* Valence angle */

    for (i = 0; i < 3; i++) 
     {
      theta[i] = acos(costheta[i]);
      sintheta[i] = 1.0/(sin(theta[i]));
      rr[i] -= req;
      theta[i] -= thetaeq;
     }

   
   for (i=0;i<4;i++)
   for (j=0;j<3;j++) f[i][j] = 0.0;
     

			/* Derivatives on local coords */

   for (i=0; i<6; i++)
   for (k=0; k<4; k++)
   for (j=0; j<3; j++) dsdr[i][k][j] = 0.0;
   for (i=0; i<6; i++) dpot[i] = 0.0;
 
 


			/*Lokale Koordinaten*/
for (i=0; i<3; i++) 
{
s[i] = rr[i];
s[i+3] = theta[i];
}


			/* Derivatives on local coords */

for (i=0; i<6; i++)
for (k=0; k<6; k++) dpot[i] += f2[i][k]*s[k];    

temp = 0.0;
for (i=0; i<6; i++)
for (j=0; j<6; j++)
for (k=0; k<6; k++) dpot[i] += 0.5*f3[i][j][k]*s[j]*s[k];


for (i=0; i<6; i++) 
for (j=0; j<6; j++)
for (k=0; k<6; k++)
for (l=0; l<6; l++) dpot[i] += 0.166667*f4[i][j][k][l]*s[j]*s[k]*s[l];


	for (i=0; i<3; i++)
	 for (j=0; j<3; j++) f[i+1][j] = -dpot[i]*n[i][j];
	 
	 for (j=0; j<3; j++) 
	 {
	 dsdr[3][2][j] = -(n[2][j] - n[1][j]*costheta[0])*fak[1]*sintheta[0];
	 f[2][j] += -dpot[3]*dsdr[3][2][j];

	  dsdr[3][3][j] = -(n[1][j] - n[2][j]*costheta[0])*fak[2]*sintheta[0];
	  f[3][j] += -dpot[3]*dsdr[3][3][j];
	  
	dsdr[4][3][j] = -(n[0][j] - n[2][j]*costheta[1])*fak[2]*sintheta[1];
	f[3][j] += -dpot[4]*dsdr[4][3][j];

	dsdr[4][1][j] = -(n[2][j] - n[0][j]*costheta[1])*fak[0]*sintheta[1];
	f[1][j] += -dpot[4]*dsdr[4][1][j];
	
	dsdr[5][2][j] = -(n[0][j] - n[1][j]*costheta[2])*fak[1]*sintheta[2];
	f[2][j] += -dpot[5]*dsdr[5][2][j];
	
	dsdr[5][1][j] = -(n[1][j] - n[0][j]*costheta[2])*fak[0]*sintheta[2];
	f[1][j] += -dpot[5]*dsdr[5][1][j];
	}

    f[0][0] = - f[1][0] - f[2][0] - f[3][0];
    f[0][1] = - f[1][1] - f[2][1] - f[3][1];
    f[0][2] = - f[1][2] - f[2][2] - f[3][2];
 


/*Potentielle Energie*/

/*Zweiter Ordnung*/     
    temp = 0.0;
    for (j=0; j<6; j++) temp += f2[j][j]*s[j]*s[j]; 
    (*pot) += 0.5*temp;    
    for (i=0; i<6; i++)
    for (j=i+1; j<6; j++) (*pot) += f2[i][j]*s[i]*s[j]; 


/*Dritter Ordnung*/ 
    temp = 0.0;
    for (k=0; k<6; k++) temp += f3[k][k][k]*s[k]*s[k]*s[k];
    (*pot) += 0.1666667*temp;
    temp = 0.0;
    for (j=0; j<6; j++)
    for (k=j+1; k<6; k++) 
      {
       temp += f3[j][k][k]*s[j]*s[k]*s[k];
       temp += f3[j][j][k]*s[j]*s[j]*s[k];
      } 
    (*pot) += 0.5*temp;
    for (i=0; i<6; i++)
    for (j=i+1; j<6; j++)
    for (k=j+1; k<6; k++) (*pot) += f3[i][j][k]*s[i]*s[j]*s[k];

/*Vierter Ordnung*/ 
    temp = 0.0;
    for (k=0; k<6; k++) temp += f4[k][k][k][k]*s[k]*s[k]*s[k]*s[k];
    (*pot) += 0.041666667*temp;
    
    temp = 0.0;
    for (j=0; j<6; j++)
    for (k=j+1; k<6; k++) 
    {
    	temp += f4[j][k][k][k]*s[j]*s[k]*s[k]*s[k];
        temp += f4[k][j][j][j]*s[k]*s[j]*s[j]*s[j];
	}
    (*pot) += 0.166666667*temp;
    
    temp = 0.0;
    for (j=0; j<6; j++)
    for (k=j+1; k<6; k++) temp += f4[j][j][k][k]*s[j]*s[j]*s[k]*s[k];
    (*pot) += 0.25*temp;
    
    temp = 0.0;
    for (i=0; i<6; i++)
    for (j=i+1; j<6; j++) 
    for (k=j+1; k<6; k++)
    {
    	temp += f4[i][j][k][k]*s[i]*s[j]*s[k]*s[k];
        temp += f4[j][k][i][i]*s[j]*s[k]*s[i]*s[i];
        temp += f4[k][i][j][j]*s[k]*s[i]*s[j]*s[j];
	}
    (*pot) += 0.5*temp;
    
    for (i=0; i<6; i++)
    for (j=i+1; j<6; j++)
    for (k=j+1; k<6; k++)
    for (l=k+1; l<6; l++) (*pot) += f4[i][j][k][l]*s[i]*s[j]*s[k]*s[l];

 break;

   default:
    f[1][0] = f[1][1] = f[1][2] = 0;
    f[2][0] = f[2][1] = f[2][2] = 0;
    f[0][0] = f[0][1] = f[0][2] = 0;
   break;
  }
}

/*****************************************************************************
 * IntraPotForce  Berechnung der INTRAmolekularen Kraefte fuer die Molekuele *
 *                einer Species.                                             *
 *****************************************************************************/
void  IntraForceCalc(site, site_force, species, h, pot, stress)
 double         **site,                 /* Site co-ordinate arrays      (in) */
                **site_force;           /* Site force arrays            (out)*/
 spec_mp	species;
 double		h;			/* Boxlaenge			(in) */
 double		*pot;			/* Potentielle Energie 		(out)*/
 mat_mt         stress;                 /* Stress virial                (out)*/
{
  int		ispec, ipot, iintra, i, j,
                issite[NINVSIT];
  vec_mt	r[NINVSIT], f[NINVSIT];
  double	par[NPOTP], pottemp,
                h_2 = 0.5 * h;
  
  for (ispec=0; ispec < species->nmols; ispec++)
  { 
   species->intrpot[ispec] = -(*pot);	/* Extracting inter-molecular part */
    for (ipot = 0; ipot < species->npotls; ipot++)
    {
      switch( (species->Pots[ipot]).pot_id)
      {
      /******************
       * TorForce       *
       ******************/
      case 0:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  
	  for(j = 0; j < 4; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  TorForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	  }
	}
	break;

      /******************
       * AngleForce     *
       ******************/		  

      case 1:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  
	  for(j = 0; j < 3; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  AngleForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	  }
	}
	break;

      /******************
       * StretchForce    *
       ******************/
      case 2:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  
	  for(j = 0; j < 2; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
		  r[j][i] -= h;
	    }
	  }

	  StretchForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	  }
          stress[0][0]  += (r[0][0] - r[1][0])*f[0][0];
          stress[0][1]  += (r[0][0] - r[1][0])*f[0][1];
          stress[0][2]  += (r[0][0] - r[1][0])*f[0][2];
          stress[1][1]  += (r[0][1] - r[1][1])*f[0][1];
          stress[1][2]  += (r[0][1] - r[1][1])*f[0][2];
          stress[2][2]  += (r[0][2] - r[1][2])*f[0][2];
	}
	break;

      /******************
       * WagForce       *
       ******************/
      case 3:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];

	  /*
	  site[0][0] = 0.0; site[1][0] = 0.0; site[2][0] = 0.0;
	  site[0][1] = 0.0; site[1][1] = 2.0; site[2][1] = 0.0;
	  site[0][2] = -sqrt(0.75); site[1][2] = -0.5; site[2][2] = 0.0;
	  site[0][3] =  sqrt(0.75); site[1][3] = -0.5; site[2][3] = 0.0;
	   */
	  
	  for(j = 0; j < 4; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  WagForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	  }
	}
	break;	

      /******************
       * NoBondForce    *
       ******************/
      case 4:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  
	  for(j = 0; j < 2; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
		  r[j][i] -= h;
	    }
	  }

	  NoBondForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	  }
	}
	break;

      /******************
       * BellForce      *
       ******************/
       case 5:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  
	  for(j = 0; j < 6; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  BellForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	    site_force[i][issite[4]] += f[4][i];
	    site_force[i][issite[5]] += f[5][i];
	  }
	}
	break;	

      /******************
       * Chedin2Force   *
       ******************/
      case 6:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  
	  for(j = 0; j < 3; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  Chedin2Force( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	  }
	}
	break;

      /******************
       * Chedin34Force  *
       ******************/
      case 7:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];
	  par[6] = (species->Pots[ipot]).p[iintra][6];
	  par[7] = (species->Pots[ipot]).p[iintra][7];
	  par[8] = (species->Pots[ipot]).p[iintra][8];
	  
	  for(j = 0; j < 3; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  Chedin34Force( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	  }
	}
	break;

      /******************
       * Chedin5Force   *
       ******************/
      case 8:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];
	  
	  for(j = 0; j < 3; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  Chedin5Force( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	  }
	}
	break;

      /******************
       * Chedin6Force   *
       ******************/
      case 9:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];
	  par[6] = (species->Pots[ipot]).p[iintra][6];
	  par[7] = (species->Pots[ipot]).p[iintra][7];
	  par[8] = (species->Pots[ipot]).p[iintra][8];
	  par[9] = (species->Pots[ipot]).p[iintra][9];
	  
	  for(j = 0; j < 3; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  Chedin6Force( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	  }
	}
	break;

      /******************
       * XY2Z2SymStr    *
       ******************/
      case 10:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];
	  par[6] = (species->Pots[ipot]).p[iintra][6];
	  
	  for(j = 0; j < 5; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  XY2Z2SymStrForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	    site_force[i][issite[4]] += f[4][i];
	  }
	}
	break;

      /******************
       * XY2Z2SymBend   *
       ******************/
      case 11:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];
	  par[6] = (species->Pots[ipot]).p[iintra][6];
	  par[7] = (species->Pots[ipot]).p[iintra][7];
	  par[8] = (species->Pots[ipot]).p[iintra][8];
	  
	  for(j = 0; j < 5; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  XY2Z2SymBendForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	    site_force[i][issite[4]] += f[4][i];
	  }
	}
	break;

      /******************
       * XY2Z2          *
       ******************/
      case 12:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
/*	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];
	  par[6] = (species->Pots[ipot]).p[iintra][6];
	  par[7] = (species->Pots[ipot]).p[iintra][7];
	  par[8] = (species->Pots[ipot]).p[iintra][8];
	  par[9] = (species->Pots[ipot]).p[iintra][9];
	  par[10]= (species->Pots[ipot]).p[iintra][10];*/
	  
	  for(j = 0; j < 5; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  XY2Z2Force( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	    site_force[i][issite[4]] += f[4][i];
	  }
	}
	break;

      /******************
       * XY4            *
       ******************/
      case 13:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];
	  par[6] = (species->Pots[ipot]).p[iintra][6];
	  
	  for(j = 0; j < 5; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  XY4Force( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	    site_force[i][issite[4]] += f[4][i];
	  }
	}
	break;

      /******************
       * Morse          *
       ******************/
      case 14:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  
	  for(j = 0; j < 2; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
			r[j][i] -= h;
	    }
	  }
	  
	  MorseForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	  }
	}
	break;

      /******************
       * Benzene        *
       ******************/
      case 15:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  par[3] = (species->Pots[ipot]).p[iintra][3];
	  par[4] = (species->Pots[ipot]).p[iintra][4];
	  par[5] = (species->Pots[ipot]).p[iintra][5];

	  for(j = 0; j < 12; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;

	    for ( i = 0; i < 3; i++)
	    {
		r[j][i] = site[i][issite[j]];
		if (j > 0)
		    if (r[0][i] > 0.0) 
		    {
			if ( r[0][i] - r[j][i] > h_2 )
			    r[j][i] += h;
		    }
		    else
			if ( r[j][i] - r[0][i] > h_2 )
			    r[j][i] -= h;
	    }

	  } 

	  BenzeneForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
    	   for (j = 0; j < 12; j++) site_force[i][issite[j]] += f[j][i];

	}
	break;

      /*********************
       * StretchCrossForce *
       *********************/
      case 16:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  par[2] = (species->Pots[ipot]).p[iintra][2];
	  
	  for(j = 0; j < 3; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
		  r[j][i] -= h;
	    }
	  }

	  StretchCrossForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	  }
	}
	break;

      /*********************
       * Water Force	   *
       *********************/
      case 17:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  par[1] = (species->Pots[ipot]).p[iintra][1];
	  
	  for(j = 0; j < 3; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
		  r[j][i] -= h;
	    }
	  }

	  WaterForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	  }
	}
	break;

      /*********************
       * Ammonia Force	   *
       *********************/
      case 18:
	for( iintra = 0; iintra < (species->Pots[ipot]).n_intact; iintra++)
	{
	  par[0] = (species->Pots[ipot]).p[iintra][0];
	  
	  for(j = 0; j < 4; j++)
	  {
	    issite[j] = ispec * species->nsites
	                + (species->Pots[ipot]).SiteArr[iintra][j] - 1;
	    
	    for ( i = 0; i < 3; i++)
	    {
	       r[j][i] = site[i][issite[j]];	  
	       if (j > 0)
		  if (r[0][i] > 0.0) 
		  {
		     if ( r[0][i] - r[j][i] > h_2 )
			r[j][i] += h;
		  }
		  else
		     if ( r[j][i] - r[0][i] > h_2 )
		  r[j][i] -= h;
	    }
	  }

	  AmmoniaForce( r, f, par, pot);
	  for (i = 0; i < 3; i++)
	  {
	    site_force[i][issite[0]] += f[0][i];
	    site_force[i][issite[1]] += f[1][i];
	    site_force[i][issite[2]] += f[2][i];
	    site_force[i][issite[3]] += f[3][i];
	  }
	}
	break;

	default : 	/* should never occur            */
	                /* HIER FEHLERMELDUNG EINBAUEN ! */
	 message(NULLI, NULLP, FATAL, INTERR, __LINE__, __FILE__);
	break;
      }
    }
   species->intrpot[ispec] += *pot;	
  }
}
