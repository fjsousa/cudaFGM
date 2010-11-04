#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fireCudaLib.h"

int F_SlopeAspect(double*, int, int, double, double);
double* F_reSize(int, int, double, double);

int createMaps(int N, int M, double W, double H)
{
	double* reSized;

	reSized = F_reSize(N, M, W, H);	
	F_SlopeAspect(reSized, N, M, W, H);

	return (0);
}

double* F_reSize(int N, int M, double W, double H)
{

	int n, m, row, col, index;
	int N1, N2, N3;
	double Dx, Dy, DX, DY, X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3, a, b, c; 
	double* original, *midStep, *reSized;
	FILE *out, *original_file;
	
	original_file = fopen("sEstrela.map", "r");	
	fscanf(original_file, "%d %d", &n, &m);

	Dx = W / (m-1);
	Dy = H / (n-1);
	DX = W / (M-1);
	DY = H / (N-1);

	original = (double*) malloc(n*m*sizeof(double));
	midStep  = (double*) malloc(n*M*sizeof(double));
	reSized  = (double*) malloc(N*M*sizeof(double));

	//Read original data
	for (row = 0; row < n; row++)
	{
		for (col = 0; col < m; col++)
		{
			fscanf(original_file, "%lf", &original[col + m*row]);
			original[col + m*row] *=4;
		}
	}
	
	//X direction
	//copy borders
	for (row = 0; row < n; row++)
	{
		midStep[0 + row*M]   = original[0 + row*m];
		midStep[M-1 + row*M] = original[m-1 + row*m];
	}


	if (m < M)
	{
		for (row = 0; row < n; row++ )
		{
			for (col = 1; col < M-1; col++ )
			{ 
				N1 = DX/Dx*col -1 ;
				N2 = N1+1;
				N3 = N1+2;

				X1 = Dx * N1;
				X2 = Dx * N2;
				X3 = Dx * N3;

				Z1 = original[N1 + m*row];
				Z2 = original[N2 + m*row];
				Z3 = original[N3 + m*row];

				a = -( -X2*Z1 + X3*Z1 + X1*Z2 - X3*Z2 - X1*Z3 + X2*Z3 )/
						(X2 - X1)/(X2 - X3)/(X3 - X1);

				//if (a < smidgen) a = 0.;
				
				b = -( -Z2*X1*X1 + Z3*X1*X1 + Z1*X2*X2 - Z1*X3*X3 + Z2*X3*X3 - X2*X2*Z3 )/
						(X1 - X2)/(X1 -  X3)/(X2 - X3);
				
				//if (b < smidgen) b = 0.;

				c = - ( X3*Z2*X1*X1 - X2*Z3*X1*X1 - X3*X3*Z2*X1 + X2*X2*X1*Z3 + X3*X3*X2*Z1 - X2*X2*X3*Z1)
						/(X2 - X3)/(X1*X1 - X2*X1 - X3*X1 + X2*X3);
				
				//if (c < smidgen) c = 0.;

				midStep[col + row*M] = a*col*DX*col*DX + b*col*DX + c;

			}
		}
	}
	else if (m > M)
	{
		for (row = 0; row < n; row++ )
		{
			for (col = 1; col < M-1; col++ )
			{
				index = (int) DX/Dx*col;
				midStep[col + row*M] = original[index + m*row];
			}
		}
	}
	else 
	{	
		for (row = 0; row < n; row++ )
			for (col = 1; col < M-1; col++ )
				midStep[col + row*M] =  original[col + row*M];
	}
	
	//Y direction
	//copy borders
	for (col = 0; col < M; col++)
	{
		reSized[col]   				 = midStep[col];
		reSized[col + (N-1)*M] = midStep[col + (n-1)*M];
	}

	if (n < N)
	{
		for (row = 1; row < N-1; row++ )
		{
			for (col = 0; col < M; col++ )
			{ 
				N1 = DY/Dy*row -1 ;
				N2 = N1+1;
				N3 = N1+2;

				Y1 = Dy * N1;
				Y2 = Dy * N2;
				Y3 = Dy * N3;

				Z1 = midStep[col + N1*M];
				Z2 = midStep[col + N2*M];
				Z3 = midStep[col + N3*M];

				a = -( -Y2*Z1 + Y3*Z1 + Y1*Z2 - Y3*Z2 + -Y1*Z3 + Y2*Z3 )/
						(Y2 - Y1)/(Y2 -  Y3)/(Y3 - Y1);
				
				//if (a < smidgen) a = 0.;

				b = -( -Z2*Y1*Y1 + Z3*Y1*Y1 + Z1*Y2*Y2 - Z1*Y3*Y3 + Z2*Y3*Y3 - Y2*Y2*Z3 )/
						(Y1 - Y2)/(Y1 -  Y3)/(Y2 - Y3);

				//if (b < smidgen) b = 0.;
				
				c = - ( Y3*Z2*Y1*Y1 - Y2*Z3*Y1*Y1 - Y3*Y3*Z2*Y1 + Y2*Y2*Y1*Z3 + Y3*Y3*Y2*Z1 - Y2*Y2*Y3*Z1)
						/(Y2 - Y3)/(Y1*Y1 - Y2*Y1 - Y3*Y1 + Y2*Y3);

				//if (c < smidgen) c = 0.;

				reSized[col + row*M] = a*row*DY*row*DY + b*row*DY + c;

			}
		}
	}
	else if (n > N)
	{
		for (row = 1; row < N-1; row++ )
		{
			for (col = 0; col < M; col++ )
			{
				index = (int) DY/Dy*row;
				reSized[col + row*M] = midStep[col + index*M];
			}
		}
	}
	else 
	{	
		for (row = 1; row < N-1; row ++ )
			for (col = 0; col < M; col ++ )
				reSized[col + row*M] =  midStep[col + row*M];
	}

	//print
	out = fopen("reSized.map", "w" );
	for (row = 0; row < N; row ++ )
	{
		for (col = 0; col < M; col ++ )
		{
			fprintf(out, " %f ", reSized[col + row*M]);
		}
			fprintf(out, "\n");
	}
	/*
	out = fopen("reSized.vtk", "w" );
	fprintf(out, "# vtk DataFile Version 2.0\n");
	fprintf(out, "Height\n");
	fprintf(out, "ASCII\n");
	fprintf(out, "DATASET POLYDATA\n");
	fprintf(out, "POINTS %d double\n",M*N); 
	for (row = 0; row < N; row ++ )
	{
		for (col = 0; col < M; col ++ )
		{
			fprintf(out, "%f %f %f", DX*col, DY*row, reSized[col + row*M]);
			fprintf(out, "\n");
		}
	}
	*/

	fclose(out);
	fclose(original_file);
	free(original);
	free(midStep);
	
	return(reSized);
	
}
	
////////////////////////////////////////////////////////////////////////////////////////////
//Slope and Aspect
////////////////////////////////////////////////////////////////////////////////////////////	
int F_SlopeAspect(double* reSized, int N, int M, double W, double H)
{
	int n, m, row, col, min_n;
	double* slope, *aspect;
	int P_row[4], P_col[4];
	double P_h[4];
	double DX, DY, Hmed, min, A, B, C, x1, x2, x3, y1, y2, y3, z1, z2, z3, angle; 
	FILE* original_file, *out;
	
	original_file = fopen("sEstrela.map", "r");	
	fscanf(original_file, "%d %d", &n, &m);

	DX = W / (M-1);
	DY = H / (N-1);

	slope  = (double*) malloc((N-1)*(M-1)*sizeof(double));
	aspect = (double*) malloc((N-1)*(M-1)*sizeof(double));

	for (row = 0; row < N-1; row++ )
	{
		for(col = 0; col < M-1; col++)
		{
			//Find for points that must be fited in the plane
			P_row[0] = row;   		P_col[0] = col;       P_h[0] = reSized[P_col[0] + M*P_row[0]];
			P_row[1] = row;   		P_col[1] = col+1;     P_h[1] = reSized[P_col[1] + M*P_row[1]];
			P_row[2] = row+1; 		P_col[2] = col;       P_h[2] = reSized[P_col[2] + M*P_row[2]];
			P_row[3] = row+1; 		P_col[3] = col+1;     P_h[3] = reSized[P_col[3] + M*P_row[3]];

			Hmed = (P_h[0] + P_h[1] + P_h[2] + P_h[3]) / 4;

			//find the one point that is close to the averadge
			min = INF;
			for (n = 0; n < 3; n++)
			{
				if (abs(P_h[n] - Hmed) < min) 
				{
					min   = abs(P_h[n] - Hmed); 
					min_n = n;
				}
			}

			//Define the plane, excluding that point - the excluded point is reseted by  4 position in the point array
			for (n = 0; n < 2; n++)
			{
				if ( n == min_n) 
				{
					P_row[n] = P_row[3];
					P_col[n] = P_col[3];
					P_h[n]   = P_h[3];
				}
			}
			
			y1 = DY*P_row[0];     x1 = DX*P_col[0];   z1 = P_h[0];
			y2 = DY*P_row[1];     x2 = DX*P_col[1];   z2 = P_h[1];
			y3 = DY*P_row[2];     x3 = DX*P_col[2];   z3 = P_h[2];

			A = y1 *(z2 - z3) + y2 *(z3 - z1) + y3 *(z1 - z2);
			B = z1 *(x2 - x3) + z2 *(x3 - x1) + z3 *(x1 - x2);
			C = x1 *(y2 - y3) + x2 *(y3 - y1) + x3 *(y1 - y2);
			//D =  - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)

			//horizontal plane is only z = 0, or c = 1
			angle = acos(C/sqrt(A*A + B*B + C*C));

			slope[col + row*(M-1)] =  fabs(tan(angle)); 

			//aspect is angle between (A,B,C) projection in the horizontal -> (A,B,0) and (0,-1,0)	
			
			//Get simetric normal vector, if vector is pointing "upward"
			if (C > 0)
			{
				C = C*-1;
				A = A*-1;
				B = B*-1;
			}
			//Special Cases
			if (A == 0 && B > 0)
				angle = 0.;
			
			else if (A > 0 && B == 0)
				angle = 270.;

			else if (A == 0 && B < 0)
				angle = 180.;
			
			else if (A < 0 && B == 0)
				angle = 90.;
			
			if ( A > 0  && B < 0) //1st quadrant 
			{	
				angle = atan( A / B);
				angle = RadToDeg( fabs( angle ) );
				angle = angle + 180;
			}

			if ( A > 0  && B > 0) //2st quadrant 
			{
				angle = atan( A / B);
				angle = 360. - RadToDeg( fabs( angle )) ;
			}

			if ( A < 0  && B > 0) //3st quadrant 
			{
				angle = atan( A / B);
				angle = RadToDeg( fabs( angle ) );
			}

			if ( A < 0  && B < 0) //4st quadrant 
			{
				angle = atan( A / B);
				angle = 180. - RadToDeg( fabs( angle ));
			}

			aspect[col + row*(M-1)] = angle;
		
		}
	}

	out = fopen("slope.map", "w" );
	for (row = 0; row < N-1; row ++ )
	{
		for (col = 0; col < M-1; col ++ )
		{
			fprintf(out, " %f ", slope[col + row*(M-1)]);
		}
			fprintf(out, "\n");
	}

	out = fopen("aspect.map", "w" );
	for (row = 0; row < N-1; row ++ )
	{
		for (col = 0; col < M-1; col ++ )
		{
			fprintf(out, " %f ", aspect[col + row*(M-1)]);
		}
			fprintf(out, "\n");
	}
	
	fclose(out);
	fclose(original_file);
	free(aspect);
	free(slope);
	free(reSized);	
	
	return(0);
}
