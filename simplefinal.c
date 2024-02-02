#include<stdio.h>
#include<math.h>
#include<omp.h> 
typedef double newvar[30][30];

//pressure correction for lid driven cavity usinf fdm at RE = 100

int main(){
    int xgridsize = 21;
    int ygridsize = 21;
    float dx;
    double dy;
    double dt = 0.01;
    double rho = 1000;
    double mew = 10;
    newvar pfeild;
    newvar ufeild;
    newvar vfeild;
    newvar unewfeild;
    newvar vnewfeild;
    newvar pcorr;
    newvar pcorrnew;
    newvar ufinal;
    newvar vfinal;
    newvar pfinal;
    int iterations = 0;
    double pnorm = 1;
    double dpnorm = 1;
    int i,j,k,n;
    double v1,v2,u1,u2;
    double a, b, c, d;
    int m = 0;
    double temp;
    double alpha;
    double t1,t2;
    dy = 1.0/(ygridsize-1);
    dx = 1.0/(xgridsize-1);
    t1 = omp_get_wtime();

    //initiallising the required arrays 
    for(i=0;i<ygridsize + 1;i++){
        for(j=0;j<xgridsize;j++){
            ufeild[i][j] = 0;
            unewfeild[i][j] = 0;
        }
    }

    for(i=0;i<ygridsize;i++){
        for(j=0;j<xgridsize + 1;j++){
            vfeild[i][j] = 0;
            vnewfeild[i][j] = 0;
        }
    }

    for(i=0;i<ygridsize + 1;i++){
        for(j=0;j<xgridsize + 1;j++){
            pcorr[i][j] = 0;
            pcorrnew[i][j] = 0;
            pfeild[i][j] = 1;
        }
    }

    for(i=1;i<xgridsize-1;i++){
            ufeild[0][i] = 2 - ufeild[1][i];
            ufeild[ygridsize][i] = - ufeild[ygridsize - 1][i];
        }

    //non linear iterative loop 
    for(;;){
        iterations = iterations + 1;
        //calculating u* 
        for(i=1;i<ygridsize;i++){
            for(j=1;j<xgridsize - 1;j++){
                v1 = 0.5 * (vfeild[i-1][j] + vfeild[i-1][j+1]);
                v2 = 0.5 * (vfeild[i][j] + vfeild[i][j+1]);
                a =  -(((rho * ufeild[i][j+1]*ufeild[i][j+1] - rho * ufeild[i][j-1]*ufeild[i][j-1] ) / (2*dx)) + (((rho * v1 * ufeild[i-1][j]) - (rho*v2*ufeild[i+1][j])) / (2*dy))) + mew*(((ufeild[i][j+1]-2*ufeild[i][j]+ufeild[i][j-1])/(dx*dx)) + ((ufeild[i+1][j]-2*ufeild[i][j]+ufeild[i-1][j])/(dy*dy)));
                unewfeild[i][j] = (1/rho)*(rho * ufeild[i][j] + a * dt - (dt/dx)*(pfeild[i][j+1]-pfeild[i][j]));
                
            }
        }

        //applying boundary conditions 

        for(i=1;i<xgridsize-1;i++){
            unewfeild[0][i] = 2 - unewfeild[1][i];
            unewfeild[ygridsize][i] = - unewfeild[ygridsize - 1][i];
        }

        //bug code
        /*for(i=0;i<ygridsize + 1;i++){
            for(j=0;j<xgridsize;j++){
                printf("%f\t",unewfeild[i][j]);
            }
            printf("\n");
        }*/

        //calculating v* 
        for(i=1;i<ygridsize - 1;i++){
            for(j=1;j<xgridsize;j++){
                u1 = 0.5 * (ufeild[i][j-1] + ufeild[i+1][j-1]);
                u2 = 0.5 * (ufeild[i][j] + ufeild[i+1][j]);
                b = - ((( (rho*u2*vfeild[i][j+1]) - (rho*u1*vfeild[i][j-1]) )/(2*dx)) + (( (rho * vfeild[i-1][j]*vfeild[i-1][j]) - (rho * vfeild[i+1][j]*vfeild[i+1][j]))/(2*dy))) + mew*(((vfeild[i][j+1] - 2*vfeild[i][j] + vfeild[i][j-1])/(dx*dx))+((vfeild[i+1][j] - 2*vfeild[i][j] + vfeild[i-1][j])/(dy*dy)));
                vnewfeild[i][j] = vfeild[i][j] + (1/rho)*(b*dt -(dt/dx)*(pfeild[i][j]-pfeild[i+1][j]));
            }
        }

        // applying boundary conditions 

        for(i=1;i<ygridsize-1;i++){
            vnewfeild[i][0] = - vnewfeild[i][1];
            vnewfeild[i][xgridsize] = - vnewfeild[i][xgridsize - 1];
        }

        /*for(i=0;i<ygridsize;i++){
            for(j=0;j<xgridsize + 1;j++){
                printf("%f\t",vnewfeild[i][j]);
            }
            printf("\n");
        }
        printf("\n");*/

        //solving for pressure correction 
        for(;;){
            m = m + 1;

            // equating pcorr to pcorrnew 
            for(i=0;i<ygridsize + 1;i++){
                for(j=0;j<xgridsize + 1;j++){
                    pcorr[i][j] = pcorrnew[i][j];
                }
            }

            //running the pressure equation 
            for(i=1;i<ygridsize;i++){
                for(j=1;j<xgridsize;j++){
                    a = 2 * (((dt)/(dx*dx)) + ((dt)/(dy*dy)));
                    b = -((dt)/(dx*dx));
                    c = -((dt)/(dy*dy));
                    d = (1.0/dx) * ((rho * unewfeild[i][j]) - (rho*unewfeild[i][j-1])) + (1.0/dy)*((rho*vnewfeild[i-1][j])- (rho*vnewfeild[i][j]));
                    pcorrnew[i][j] = (-1/a) * ( b * pcorrnew[i][j+1] + b * pcorrnew[i][j-1] + c * pcorrnew[i+1][j] + c*pcorrnew[i-1][j] + d);
                }
            }
            //calculating the dpnorm 
            dpnorm = 0;
            for(i=0;i<ygridsize + 1;i++){
                for(j=0;j<xgridsize + 1;j++){
                    temp = (pcorrnew[i][j] - pcorr[i][j]) * (pcorrnew[i][j] - pcorr[i][j]);
                    dpnorm = dpnorm + temp;
                }
            }
            pnorm = 0;
            for(i=0;i<ygridsize + 1;i++){
                for(j=0;j<xgridsize + 1;j++){
                    temp = pcorrnew[i][j] * pcorrnew[i][j] ;
                    pnorm = pnorm + temp;
                }
            }
            dpnorm = pow(dpnorm,0.5);
            pnorm = pow(pnorm,0.5);
            dpnorm = dpnorm / pnorm ; 
            //printf("The dpnorm is %f \n",dpnorm);


            //break condition 
            if(dpnorm < 0.000001 ){
                break;
            }
        }
        //update the pressure values 
        alpha = 0.1;
        for(i=0;i<ygridsize + 1;i++){
            for(j=0;j<xgridsize + 1;j++){
                pfeild[i][j] = pfeild[i][j] + alpha*pcorrnew[i][j];
            }
        }

        //setting newman boundary conditions 
        for(i=0;i<xgridsize + 1;i++){
            pfeild[0][i] = pfeild[1][i];
            pfeild[ygridsize][i] = pfeild[ygridsize - 1][i];
        }

        for(i=0;i<ygridsize + 1;i++){
            pfeild[i][0] = pfeild[i][1];
            pfeild[i][xgridsize] = pfeild[i][xgridsize - 1];
        }

        //updating u velocity 
        for(i=0;i<ygridsize + 1;i++){
            for(j=0;j<xgridsize;j++){
                ufeild[i][j] = unewfeild[i][j];
            }   
        }
        //updating v velocity 
        for(i=0;i<ygridsize;i++){
            for(j=0;j<xgridsize + 1;j++){
                vfeild[i][j] = vnewfeild[i][j];
            }
        }

        pnorm = 0;
        for(i=0;i<ygridsize + 1;i++){
            for(j=0;j<xgridsize + 1;j++){
                temp = pcorrnew[i][j] * pcorrnew[i][j] ;
                pnorm = pnorm + temp;
            }
        }
        pnorm = pow(pnorm,0.5);
        //printf("pnorm - %f\n",pnorm);
        if(pnorm < 0.0001){
            break;
        }
       
    }

    //updating v,u and p on the normal grid 
    for(i=0;i<ygridsize;i++){
        for(j=0;j<xgridsize;j++){
            ufinal[i][j] = 0.5 * (ufeild[i][j] + ufeild[i+1][j]);
            vfinal[i][j] = 0.5 * (vfeild[i][j] + vfeild[i][j+1]);
            pfinal[i][j] = 0.25 * (pfeild[i][j] + pfeild[i+1][j] + pfeild[i][j+1] + pfeild[i+1][j+1]);
        }
    }
    t2 = omp_get_wtime();
    printf("the final Pressure values -\n");
    for(i=0;i<ygridsize;i++){
        for(j=0;j<xgridsize;j++){
            printf("%f\t",pfinal[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("the final u velocity values -\n");
    for(i=0;i<ygridsize;i++){
        for(j=0;j<xgridsize;j++){
            printf("%f\t",ufinal[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("the final v velocity values -\n");
    for(i=0;i<ygridsize;i++){
        for(j=0;j<xgridsize;j++){
            printf("%f\t",vfinal[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("Total time taken is %g seconds\n",t2-t1);

    return 0;
}