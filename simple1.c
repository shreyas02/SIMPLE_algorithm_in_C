#include<stdio.h>
#include<math.h>
typedef double newvar[120][120];

//pressure correction for lid driven cavity usinf fdm 

void main(){
    int xgridsize,ygridsize;
    float dx,dy,dt,rho,mew;
    newvar u,v,p,dp,unew,vnew,newdp;
    int i,j,k;
    float a,b,c,d,vdash,vavg,uavg,udash;
    double unorm, vnorm, dpnorm, pnorm;
    float temp;
    float tolerance = 0.000001;

    //user input

    dx = 0.1;

    dy = 0.1;

    dt = 0.0001;

    xgridsize = 10;

    ygridsize = 10;

    rho = 1000;

    mew = 100;

    // initial conditions for pressure feild 
    for(i=0;i<ygridsize;i++){
        for(j=0;j<xgridsize;j++){
            p[i][j] = 0;
            dp[i][j]= 0;
        }
    }

    for(i = 0; i<ygridsize; i++){
        for(j=0; j<xgridsize + 1; j++){
            u[i][j] = 0;
            unew[i][j]=0;
        }
    }

    for(i = 0; i<ygridsize +1; i++){
        for(j=0; j<xgridsize; j++){
            v[i][j] = 0;
            vnew[i][j]=0;
        }
    }

    // setting boundary conditions 

    for(j=1;j<xgridsize ;j++){
            u[0][j] = 1;
            u[ygridsize -1][j] = 0;
            unew[0][j] = 1;
            unew[ygridsize -1][j] = 0;
    }

    for(j=1;j<ygridsize ;j++){
            v[j][0] = 0;
            v[j][xgridsize - 1] = 0;
            vnew[j][0] = 0;
            vnew[j][xgridsize - 1] = 0;
    }

    // code to be repeated after every update of v 
    for(i=0; i<xgridsize;i++){
            v[0][i] = - v[1][i];
            vnew[0][i] = - vnew[1][i];
    }
    for(i=0; i<xgridsize;i++){
            v[ygridsize][i] = - v[ygridsize - 1][i];
            vnew[ygridsize][i] = - vnew[ygridsize - 1][i];
    }

    // code to be repeated after every update of u 
    for(i=0; i<ygridsize;i++){
            u[i][0] = - u[i][1];
    }
    for(i=0; i<ygridsize;i++){
            u[i][xgridsize] = - u[i][xgridsize - 1];
    }

    do{

    //calculate u from pressure feild
    do{
    for(i=1;i<ygridsize - 1;i++){
        for(j=1;j<xgridsize;j++){
        vdash = 0.5*(v[i][j-1]+v[i][j]); // at point a 
        vavg = 0.5*(v[i+1][j-1]+v[i + 1][j]);  // at point b 
        a = - (((rho * u[i][j+1]*u[i][j+1] - rho * u[i][j-1]*u[i][j-1] ) / (2*dx)) + (((rho*vdash*u[i-1][j]) - (rho*vavg*u[i+1][j])) / (2*dy))) + mew*(((u[i][j+1]-2*u[i][j]+u[i][j-1])/(dx*dx)) + ((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dy*dy)));
        unew[i][j] = 1/rho*(rho * u[i][j] + a * dt - (dt/dx)*(p[i][j]-p[i][j-1])); 
        }
    }

    // code to be repeated after every update of u 
    for(i=0; i<ygridsize;i++){
            unew[i][0] = - unew[i][1];
    }
    for(i=0; i<ygridsize;i++){
            unew[i][xgridsize] = - unew[i][xgridsize - 1];
    }

    //unorm calculation
    unorm = 0;
    for(i = 0; i<ygridsize; i++){
    for(j=0; j<xgridsize + 1; j++){
        temp = unew[i][j] - u[i][j];
        unorm = unorm + temp*temp;
        }
    }
    unorm = pow(unorm,0.5);

    //update u velocity
    for(i = 0; i<ygridsize; i++){
        for(j=0; j<xgridsize + 1; j++){
            //u[i][j] = unew[i][j];
            u[i][j] = u[i][j] + 0.8*(unew[i][j] - u[i][j]); //under relaxation 
            //printf("%f\t",u[i][j]);
        }
        //printf("\n");
    }
    //printf("\n");
    printf("unorm is %0.8lf\n",unorm);
    }while(unorm > tolerance);

    //calculate v from pressure feild
    do{
        for(i=1;i<ygridsize;i++){
        for(j=1;j<xgridsize - 1;j++){
            uavg = 0.5*(u[i-1][j]+u[i][j]); //at point c
            udash = 0.5*(u[i-1][j+1]+u[i][j+1]); // at point d 
            b = - (((rho * v[i][j+1]*udash - rho * u[i][j-1]*uavg ) / (2*dx)) + (((rho*v[i-1][j]*v[i-1][j]) - (rho*v[i+1][j]*v[i+1][j])) / (2*dy))) + mew*(((v[i][j+1]-2*v[i][j]+v[i][j-1])/(dx*dx)) + ((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dy*dy)));
            vnew[i][j] = 1/rho*(rho * v[i][j] + b * dt - (dt/dy)*(p[i-1][j]-p[i][j])); 
        }
        }

    // code to be repeated after every update of v 
    for(i=0; i<xgridsize;i++){
            vnew[0][i] = - vnew[1][i];
    }
    for(i=0; i<xgridsize;i++){
            vnew[ygridsize][i] = - vnew[ygridsize - 1][i];
    }

    //vnorm calculation
    vnorm = 0;
    for(i = 0; i<ygridsize + 1; i++){
    for(j=0; j<xgridsize; j++){
        temp = vnew[i][j] - v[i][j];
        vnorm = vnorm + temp*temp;
        }
    }
    vnorm = pow(vnorm,0.5);

    //update u velocity
    for(i = 0; i<ygridsize + 1; i++){
        for(j=0; j<xgridsize; j++){
            //v[i][j] = vnew[i][j];
            v[i][j] = v[i][j] + 0.8*(vnew[i][j] - v[i][j]);
            //printf("%f\t",v[i][j]);
        }
        //printf("\n");
    }
    //printf("\n");
    printf("vnorm is %0.9lf\n",vnorm);
    }while(vnorm > tolerance);

    //pressure correction 

    //newman boundary conditions would be used for ghost cells 

    //a pressure point at top left node is taken and then pinned to 0 

    for(i = 0; i<ygridsize + 2; i++){
        for(j=0; j<xgridsize + 2; j++){
            dp[i][j] = 0;
            newdp[i][j] = 0;
        }
    }

    //code to be repeated after every dp update
    for(i=1;i<ygridsize + 1; i++){
        dp[i][0] = dp[i][1];
        dp[i][xgridsize + 1] = dp[i][xgridsize];
        newdp[i][0] = newdp[i][1];
        newdp[i][xgridsize + 1] = newdp[i][xgridsize];
    }
    for(j=1;j<xgridsize + 1; j++){
        dp[0][j] = dp[1][j];
        dp[ygridsize+1][j] = dp[ygridsize][j];
        newdp[0][j] = newdp[1][j];
        newdp[ygridsize+1][j] = newdp[ygridsize][j];
    }

    do{
    for(i=1;i<ygridsize + 1;i++){
        for(j=1;j<xgridsize + 1;j++){
            if (i != 1 || j != 1){
            a = 2 * ((dt / (dx*dx)) + (dt / (dy*dy)));
            b = - dt / (dx*dx) ;
            c = - dt / (dy*dy) ;
            d = (1/dx) * ((rho*u[i-1][j]) - (rho*u[i-1][j-1]) ) + (1/dy) * ((rho*v[i-1][j-1]) - (rho*v[i][j-1])); 
            newdp[i][j] = (-1/a) * (b*dp[i][j+1] + b*dp[i][j-1] + c*dp[i+1][j] + c*dp[i-1][j] +d);
            }
            else{
                newdp[i][j] = 0 ;
            }
        }
    }

    //code to be repeated after every dp update
    for(i=1;i<ygridsize + 1; i++){
        newdp[i][0] = newdp[i][1];
        newdp[i][xgridsize + 1] = newdp[i][xgridsize];
    }
    for(j=1;j<xgridsize + 1; j++){
        newdp[0][j] = newdp[1][j];
        newdp[ygridsize+1][j] = newdp[ygridsize][j];
    }

    /*for(i = 0; i<ygridsize + 2; i++){
        for(j=0; j<xgridsize + 2; j++){
            printf("%f\t",newdp[i][j]);
        }
        printf("\n");
    }*/

    dpnorm = 0;
    for(i = 0; i<ygridsize + 2; i++){
        for(j=0; j<xgridsize + 2; j++){
            temp = (newdp[i][j] - dp[i][j]);
            dpnorm = dpnorm + temp*temp;
        }
    }
    dpnorm = pow(dpnorm,0.5);

    //update dp 

    for(i = 0; i<ygridsize + 2; i++){
        for(j=0; j<xgridsize + 2; j++){
            dp[i][j] = newdp[i][j];
        }
    }
    printf("dpnorm is %0.8f\n",dpnorm);
    }while(dpnorm > tolerance);

    //update pressure 

    for(i = 0; i<ygridsize; i++){
        for(j=0; j<xgridsize; j++){
            p[i][j] = p[i][j] + 0.8 * dp[i+1][j+1]; //under relaxation
            printf("%f\t",p[i][j]);
        }
        printf("\n");
    }

    pnorm = 0;
    for(i = 0; i<ygridsize; i++){
        for(j=0; j<xgridsize; j++){
            pnorm = pnorm + dp[i+1][j+1]*dp[i+1][j+1];
        }
    }
    pnorm = pow(pnorm,0.5);
    printf("pnorm is %f\n",pnorm);
    }while(pnorm > tolerance);
}