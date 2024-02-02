#include<stdio.h>
#include<math.h>
typedef float newvar[7][7];

void main(){
    int xgridsize,ygridsize; 
    float dx,dy,dt;
    float a,b,c,d,rho;
    int i,j,k,n;
    newvar u,v ,p , dp , newu, newv, olddp, deltadp;
    float vdash, vavg, udash, uavg, mew, normdeltadp, normdp;

    //Set values of dx, dy and dt

    printf("Enter dx\n");
    scanf("%f",&dx);

    printf("Enter dy\n");
    scanf("%f",&dy);

    printf("Enter dt\n");
    scanf("%f",&dt);

    //accept values of size from user 

    printf("Enter the size of the grid in x direction\n");
    scanf("%d",&xgridsize);

    printf("Enter the size of the grid in y direction\n");
    scanf("%d",&ygridsize);

    printf("enter density\n");
    scanf("%f",&rho);

    printf("enter dynamic viscosity\n");
    scanf("%f",&mew);

    //enter pressure feild 

    printf("Enter the pressure feild \n");

    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize;i++){
            scanf("%f",&p[j][i]);
        }
    }

    printf("The pressure feild is \n");

    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize;i++){
            printf("%f\t",p[j][i]);
        }
        printf("\n");
    }

    //inital conditions 

    printf("Enter initial values for u\n");
    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize-1;i++){
            if(j==ygridsize-1 || i==xgridsize-1 || i==0 || j==0){
            printf("for x = %d and y = %d \n",i,j);
                    scanf("%f", &u[j][i]);                      // u will be i+1/2 , j 
                    newu[j][i] = u[j][i];
            }
            else{
                u[j][i] = 0 ;                 
                newu[j][i] = 0;
            }
        }
    }

    printf("Enter initial values for v\n");
    for(j=0;j<ygridsize-1;j++){
        for(i=0;i<xgridsize;i++){
            if(j==ygridsize-1 || i==xgridsize-1 || i==0 || j==0){
            printf("for x = %d and y = %d \n",i,j);         // v will be i,j+1/2
                    scanf("%f", &v[j][i]);
                    newv[j][i] = v[j][i];
            }
            else{
                v[j][i] = 0 ;                 
                newv[j][i] = 0;
            }
        }
    }

    printf("initial velocity u is -\n");
    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize -1 ;i++){
            printf("%f\t",u[j][i]);
        }
        printf("\n");
    }

    printf("Iitial velocity v is -\n");
    for(j=0;j<ygridsize-1;j++){
        for(i=0;i<xgridsize  ;i++){
            printf("%f\t",v[j][i]);
        }
        printf("\n");
    }

    //setting dp to 0 
    for(j=1;j<ygridsize-1;j++){
        for(i=1;i<xgridsize - 1 ;i++){
            olddp[j][i] = 0;
            dp[j][i] = 0;
        }
    }

    // starting the do while loop 

    do{
    
    // calculate u from pressure feild  

    for(k=0; k< 1000; k++){
    for(j=1;j<ygridsize-1;j++){
        for(i=1;i<xgridsize -2 ;i++){
            //jacobi being used here
            vdash = 0.5*(v[j][i]+v[j][i+1]);
            vavg = 0.5*(v[j-1][i]+v[j-1][i+1]);
            a = - (((rho * u[j][i+1]*u[j][i+1] - rho * u[j][i-1]*u[j][i-1] ) / (2*dx)) - (((rho*vdash*u[j+1][i]) - (rho*vavg*u[j-1][i])) / (2*dy))) + mew*(((u[j][i+1]-2*u[j][i]+u[j][i-1])/(dx*dx)) + ((u[j+1][1]-2*u[i][j]+u[j-1][i])/(dy*dy)));
            newu[j][i]= u[j][i] + ((a*dt)/rho) - ((dt/dx)*(p[j][i+1]-p[j][i]))/rho; 
        }
    }
    }

    // calculate v from pressure feild  

    for(k=0; k< 1000; k++){
    for(j=1;j<ygridsize-2;j++){
        for(i=1;i<xgridsize - 1 ;i++){
            udash = 0.5*(v[j][i]+v[j+1][i]);
            uavg = 0.5*(v[j][i-1]+v[j+1][i-1]);
            b = - (((rho * v[j][i+1]*udash - rho * u[j][i-1]*uavg ) / (2*dx)) - (((rho*v[j+1][i]*v[j+1][i]) - (rho*v[j-1][i]*v[j-1][i])) / (2*dy))) + mew*(((v[j][i+1]-2*v[j][i]+v[j][i-1])/(dx*dx)) + ((v[j+1][1]-2*v[i][j]+v[j-1][i])/(dy*dy)));
            newv[j][i]= v[j][i] + ((b*dt)/rho) - ((dt/dx)*(p[j][i+1]-p[j][i]))/rho; 
        }
    }
    }

    // pressure correction 
    do{
    for(j=1;j<ygridsize-1;j++){
        for(i=1;i<xgridsize - 1 ;i++){
            // gausssidel method used here 
            a = 2 * ((dt / (dx*dx)) + (dt / (dy*dy)));
            b = - dt / (dx*dx) ;
            c = - dt / (dy*dy) ;
            d = (1/dx) * ((rho*u[j][i]) - (rho*u[j][i-1]) ) + (1/dy) * ((rho*v[j][i]) - (rho*v[j-1][i])); 
            olddp[j][i] = dp[j][i];
            dp[j][i] = (-1/a) * (b*dp[j][i+1] + b*dp[j][i-1] + c*dp[j+1][i] + c*dp[j-1][i] +d);
        }
    }

    for(j=1;j<ygridsize-1;j++){
        for(i=1;i<xgridsize - 1 ;i++){
            deltadp[j][i] = dp[j][i] - olddp[j][i];
        }
    }

    // norm calculation
    normdeltadp = 0;
    for(j=1;j<ygridsize-1;j++){
        for(i=1;i<xgridsize - 1 ;i++){
            normdeltadp = normdeltadp + deltadp[j][i]*deltadp[j][i];
        }
    }
    normdeltadp = pow(normdeltadp,0.5);

    //printing norm dp
    printf("Norm deltaDp is %f\n",normdeltadp);
    } while(normdeltadp > 0.000000000001);


    //update the pressure values 

    for(j=0;j<ygridsize-1;j++){
        for(i=0;i<xgridsize - 1 ;i++){
            p[j][i] = p[j][i] + dp[j][i];
        }
    }

    //print pressure values 

    printf("Intermediate Pressure values\n");
    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize ;i++){
            printf("%f\t",p[j][i]);
        }
        printf("\n");
    }

    //update the velocity values 

    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize - 1 ;i++){
            u[j][i] = newu[j][i]; 
        }
    }

    printf("intermediate u is -\n");
    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize -1 ;i++){
            printf("%f\t",u[j][i]);
        }
        printf("\n");
    }

    for(j=0;j<ygridsize - 1 ;j++){
        for(i=0;i<xgridsize ;i++){
            v[j][i] = newv[j][i]; 
        }
    }

    // print velocity v 

    printf("intermediate v is -\n");
    for(j=0;j<ygridsize-1;j++){
        for(i=0;i<xgridsize ;i++){
            printf("%f\t",v[j][i]);
        }
        printf("\n");
    }

    //calculate dpnorm
    normdp = 0;
    for(j=1;j<ygridsize-1;j++){
        for(i=1;i<xgridsize-1 ;i++){
            normdp = normdp + dp[j][i]*dp[j][i];
        }
    }
    normdp = pow(normdp,0.5);
    printf("norm dp is %f \n",normdp*10000000000);
    }while(normdp > 0.0000000000000001);

    // print velocity u 

    printf("Final velocity u is -\n");
    for(j=0;j<ygridsize;j++){
        for(i=0;i<xgridsize -1 ;i++){
            printf("%f\t",u[j][i]);
        }
        printf("\n");
    }

    // print velocity v 

    printf("Final velocity v is -\n");
    for(j=0;j<ygridsize-1;j++){
        for(i=0;i<xgridsize  ;i++){
            printf("%f\t",v[j][i]);
        }
        printf("\n");
    }
}
