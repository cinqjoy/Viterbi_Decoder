#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long *idum;

void normal(double *n1,double *n2,double sigma);
double ran1(void);

int main(int argc, char *argv[])
{
    int N,SEED,HorS;
    float SNRdB;
    
    printf("N= ");
    scanf("%d", &N);
    printf("SEED= ");
    scanf("%d", &SEED);
    printf("SNR= ");
    scanf("%f", &SNRdB);
    printf("1)SOFT 2)HARD DECISION? ");
    scanf("%d", &HorS);

    idum=(long*)malloc(sizeof(long));
    *idum=SEED; 
    double sigma=1.0/sqrt(pow(10,SNRdB/10));
    
    //placing the state
    int s[64][6],i; 
    for(i=0;i<64;i++){ 
        s[i][0]=i%2; 
        s[i][1]=(i/2)%2; 
        s[i][2]=(i/4)%2; 
        s[i][3]=(i/8)%2; 
        s[i][4]=(i/16)%2; 
        s[i][5]=(i/32)%2; 
        }
    
    //generate the information bits
    int *u;
    u=(int*)malloc(sizeof(int)*(N+32));
    *u=1;
    *(u+1)=*(u+2)=*(u+3)=*(u+4)=*(u+5)=0;
    for(i=0;i<(N+32-6);i++)
        u[i+6]=u[i+1]^u[i];
    
    //encode the information sequence using the generator matrix
    int *x1,*x2,st=0;
    x1=(int*)malloc(sizeof(int)*(N+32));
    x2=(int*)malloc(sizeof(int)*(N+32));
    for(i=0;i<N+32;i++){
        x1[i]=u[i]^s[st][1]^s[st][2]^s[st][4]^s[st][5]; 
        x2[i]=u[i]^s[st][0]^s[st][1]^s[st][2]^s[st][5];
        st=u[i]+s[st][0]*2+s[st][1]*4+s[st][2]*8+s[st][3]*16+s[st][4]*32;
        }   
        
    int time,st0,st1,v01,v02,v11,v12,j,k,output,*decodeu,error=0;
    decodeu=(int*)malloc(sizeof(int)*(N)); 
    int surv[64][32]={-1};
    double y1,y2,n1,n2;  
  
    if(HorS==1){
    //soft decision
    double metric[64][32]={-1};
    for(time=0;time<6;time++){                 
        //channel input:map 0 to +1 and 1 to -1
        if(x1[time]==0) y1=1.0;
           else y1=-1.0;
        if(x2[time]==0) y2=1.0;
           else y2=-1.0;
        
        //add normal r.v. to the +1,-1
        normal(&n1,&n2,sigma);
        y1+=n1;
        y2+=n2; 
           
        for(i=0;i<(pow(2,time));i++){                       
            //next state
            st0=0+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
            st1=1+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
               
            //encode msg 0
            v01=0^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v02=0^s[i][0]^s[i][1]^s[i][2]^s[i][5];
            if(v01==0) v01=1;
               else v01=-1;
            if(v02==0) v02=1;
               else v02=-1;
               
            //encode msg 1
            v11=1^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v12=1^s[i][0]^s[i][1]^s[i][2]^s[i][5];
            if(v11==0) v11=1;
               else v11=-1;
            if(v12==0) v12=1;
               else v12=-1;
               
            //save the metric and survivor
            metric[st0][time]=v01*y1+v02*y2;
            metric[st1][time]=v11*y1+v12*y2;
            surv[st0][time]=i;
            surv[st1][time]=i;   
            if(time!=0){
               metric[st0][time]+=metric[i][time-1];
               metric[st1][time]+=metric[i][time-1];    
               }
            }                    
        }//for(time)
        
    for(time=6;time<31;time++){
        //channel input:map 0 to +1 and 1 to -1
        if(x1[time]==0) y1=1.0;
           else y1=-1.0;
        if(x2[time]==0) y2=1.0;
           else y2=-1.0;
        
        //add normal r.v. to the +1,-1
        normal(&n1,&n2,sigma);
        y1+=n1;
        y2+=n2;

        for(i=0;i<64;i++){           
            //next state
            st0=0+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
            st1=1+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
               
            //encode msg 0
            v01=0^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v02=0^s[i][0]^s[i][1]^s[i][2]^s[i][5];
            if(v01==0) v01=1;
               else v01=-1;
            if(v02==0) v02=1;
               else v02=-1;
                  
            //encode msg 1
            v11=1^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v12=1^s[i][0]^s[i][1]^s[i][2]^s[i][5];
            if(v11==0) v11=1;
               else v11=-1;
            if(v12==0) v12=1;
               else v12=-1;
               
            //save the metric and survivor(add)
            //compare and chose one to be the final survivor 
            if(v01*y1+v02*y2+metric[i][time-1]>metric[st0][time]){
               metric[st0][time]=v01*y1+v02*y2+metric[i][time-1];
               surv[st0][time]=i;
               }
            if(v11*y1+v12*y2+metric[i][time-1]>metric[st1][time]){
               metric[st1][time]=v11*y1+v12*y2+metric[i][time-1];
               surv[st1][time]=i;
               } 
            }
        }//for(time)
        
    double max;
    for(time=31;time<N+32;time++){
        //channel input:map 0 to +1 and 1 to -1
        if(x1[time]==0) y1=1.0;
           else y1=-1.0;
        if(x2[time]==0) y2=1.0;
           else y2=-1.0;
        
        //add normal r.v. to the +1,-1
        normal(&n1,&n2,sigma);
        y1+=n1;
        y2+=n2;
        
        for(i=0;i<64;i++){
            //next state
            st0=0+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
            st1=1+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
               
            //encode msg 0
            v01=0^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v02=0^s[i][0]^s[i][1]^s[i][2]^s[i][5];
            if(v01==0) v01=1;
               else v01=-1;
            if(v02==0) v02=1;
               else v02=-1;
                  
            //encode msg 1
            v11=1^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v12=1^s[i][0]^s[i][1]^s[i][2]^s[i][5];
            if(v11==0) v11=1;
               else v11=-1;
            if(v12==0) v12=1;
               else v12=-1;
               
            //compare and save the final survivor 
            if(v01*y1+v02*y2+metric[i][30]>metric[st0][31]){
               metric[st0][31]=v01*y1+v02*y2+metric[i][30];
               surv[st0][31]=i;
               }
            if(v11*y1+v12*y2+metric[i][30]>metric[st1][31]){
               metric[st1][31]=v11*y1+v12*y2+metric[i][30];
               surv[st1][31]=i;
               } 
            }
                      
        //find the best path 
        max=metric[0][31];
        for(i=0;i<64;i++)   
            if(metric[i][31]>max) 
               max=metric[i][31];
                      
        //trace back, save output msg   
        for(i=0;i<64;i++){
            if(metric[i][31]==max){
               k=i;
               for(j=31;j>0;j--)
                   k=surv[k][j];
               output=s[k][0];                      
               break;
               }
            }
        /*
        //fixed state
        k=0;
        for(j=31;j>0;j--)
            k=surv[k][j];
        output=s[k][0]; 
        
        //majority vote
        int outputtemp[64]={-1};    
        for(i=0;i<64;i++){ 
            k=i;
            for(j=31;j>0;j--)
                k=surv[k][j];
                outputtemp[i]=s[k][0];                      
            }
        int op0=0,op1=0;
        for(i=0;i<64;i++){
            if(outputtemp[i]==0) op0++;
            if(outputtemp[i]==1) op1++;
            }
        if(op0>=op1) output=0;
        else output=1;
        */
        //clean and update the metric   
        for(j=1;j<32;j++){ 
            for(i=0;i<64;i++){ 
                metric[i][j-1]=metric[i][j]; 
                surv[i][j-1]=surv[i][j];
                if(j==31){
                   metric[i][31]=0; 
                   surv[i][31]=-1;
                   }  
                } 
            }
           
        //ouput msg   
        *(decodeu+time-31)=output;
        
        //cal error
        error+=*(decodeu+time-31)^*(u+time-31);
        }//for(time)
    }
-----------------------------------------------------------------------------------------------------------------
    if(HorS==2){
    //hard decision
    int Y[2],metric[64][32]={0},reference[64]={0};
    for(time=0;time<6;time++){                 
        //channel input:map 0 to +1 and 1 to -1
        if(x1[time]==0) y1=1.0;
           else y1=-1.0;
        if(x2[time]==0) y2=1.0;
           else y2=-1.0;
        
        //add normal r.v. to the +1,-1
        normal(&n1,&n2,sigma);
        y1+=n1;
        y2+=n2; 
        
        //hard decision
        if(y1>0) Y[0]=1; 
           else Y[0]=-1;
        if(y2>0) Y[1]=1;
           else Y[1]=-1;
        
        //map 1 to 0 and -1 to 1        
        if(Y[0]==1) Y[0]=0;
           else Y[0]=1;
        if(Y[1]==1) Y[1]=0;
           else Y[1]=1;
           
        for(i=0;i<(pow(2,time));i++){                       
            //next state
            st0=0+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
            st1=1+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
               
            //encode msg 0
            v01=0^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v02=0^s[i][0]^s[i][1]^s[i][2]^s[i][5];
               
            //encode msg 1
            v11=1^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v12=1^s[i][0]^s[i][1]^s[i][2]^s[i][5];
               
            //save the metric and survivor
            metric[st0][time]=(v01^Y[0])+(v02^Y[1]);
            metric[st1][time]=(v11^Y[0])+(v12^Y[1]);
            surv[st0][time]=i;
            surv[st1][time]=i;  
            if(time!=0){
               metric[st0][time]+=metric[i][time-1]; 
               metric[st1][time]+=metric[i][time-1];     
               }
            }                    
        }//for(time)
      
    for(time=6;time<31;time++){
        //channel input:map 0 to +1 and 1 to -1
        if(x1[time]==0) y1=1.0;
           else y1=-1.0;
        if(x2[time]==0) y2=1.0;
           else y2=-1.0;
        
        //add normal r.v. to the +1,-1
        normal(&n1,&n2,sigma);
        y1+=n1;
        y2+=n2;

        //hard decision
        if(y1>0) Y[0]=1; 
           else Y[0]=-1;
        if(y2>0) Y[1]=1;
           else Y[1]=-1;
        
        //map 1 to 0 and -1 to 1        
        if(Y[0]==1) Y[0]=0;
           else Y[0]=1;
        if(Y[1]==1) Y[1]=0;
           else Y[1]=1;
          
        for(i=0;i<64;i++){           
            //next state
            st0=0+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
            st1=1+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
               
            //encode msg 0
            v01=0^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v02=0^s[i][0]^s[i][1]^s[i][2]^s[i][5];
                  
            //encode msg 1
            v11=1^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v12=1^s[i][0]^s[i][1]^s[i][2]^s[i][5];
             
            //save the metric and survivor(add)
            //compare and chose one to be the final survivor 
            if(reference[st0]==0){
               metric[st0][time]=(v01^Y[0])+(v02^Y[1])+metric[i][time-1];
               surv[st0][time]=i;
               reference[st0]++;
               } 
            else if((v01^Y[0])+(v02^Y[1])+metric[i][time-1]<metric[st0][time]){
               metric[st0][time]=(v01^Y[0])+(v02^Y[1])+metric[i][time-1];
               surv[st0][time]=i;
               }
            if(reference[st1]==0){
               metric[st1][time]=(v11^Y[0])+(v12^Y[1])+metric[i][time-1];
               surv[st1][time]=i;
               reference[st1]++;
               }
            else if((v11^Y[0])+(v12^Y[1])+metric[i][time-1]<metric[st1][time]){
               metric[st1][time]=(v11^Y[0])+(v12^Y[1])+metric[i][time-1];
               surv[st1][time]=i;
               } 
            }
            
        //clean reference   
        for(k=0;k<64;k++)
            reference[k]=0;
        }//for(time)
        
    int min;
    for(time=31;time<N+32;time++){
        //channel input:map 0 to +1 and 1 to -1
        if(x1[time]==0) y1=1.0;
           else y1=-1.0;
        if(x2[time]==0) y2=1.0;
           else y2=-1.0;
        
        //add normal r.v. to the +1,-1
        normal(&n1,&n2,sigma);
        y1+=n1;
        y2+=n2;
        
        //hard decision
        if(y1>0) Y[0]=1; 
           else Y[0]=-1;
        if(y2>0) Y[1]=1;
           else Y[1]=-1;
        
        //map 1 to 0 and -1 to 1        
        if(Y[0]==1) Y[0]=0;
           else Y[0]=1;
        if(Y[1]==1) Y[1]=0;
           else Y[1]=1;
           
        for(i=0;i<64;i++){
            //next state
            st0=0+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
            st1=1+s[i][0]*2+s[i][1]*4+s[i][2]*8+s[i][3]*16+s[i][4]*32;
              
            //encode msg 0
            v01=0^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v02=0^s[i][0]^s[i][1]^s[i][2]^s[i][5];
                  
            //encode msg 1
            v11=1^s[i][1]^s[i][2]^s[i][4]^s[i][5]; 
            v12=1^s[i][0]^s[i][1]^s[i][2]^s[i][5];
     
            //save the metric and survivor(add)   
            //compare and save the final survivor 
            if(reference[st0]==0){
               metric[st0][31]=(v01^Y[0])+(v02^Y[1])+metric[i][30];
               surv[st0][31]=i;
               reference[st0]++;
               } 
            else if((v01^Y[0])+(v02^Y[1])+metric[i][30]<metric[st0][31]){
               metric[st0][31]=(v01^Y[0])+(v02^Y[1])+metric[i][30];
               surv[st0][31]=i;
               }
            if(reference[st1]==0){
               metric[st1][31]=(v11^Y[0])+(v12^Y[1])+metric[i][30];
               surv[st1][31]=i;
               reference[st1]++;
               }
            else if((v11^Y[0])+(v12^Y[1])+metric[i][30]<metric[st1][31]){
               metric[st1][31]=(v11^Y[0])+(v12^Y[1])+metric[i][30];
               surv[st1][31]=i;
               }  
            }
         
        //find the best path 
        min=metric[0][31];
        for(i=0;i<64;i++)  
            if(metric[i][31]<min) 
               min=metric[i][31];
                 
        //trace back, save output msg   
        for(i=0;i<64;i++){ 
            if(metric[i][31]==min){
               k=i;
               for(j=31;j>0;j--)
                   k=surv[k][j];
               output=s[k][0];                      
               break;
               }
            }
        /*
        //fixed state
        k=0;
        for(j=31;j>0;j--)
            k=surv[k][j];
        output=s[k][0]; 
        
        //majority vote
        int outputtemp[64]={-1};    
        for(i=0;i<64;i++){ 
            k=i;
            for(j=31;j>0;j--)
                k=surv[k][j];
                outputtemp[i]=s[k][0];                      
            }
        int op0=0,op1=0;
        for(i=0;i<64;i++){
            if(outputtemp[i]==0) op0++;
            if(outputtemp[i]==1) op1++;
            }
        if(op0>=op1) output=0;
        else output=1;
        */
        //clean and update the metric   
        for(j=1;j<32;j++){ 
            for(i=0;i<64;i++){ 
                metric[i][j-1]=metric[i][j]; 
                surv[i][j-1]=surv[i][j];
                reference[i]=0;
                if(j==31){
                   metric[i][31]=0; 
                   surv[i][31]=-1;
                   }  
                } 
            }
        
        //ouput msg   
        *(decodeu+time-31)=output;
       
        //cal error
        error+=*(decodeu+time-31)^*(u+time-31);
        }//for(time)
    }
    printf("\n# Error: %d\n", error);
    printf("BER: %f\n\n", (float)(error*1.0/N));

    system("PAUSE");	
    return 0;
}

void normal(double *n1,double *n2,double sigma){
	double x1,x2,s;
	do{
	  x1=ran1();
	  x2=ran1();
	  x1=2*x1-1;
	  x2=2*x2-1;
	  s=x1*x1+x2*x2;
	}while(s>=1.0);
	*n1 = sigma*x1*sqrt(-2*log(s)/s);
    *n2 = sigma*x2*sqrt(-2*log(s)/s);
}

double ran1(void){
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy){
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--){
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
