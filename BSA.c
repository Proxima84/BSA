#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#define NPLMAX 2  /* max number of stars and planets*/  
#define NTPMAX 72000 /*max number of test particles*/

double eps=1e-14; //1e-10 Morrison & Malhotra ;  
int ThisTask;
int NTask;

void io(int num,double x,double y,double vx,double vy,int id,double Cj)
{
 char name[50];
 sprintf(name, "%s_%04d_%d","out", num,ThisTask); 
 FILE *outfile;  
 outfile=fopen(name,"a");
 fprintf(outfile, "%1.16lf\t%1.16lf\t%1.16lf\t%1.16lf\t%i\t%1.16lf\n",x,y,vx,vy,id,Cj);
 fclose(outfile);
}

void iorem(double tout,double x,double y,double vx,double vy,int id,double Cj,int type)
{
 char name[50];
 sprintf(name, "remove/remove_%d",ThisTask); 
 FILE *outfile;
 outfile=fopen(name,"a");
 fprintf(outfile, "\n%lf\t%lf\t%lf\t%e\t%e\t%i\t%e\t%i",tout,x,y,vx,vy,id,Cj,type);
 fclose(outfile);
}

void createFiles(int TaskNum, int outNum)
{
 char name[50];
 int i,j;
 FILE *outfile;  
 for(i=0;i<TaskNum;i++)
 {
  for(j=0;j<=outNum;j++)
  {
   sprintf(name, "%s_%04d_%d","out", j,i); 
   outfile=fopen(name,"wt"); 
   fclose(outfile);
  }
  sprintf(name, "remove/remove_%d",i); 
  outfile=fopen(name,"wt"); 
  fclose(outfile);
 }
}

void acc_part(double *mass,double *y,double *dy, int nbod, int n)
{
 int l;
 double xx,yy,rr2,fac;
 l=4*nbod;
 dy[l]=y[l+2];
 dy[l+1]=y[l+3];  
 xx = y[l] - y[0];
 yy = y[l+1] - y[1];
 rr2 = xx*xx + yy*yy;
 fac = mass[0]/sqrt(rr2)/rr2;  
 dy[l+2]=-fac*xx;
 dy[l+3]=-fac*yy; 
 xx = y[l] - y[4];
 yy = y[l+1] - y[5];
 rr2 = xx*xx + yy*yy;
 fac = mass[1]/sqrt(rr2)/rr2;  
 dy[l+2]-=fac*xx;
 dy[l+3]-=fac*yy;  
}

/*restricted three body problem*/
void acc_bod(double *mass,double *y, double *dy,int nbod)
{
 int i,j;
 double mi,mj,xx,yy,rr2,axx,ayy,fac; 
 for(i=0;i<4*nbod;i+=4)
 {
  dy[i]=y[i+2];
  dy[i+1]=y[i+3];
  dy[i+2]=0.0;
  dy[i+3]=0.0;
 }
 i = 0;
 j = 4;     
 mi = mass[i];
 mj = mass[j/4];
 xx = y[i]-y[j];
 yy = y[i+1]-y[j+1];
 rr2 = xx*xx + yy*yy;
 fac = 1./sqrt(rr2)/rr2;  
 axx = xx*fac;
 ayy = yy*fac;
 dy[i+2] = -axx*mj;
 dy[i+3] = -ayy*mj;
 dy[j+2] = axx*mi;
 dy[j+3] = ayy*mi; 
}

void bs_der(double *mass,double *y,double *dy,int nbod,int n)
{
 acc_bod(mass,y,dy,nbod);
 acc_part(mass,y,dy,nbod,n);
}

void io_thread(double *,int,int,int,int*);

double integral(double *mass,double x,double h0,double *y,int nbod,int n,int I, int *remove,int id)
{ 
 double xa,*dy,*tp;
 tp =(double*)malloc(12*(n+4)*sizeof(double));
 dy =(double*)malloc(4*n*sizeof(double)); 
 int i,ii;
 int lt[11]={0,1,2,3,4,6,8,12,16,24,32};
 double alt[11]={0.0,1.0,2.0,3.0,4.0,6.0,8.0,12.0,16.0,24.0,32.0}; 
 int lbreak,m,l,m1,mmk,i1max,i1,ik,k,idiv;
 double fl,flt,xb,varm,varma,d[6],h,hd,eta2,dta,yb,c,b1,den,b,dtn,var;
 xa=x;
 varma=0;
 for(i=0;i<6;i++)
  d[i]=0;
 for(i=0;i<n;i++)
 {
  dy[i]=0.0;
 }  
 bs_der(mass,y,dy,nbod,n); 
 for(i=0;i<n;i++)
 {
  ii=12*(i+1);
  tp[ii-1]=fabs(y[i]);
  if(tp[ii-1]< eps)
   tp[ii-1]=eps;
  tp[ii-4]=dy[i]; 
  tp[ii]=y[i];
  tp[ii-2]=0;
  tp[ii-3]=0;
  tp[ii-11]=0;    
 }
 for(idiv=0;idiv<=100;idiv++)
 {
  xb=h0+xa; 
  m = 1;
  lbreak=1;
  while((m<=10)&&lbreak)
  { 
   l=lt[m];
   fl=alt[m];
   varm=0.0;
   m1= ((m-1)<6) ? (m-1) : 6;
   if(m1!=0)
   {
    for(k=1;k<=m1;k++)
    {               
     mmk=m-k;
     flt=alt[mmk];   
     d[k-1]=(fl*fl/flt/flt);
    }
   }         
   h=h0/fl;
   hd=0.5*h;
   for(i=0;i<n;i++)
   {
    ii=12*(i+1);
    tp[ii-3]=tp[ii]; 
    y[i]=tp[ii]+hd*tp[ii-4];   
   }
   i1max=2*l-1;
   x=xa;   
   for(i1=0;i1<i1max;i1++)    
   {
    x+=hd;   
    bs_der(mass,y,dy,nbod,n); 
    for(i=0;i<n;i++)       
    {     
     ii=12*(i+1);
     tp[ii-1]=(tp[ii-1]>fabs(y[i]))?tp[ii-1]:fabs(y[i]);
     eta2=tp[ii-3]+h*dy[i];
     tp[ii-3]=y[i];
     y[i]=eta2;       
    }   
   } 
   bs_der(mass,y,dy,nbod,n);
   for(i=0;i<n;i++)
   {    
    ii=12*(i+1);
    dta=tp[ii-11];
    yb=(tp[ii-3]+y[i]+hd*dy[i])/2.0;    
    c=yb;            
    tp[ii-11]=yb;     
    if(m1!=0)
    {
     for(k=1;k<=m1;k++)
     { 
      b1=d[k-1]*dta;
      den=b1-c;
      dtn=dta;
      if(fabs(den)>0.)                    
      {
       b=(c-dta)/den;
       dtn=c*b;       
       c=b1*b;        
      }            
      ik=ii-11+k;
      dta=tp[ik];
      tp[ik]=dtn; 
      yb=yb+dtn; 
     }      
     var=fabs(tp[ii-2]-yb)/tp[ii-1];     
     varm=(varm>var)?varm:var;    
    }
    tp[ii-2]=yb;
   }
   if(m>3)
   {                 
    if(varm<=eps)
    {       
     x=xb;    
     for(i=0;i<n;i++)
     {
      ii=12*(i+1);
      y[i]=tp[ii-2];
     }        
     free(tp);
     free(dy);   
    
     h0=h0*1.5*pow(0.6,(m-1-m1)); 
    
     return x;
    }     
    if(varm>=varma)      
     lbreak =0;              
   }
   varma=varm;
   m++;
  }
  h0=h0/2.0;
 } 
 printf("ERROR integral: lack of convergence !!!!"); 
 *remove=1;
 free(tp);
 free(dy);
 return x;
}
 
int main(int argc, char **argv)
{
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
 MPI_Comm_size(MPI_COMM_WORLD, &NTask);
 MPI_Status status;
 double t0,dt,tstop,dtout,kappa;
 
  //чтение файла с параметрами*/
 FILE *inparfile;
 inparfile=fopen("param.in","r");
 fscanf(inparfile,"%lf %lf %lf %lf %lf",&t0,&tstop,&dt,&dtout,&kappa);  
 fclose(inparfile);
 int Stps;

 double *y,*mass;
 int n,i,nbod,*Id,id,j;
 if(ThisTask==0)
  createFiles(NTask,floor(tstop/dtout));
 double *yh,*yht;
 int *ID;
  //массивы для координат массивных тел и частиц и масс звезд и планет
 yh = (double*)malloc(4*NPLMAX*sizeof(double));
 yht = (double*)malloc(4*NTPMAX*sizeof(double));
 ID = (int*)malloc(NTPMAX*sizeof(int));
 mass = (double*)malloc(NPLMAX*sizeof(double));
 double ms,x1,x2,v1,v2;
 int ntp;
 /*чтение начальных данных*/
 FILE *insfile;
 insfile=fopen("stars.in","rt");
 i=0;
 while((!feof(insfile))&&(i<4*NPLMAX))	
 {
  fscanf(insfile,"%lf %lf %lf %lf %lf\n",&ms,&x1,&x2,&v1,&v2);  
  mass[i/4]=ms;
  yh[i]=x1; 
  yh[i+1]=x2;
  yh[i+2]=v1; 
  yh[i+3]=v2;
  i+=4;
 }
 fclose(insfile);
 nbod=i/4;
  
 /* Get data for the run and the test particles*/
 FILE *intpfile;
 intpfile=fopen("particles.in","rt");
 i=0; 
 while((!feof(intpfile))&&(i<4*NTPMAX))	
 {
  fscanf(intpfile,"%lf %lf %lf %lf %i\n",&x1,&x2,&v1,&v2,&id);
  yht[i]=x1; 
  yht[i+1]=x2;
  yht[i+2]=v1; 
  yht[i+3]=v2;    
  ID[i/4]=id;
  i+=4;       
 } 
 fclose(intpfile); 
 ntp=i/4; 
 int partSize,shift;  
 partSize = ntp/NTask;
 shift = ntp%NTask;        
 int k;
 k=0; 
 if(shift>ThisTask)       
   n=4*(nbod+partSize+1);    
 else
   n=4*(nbod+partSize); 
 y=(double*)malloc(n*sizeof(double));  
 Id = (int*)malloc(n/4*sizeof(int));    
 for(j=0;j<4*nbod;j++)
 {
  y[j] = yh[j]; 
  Id[j/4]=j/4;  
  }
 int *temp;
 temp = malloc(NTask * sizeof(int));
 MPI_Allgather(&n, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
 k=0;
 for(i=0;i<ThisTask;i++)
  k+=(temp[i]-4*nbod);
 while(j<n)
 {    
   y[j] = yht[k];    
   Id[j/4]=ID[k/4];   
   j++;  
   k++;   
  } 
 free(temp);
 double *yr, *ys;
 yr=(double*)malloc(4*(nbod+1)*sizeof(double)); 
 ys=(double*)malloc(4*nbod*sizeof(double)); 
 int *remove,I,num,s; 
 s=4*nbod;
 remove=(int*)malloc(sizeof(int)); 
 double tout,dttmp,x,Rp,Rs,Cj,Cr;
 double Rh,r;
 Rh=pow((mass[1]/3./mass[0]),1./3.);
 Stps=floor(dtout/dt);
 for(I=s;I<n;I+=4)
 {
  for(i=0;i<s;i++)
   yr[i]=y[i];  
  yr[s]=y[I];
  yr[s+1]=y[I+1];
  yr[s+2]=y[I+2];
  yr[s+3]=y[I+3];  
  Rs=sqrt((yr[s]-yr[0])*(yr[s]-yr[0])+(yr[s+1]-yr[1])*(yr[s+1]-yr[1])); 
  Rp=sqrt((yr[s]-yr[4])*(yr[s]-yr[4])+(yr[s+1]-yr[5])*(yr[s+1]-yr[5]));
  if(Rp<0.75*Rh)
   continue;
  v2=yr[s+2]*yr[s+2]+yr[s+3]*yr[s+3];
  Cj=-v2+4*M_PI*(yr[s]*yr[s+3]-yr[s+1]*yr[s+2])+2.*(mass[0]/Rs+mass[1]/Rp);
  Cr=Cj;
  io(0,yr[s],yr[s+1],yr[s+2],yr[s+3],Id[I/4],Cj/4./M_PI/M_PI);  
  *remove=0;
  num=1;
  tout =t0+dtout;    
  while((tout<=tstop)&&(!*remove))
  { 
   j=0;   
   while((j<Stps)&&(!*remove)) 
   { 

    for(i=0;i<s;i++)
      ys[i]=yr[i]; 
    y[I]=yr[s];
    y[I+1]=yr[s+1];
    y[I+2]=yr[s+2];
    y[I+3]=yr[s+3]; 
     for(i=0;i<s;i++)
      yr[i]=ys[i];  
     yr[s]=y[I];
     yr[s+1]=y[I+1];
     yr[s+2]=y[I+2];
     yr[s+3]=y[I+3];      
    dttmp = dt;        
    x = 0.0;
   
    while(((fabs(x-dt)/dt) > 1.0e-7)&&(!*remove))
    {    
     x=integral(mass,x,dttmp,yr,nbod,4*(nbod+1),I,remove,Id[I/4]);
     dttmp = dt - x;        
     Rs=sqrt((yr[s]-yr[0])*(yr[s]-yr[0])+(yr[s+1]-yr[1])*(yr[s+1]-yr[1])); 
     Rp=sqrt((yr[s]-yr[4])*(yr[s]-yr[4])+(yr[s+1]-yr[5])*(yr[s+1]-yr[5]));
     v2=yr[s+2]*yr[s+2]+yr[s+3]*yr[s+3];
     Cr=-v2+4*M_PI*(yr[s]*yr[s+3]-yr[s+1]*yr[s+2])+2.*(mass[0]/Rs+mass[1]/Rp);  
     if(Rp<kappa*Rh)
     {
       iorem(tout-dtout+dt*(j+1)-dttmp,yr[s],yr[s+1],yr[s+2],yr[s+3],Id[I/4],fabs(Cr-Cj)/Cj,0);
       *remove=1; 
     }       
    }                    
    for(i=0;i<s;i++)
      ys[i]=yr[i]; 
    y[I]=yr[s];
    y[I+1]=yr[s+1];
    y[I+2]=yr[s+2];
    y[I+3]=yr[s+3];          
    r=sqrt(yr[s]*yr[s]+yr[s+1]*yr[s+1]);
    if(r>4.)
    {
      iorem(tout-dtout+dt*j,yr[s],yr[s+1],yr[s+2],yr[s+3],Id[I/4],fabs(Cr-Cj)/Cj,1);
      *remove=1;                  
    }   
     j++;   
   }
   if(!*remove)       
    io(num,yr[s],yr[s+1],yr[s+2],yr[s+3],Id[I/4],fabs(Cr-Cj)/Cj);  
   
    num++;
    tout = tout + dtout;       
  }  
 } 
 free(yr); 
 free(mass);
 free(y); 
 free(Id);
 free(remove);
 MPI_Finalize();
}



