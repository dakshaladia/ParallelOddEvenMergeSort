%%cu
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<stdio.h>
#include<stdlib.h>
#define CUDA_ERROR_CHECK

__device__ void prez(int i,int *d,int s,int p,int *r)
{   
    int index=i*(s/p);
    int index2=(i+1)*(s/p);
    int c11=0,c21=0,c12=s/p-1,c22=s/p-1;

    for(int k=0;k<s/p;k++)
    {
        if(d[index+c11]<=d[index2+c21])
        {
            r[k]=d[index+c11];
            c11++;
        }
        else
        {
            r[k]=d[index2+c21];
            c21++;
        }
    }
    __syncthreads();
    for(int k=0;k<s/p;k++)
    {
        if(d[index2+c22]<d[index2+c12])
        {
            r[(2*(s/p))-1-k]=d[index2+c22];
            c22--;
        }
        else
        {
            r[(2*(s/p))-1-k]=d[index+c12];
            c12--;
        }
    }
    __syncthreads();
    
    int rt=0;
        for(int y=0;y<(s/p);y++)
        {
            d[index+y]=r[rt++];
        }
        for(int y=0;y<(s/p);y++)
        {
            d[index2+y]=r[rt++];
        }
    __syncthreads();
}

__global__ void pcm (int *d,int s,int p,int x)
{   
    int i=blockIdx.x;
    int index=i*blockDim.x;
    int len=((s/p)*2);
    int tid=threadIdx.x;
    int size=sizeof(int)*len;

    if(x==0 && ((i*2+1)<p))
    { 
        int len=((s/p)*2);
        int size=sizeof(int)*len;
        int *r=(int *)malloc(sizeof(int)*len);
        prez(i*2,d,s,p,r);
        __syncthreads();
    }
    else if(x==1 && ((i*2+2)<p))
    {
        int len=((s/p)*2);
        int size=sizeof(int)*len;
        int *r=(int *)malloc(sizeof(int)*len);
        prez(i*2+1,d,s,p,r);
        __syncthreads();
    }
}  
   
__device__ void quicksort(int *number,int first,int last)
{
   int i, j, pivot, temp;
   if(first<last)
   {
        pivot=first;
        i=first;
        j=last;
    while(i<j)
    {
    while(number[i]<=number[pivot]&&i<last)
        i++;
        while(number[j]>number[pivot])
        j--;
        if(i<j)
        {
            temp=number[i];
            number[i]=number[j];
            number[j]=temp;
        }
    }

    temp=number[pivot];
    number[pivot]=number[j];
    number[j]=temp;
    quicksort(number,first,j-1);
    quicksort(number,j+1,last);
    }
}

__global__ void qsort(int *d,int s,int p,int *fg)
{
    int i=blockIdx.x;
    int index=i*blockDim.x;
    
    for(int j=0;j<(s/p);j++)
    {
       fg[index+j]=d[index+j];
    }
    
    quicksort(fg,index,index+(s/p)-1);
     
    for(int j=0;j<(s/p);j++)
    {
       d[index+j]=fg[index+j]; 
    }
  
    for(int e=0;e<s/p;e++)
    {
        
        printf("%d ",d[e]);
           
    }
    printf("\n");*/
}

int main()
{   
    int p=2,s=12;
    cudaEvent_t start, stop;     		
	float elapsed_time_ms;  
    int *d_a,*d_r,*d_f;
    int fg[s]={0};
    int d[12]={12,11,10,9,8,7,6,5,4,3,2,1};
    printf("Original Array:\n");
    for(int i=0;i<p*(s/p);i++)
    {
        printf("%d ",d[i]);
    }
    printf("\n\n\n");
    int size=(sizeof(int)*s);
    cudaMalloc((void**)&d_a,size);
    cudaMalloc((void**)&d_f,size);
    int lenn=s/p;
    cudaEventCreate(&start);     		
	cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    cudaMemcpy(d_a,&d,size,cudaMemcpyHostToDevice);
    cudaMemcpy(d_f,&fg,size,cudaMemcpyHostToDevice);
    qsort<<<p,lenn>>>(d_a,s,p,d_f);
    cudaDeviceSynchronize();
    cudaError_t error =cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("CUDA Error1: %s\n", cudaGetErrorString(error));
    }

    cudaMemcpy(&fg,d_f,size,cudaMemcpyDeviceToHost);
    cudaMemcpy(&d,d_a,size,cudaMemcpyDeviceToHost);
    cudaMemcpy(d_a,&d,size,cudaMemcpyHostToDevice);
    printf("QQSorted Array\n");
    for(int e=0;e<s;e++)
    {
        printf("%d ",d[e]);
    }
    printf("\n");
    for(int pp=0;pp<p;pp++)
        pcm<<<p/2,1>>>(d_a,s,p,pp%2);
 
    cudaDeviceSynchronize();
    error =cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("CUDA Error2: %s\n", cudaGetErrorString(error));
    }
    cudaMemcpy(&d,d_a,size,cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
 
    printf("Sorted Array: \n");
    for(int e=0;e<s;e++)
    {
        printf("%d ",d[e]);
    }
    printf("\n");
    cudaEventRecord(stop, 0);     	
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time_ms, start, stop );
    printf("\nTime taken for the entire computation: %f ms.\n", elapsed_time_ms);
    cudaFree(d_a);
}
