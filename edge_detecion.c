#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
typedef struct PGMImage {
    char pgmType[3];
    int **data;
    int row;
    int col;
    int maxValue;
}PGM;
typedef struct GaussFilter {
    double **data;
    int kernel_size;
}GAUSS;
void print_mat(int**,int,int);
PGM* create_pgm(char[],int,int,int);
PGM* read_from_p2(char[]);
PGM* read_from_p5(char[]);
void print_pgm_info(PGM*);
PGM* pgm_type_check(char[]);
GAUSS* create_gaussFilter(int,float);
int** allocateMatrix(int,int);
GAUSS* create_gauss(int);
int** apply_gauss(PGM*,GAUSS*);
void img_to_file(int**,PGM*,char[]);
int** sobel_x(PGM*);
int** sobel_y(PGM*);
int** sobel_xy(PGM*);
int** laplacian(PGM*,int);
void normalize_img(int**,PGM*);

int main(){

	FILE* file;
	PGM* pgm;
	PGM* blur_pgm_yedek;
	PGM* blur_pgm;
	GAUSS* gauss;
	char file_name[20];
	char test_name[30];
	int kernel,i,j;
	int** blurred;
	int** sobel;
	int** laplacian_img;
	
	printf("Lutfen pgm dosya giriniz (file.pgm) : ");
	scanf("%s", file_name);
			
	pgm=pgm_type_check(file_name);
	
	sobel =sobel_x(pgm);
	//normalize_img(sobel,pgm);
	img_to_file(sobel,pgm,"sobel_x.pgm");
	sobel =sobel_y(pgm);
	//normalize_img(sobel,pgm);
	img_to_file(sobel,pgm,"sobel_y.pgm");
	sobel =sobel_xy(pgm);
	//normalize_img(sobel,pgm);
	img_to_file(sobel,pgm,"sobel_xy.pgm");	
	
	
	
	int kernels[] = {3,5,7};
	float sigma[] = {1.0,2.0,4.0};
	char knum[2];
	char snum[2];
	
	char tut[20];
	
	strcpy(tut,file_name);
	
	for (i=0;i<3;i++){
		
		itoa(kernels[i], knum, 10);
		
		for(j=0;j<3;j++){	
		strcpy(file_name,tut);
		strcat(file_name,knum);
		strcat(file_name,"_");
		gcvt(sigma[j], 10, snum);
		strcat(file_name,snum);
		strcat(file_name,".pgm");
		gauss =create_gaussFilter(kernels[i],sigma[j]);
		blurred=apply_gauss(pgm,gauss);
		img_to_file(blurred,pgm,file_name);	
		
		blur_pgm = pgm_type_check(file_name);
		blur_pgm_yedek = pgm_type_check(file_name);
		
		
		sobel =sobel_x(blur_pgm);
		strcpy(file_name,tut);
		strcat(file_name,knum);
		strcat(file_name,"_");
		gcvt(sigma[j], 10, snum);
		strcat(file_name,snum);
		strcat(file_name,"_sobel_blur_x.pgm");
		//normalize_img(sobel,blur_pgm);
		img_to_file(sobel,blur_pgm,file_name);
		
		sobel =sobel_y(blur_pgm);
		strcpy(file_name,tut);
		strcat(file_name,knum);
		strcat(file_name,"_");
		gcvt(sigma[j], 10, snum);
		strcat(file_name,snum);
		strcat(file_name,"_sobel_blur_y.pgm");
		//normalize_img(sobel,blur_pgm);
		img_to_file(sobel,blur_pgm,file_name);
		
		
		sobel =sobel_xy(blur_pgm);
		strcpy(file_name,tut);
		strcat(file_name,knum);
		strcat(file_name,"_");
		gcvt(sigma[j], 10, snum);
		strcat(file_name,snum);
		strcat(file_name,"_sobel_blur_xy.pgm");
		//normalize_img(sobel,blur_pgm);
		img_to_file(sobel,blur_pgm,file_name);
		
		
		
		laplacian_img=laplacian(blur_pgm_yedek,1);
		strcpy(file_name,tut);
		strcat(file_name,knum);
		strcat(file_name,"_");
		gcvt(sigma[j], 10, snum);
		strcat(file_name,snum);
		strcat(file_name,"_laplacian_1.pgm");
		//normalize_img(laplacian_img,blur_pgm_yedek);
		img_to_file(laplacian_img,blur_pgm_yedek,file_name);
		
		
		laplacian_img=laplacian(blur_pgm_yedek,2);
		strcpy(file_name,tut);
		strcat(file_name,knum);
		strcat(file_name,"_");
		gcvt(sigma[j], 10, snum);
		strcat(file_name,snum);
		strcat(file_name,"_laplacian_2.pgm");
		//normalize_img(laplacian_img,blur_pgm_yedek);
		img_to_file(laplacian_img,blur_pgm_yedek,file_name);
		
		}
		
	}
			
	return 0;
}
PGM* pgm_type_check(char file_name[]){
	
	char type[2];
	FILE* file;
	PGM* pgm;
	
	file = fopen(file_name,"r");
	
	if(file == NULL){
		printf("Can't open file to read!");
		return NULL;
		
	}
	
	fscanf(file, "%s\n", &type);
	fclose(file);
	
	if (strcmp(type,"P2")==0){
		
		//printf("P2 Dosyasi\n");
		pgm=read_from_p2(file_name);
		return pgm;
		
	}else if (strcmp(type,"P5")==0){
		
		//printf("P5 Dosyasi\n");
		pgm=read_from_p5(file_name);
		return pgm;
		
	}
	
	
}
PGM* read_from_p2(char file_name[]){
	
	int row,col,maxVal;
	char type[2];
	FILE* file;
	PGM* pgm;
	
	file = fopen(file_name,"r");
	
	if(file == NULL){
		printf("Can't open file to read!");
		return NULL;
		
	}else
	//P2 veya P5 olabilir. Kontrol edilmeli.
	fscanf(file, "%s\n", &type);

	
	//Bos okuma yap ve satiri gec
	fscanf(file, "%*[^\n]\n");
	
	fscanf(file, "%d %d\n", &col,&row);

	fscanf(file, "%d\n", &maxVal);
	
	pgm = create_pgm(type,row,col,maxVal);
	
	int i,j;
	
    for(i=0;i<pgm->row;i++){
    	for(j=0;j<pgm->col;j++){
            if(fscanf(file, "%d", &pgm->data[i][j]) != 1){
            	printf("Dosya bitti.\n");
				return NULL;	
			}
		}
	}
	
	return pgm;
}
PGM* read_from_p5(char file_name[]){
	
	int row,col,maxVal;
	char type[2];
	FILE* file;
	PGM* pgm;
	
	file = fopen(file_name,"rb");
	
	if(file == NULL){
		printf("Can't open file to read!");
		return NULL;
		
	}
	//P2 veya P5 olabilir. Kontrol edilmeli.
	fscanf(file, "%s\n", &type);
	
	//Bos okuma yap ve satiri gec
	fscanf(file, "%*[^\n]\n");
	
	fscanf(file, "%d %d\n", &col,&row);
	
	fscanf(file, "%d", &maxVal);
	
	pgm = create_pgm(type,row,col,maxVal);
	
	int i,j,a;
	
    for(i=0;i<pgm->row;i++){
    	for(j=0;j<pgm->col;j++){
    		
			a= fgetc(file);
            pgm->data[i][j]=a;	
		}
	}
	
	return pgm;
}
PGM* create_pgm(char type[],int row,int col,int maxVal){
	
	int i;
	
	PGM* pgm = (PGM*)malloc(sizeof(PGM));

	pgm->row = row;
	pgm->col = col;
	pgm->maxValue = maxVal;
	strcpy(pgm->pgmType, type);
	
	pgm->data = (int**)malloc((pgm->row)*sizeof(int*));
	
    for(int i = 0; i < pgm->row; i++)
    {
        pgm->data[i] = (int*)malloc((pgm->col)*sizeof(int));
    }
	return pgm;
}
void print_pgm_info(PGM* pgm){
	
	int i,j;
	printf("-------------INFO----------");
	printf("About PGM File:\n");
	
	printf("Type: %s\n",pgm->pgmType);
	printf("Rows and cols: %d  %d\n",pgm->row,pgm->col);
	printf("Max Value: %d\n",pgm->maxValue);
	
	printf("Image: \n ");
	for(i=0;i<pgm->row;i++){
    	for(j=0;j<pgm->col;j++){
            
            printf("%d ",pgm->data[i][j]);
            
		}
		printf("\n");
	}
	
	
}
GAUSS* create_gauss(int kernel){
	
	int i;
	
	GAUSS* gauss = (GAUSS*)malloc(sizeof(GAUSS));

	gauss->kernel_size = kernel;
	
	gauss->data = (double**)malloc((gauss->kernel_size)*sizeof(double*));
	
    for(i = 0; i < kernel; i++)
    {
        gauss->data[i] = (double*)malloc((kernel)*sizeof(double));
    }
	return gauss;
}
GAUSS* create_gaussFilter(int smooth_kernel_size,float sigma){
	
    GAUSS* gauss=create_gauss(smooth_kernel_size);
    double sum = 0;
    int i, j;
    float  K = 1/(2*(3.14)* pow(sigma,2));
    
    //printf("K: %f",K);
    for (i = 0; i < smooth_kernel_size; i++) {
        for (j = 0; j < smooth_kernel_size; j++) {
            double x = i - (smooth_kernel_size - 1) / 2.0;
            double y = j - (smooth_kernel_size - 1) / 2.0;
            gauss->data[i][j] = K * exp(((pow(x, 2) + pow(y, 2)) / ((2 * pow(sigma, 2)))) * (-1));
            sum += gauss->data[i][j];
        }
    }

    for (i = 0; i < smooth_kernel_size; i++) {
        for (j = 0; j < smooth_kernel_size; j++) {
            gauss->data[i][j] /= sum;
        }
    }
    
    for (i = 0; i < smooth_kernel_size; i++) {
        for (j = 0; j < smooth_kernel_size; j++) {
            printf("%f ", gauss->data[i][j]);
        }
        printf("\n");
    }
    return gauss;	
	
}
int** apply_gauss(PGM*pgm,GAUSS* gauss){
	
	int i,j,k,l,i_tut,j_tut;
	double sum;
	int** blur_img= allocateMatrix(pgm->row,pgm->col);
	
	
	memcpy(blur_img, pgm->data, sizeof(pgm->data));
	
	int bound = gauss->kernel_size-1;

	for (i = 0; i < (pgm->row)-bound; i++) {
		
        for (j = 0; j < (pgm->col)-bound; j++) {
        	
        	sum = 0;
			i_tut = i + (gauss->kernel_size-1) /2;
			j_tut = j + (gauss->kernel_size-1) /2;
 
			for (k=0;k<gauss->kernel_size;k++){
        		
        		for(l=0;l<gauss->kernel_size;l++){
        			
        			sum += (pgm->data[i+k][j+l])*(gauss->data[k][l]);
        			//printf("PGM[%d][%d] ile Gauss[%d][%d] yi carptim\n",i+k,j+l,k,l);
				}
			}

			blur_img[i_tut][j_tut] = (int)sum ;
            
        }
    }
	

	return blur_img;
	
}
int** allocateMatrix(int row, int col)
{
    int **matrix;
    int i;
 
    matrix = (int **)malloc(sizeof(int*) * row);
    if (matrix == NULL) {
        printf("memory allocation failure");
        return NULL;
    }
 
    for (i = 0; i < row; ++i) {
        matrix[i] = (int *)malloc(sizeof(int) * col);
        if (matrix[i] == NULL) {
            printf("memory allocation failure");
        	return NULL;
        }
    } 
    return matrix;
}	
void img_to_file(int** mat, PGM* pgm,char file_name[]){
	
	
	FILE* file;
	
	if (strcmp(pgm->pgmType,"P2")==0){
		
		file = fopen(file_name,"w");
		
		fprintf(file,"%s\n",pgm->pgmType);
		fprintf(file,"# Created by Bengi\n");
		fprintf(file,"%d %d\n",pgm->col,pgm->row);
		fprintf(file,"%d\n",pgm->maxValue);
		
		int i,j;
		for(i=0;i<pgm->row;i++){
			
			for(j=0;j<pgm->col;j++){
				
				fprintf(file,"%d ",mat[i][j]);
				
			}

		}
		
		
	}else{
		
		file = fopen(file_name,"wb");
		
		fprintf(file,"%s\n",pgm->pgmType);
		fprintf(file,"# Created by Bengi\n");
		fprintf(file,"%d %d\n",pgm->col,pgm->row);
		fprintf(file,"%d",pgm->maxValue);
		
		int i,j;
		for(i=0;i<pgm->row;i++){
			
			for(j=0;j<pgm->col;j++){
				
				fputc(mat[i][j], file);
				
			}
		}
		
	}
	

	
	fclose(file);
	
}
void print_mat(int**mat,int row,int col){
	
	int i,j;
	
	for(i=0;i<row;i++){
		
		for(j=0;j<col;j++){
			
			printf("%d ",mat[i][j]);
			
		}
		printf("\n");
	}
	
}
int** sobel_x(PGM* pgm){
	
	int sobel_x_filter[][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
	int i,j,k,l,i_tut,j_tut;
	double gx;
	int** sobel_img= allocateMatrix(pgm->row,pgm->col);
	
	memcpy(sobel_img, pgm->data, sizeof(pgm->data));
			
	for (i = 0; i < (pgm->row)-2; i++) {
		
        for (j = 0; j < (pgm->col)-2; j++) {
        	
        	gx = 0;
			i_tut = i + 1;
			j_tut = j + 1;
 
			for (k=0;k<3;k++){
        		
        		for(l=0;l<3;l++){
        			
        			gx += (pgm->data[i+k][j+l])*(sobel_x_filter[k][l]);
        			
				}
			}
			sobel_img[i_tut][j_tut] = (int)gx ;            
        }
    }
	
	return sobel_img;
}
int** sobel_y(PGM* pgm){
	
	int sobel_y_filter[][3]= {{-1,-2,-1},{0,0,0},{1,2,1}};
	
	int i,j,k,l,i_tut,j_tut;
	double gy;
	int** sobel_img= allocateMatrix(pgm->row,pgm->col);
	
	memcpy(sobel_img, pgm->data, sizeof(pgm->data));
			
	for (i = 0; i < (pgm->row)-2; i++) {
		
        for (j = 0; j < (pgm->col)-2; j++) {
        	
        	gy = 0;
			i_tut = i + 1;
			j_tut = j + 1;
 
			for (k=0;k<3;k++){
        		
        		for(l=0;l<3;l++){
        			
        			gy += (pgm->data[i+k][j+l])*(sobel_y_filter[k][l]);
        			
				}
			}
			sobel_img[i_tut][j_tut] = (int)gy ;
        }
    }
	
	return sobel_img;
}
int** sobel_xy(PGM* pgm){
	
	int sobel_x_filter[][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
	int sobel_y_filter[][3]= {{-1,-2,-1},{0,0,0},{1,2,1}};
	
	int i,j,k,l,i_tut,j_tut;
	double gxy,gx,gy;
	int** sobel_img= allocateMatrix(pgm->row,pgm->col);
	
	memcpy(sobel_img, pgm->data, sizeof(pgm->data));
			
	for (i = 0; i < (pgm->row)-2; i++) {
		
        for (j = 0; j < (pgm->col)-2; j++) {
        	
        	gxy = gx= gy=  0;
			i_tut = i + 1;
			j_tut = j + 1;
 
			for (k=0;k<3;k++){
        		
        		for(l=0;l<3;l++){
        			
        			gx += (pgm->data[i+k][j+l])*(sobel_x_filter[k][l]);
        			gy += (pgm->data[i+k][j+l])*(sobel_y_filter[k][l]);
        			gxy = sqrt(gx*gx +gy*gy);
				}
			}
			sobel_img[i_tut][j_tut] = (int)gxy ;
        }
    }
	
	return sobel_img;
}
int** laplacian(PGM* pgm,int index){
	
	int laplacian_1[][3] = {{0,-1,0},{-1,4,-1},{0,-1,0}};
	int laplacian_2[][3]= {{-1,-1,-1},{-1,8,-1},{-1,-1,-1}};
	
	int i,j,k,l,i_tut,j_tut,temp;
	double sum;
	int** laplacian_img= allocateMatrix(pgm->row,pgm->col);
	
	memcpy(laplacian_img, pgm->data, sizeof(pgm->data));
			
	for (i = 0; i < (pgm->row)-2; i++) {
		
        for (j = 0; j < (pgm->col)-2; j++) {
        	
        	sum=  0;
			i_tut = i + 1;
			j_tut = j + 1;
 
			for (k=0;k<3;k++){
        		
        		for(l=0;l<3;l++){
        			
        			
        			if (index==1){
        				
        				temp=laplacian_1[k][l];
        				
					}else{
						temp=laplacian_2[k][l];
					}
        			sum += (pgm->data[i+k][j+l])*(temp);

				}
			}
			laplacian_img[i_tut][j_tut] = (int)sum ;
        }
    }
	
	return laplacian_img;
}
void normalize_img(int** mat,PGM* pgm){
	
	int max = -1, min = 100000;
	int i,j;
	
	for(i=0;i<pgm->row;i++){
		
		for(j=0;j<pgm->col;j++){
			
			if(mat[i][j]>max){
				max= mat[i][j];
			}
			
			if(mat[i][j]<min){
				
				min = mat[i][j];
				
			}	
		}
		
	}
	float tut;
	int pay,payda;
	for(i=0;i<pgm->row;i++){
		
		for(j=0;j<pgm->col;j++){
			
			pay = mat[i][j] - min;
			payda = max-min;
			tut =  pay/(float)payda;
			tut=(int)(tut*255);
			mat[i][j] = tut;
			
		}
		
	}
	
	
}
