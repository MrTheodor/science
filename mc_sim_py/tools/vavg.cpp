#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

char bg_file[255] = "";
char k_file[255] = "";
char xyz_file[255] = "";
char cfg_file[255] = "";
char data_file[255] = "";
char avg_file[255] = "";

//char array_file[255] = "";


class atom
{
        public:
                char type;
                float x,y,z;
                atom():type(0),x(0),y(0),z(0)
                {}
                
} ;

class box
{
        public:
                atom *atoms;
                int N;
                int allN;
                int size;
				int L;

                int K;
                double k_min, k_max, scale;
                char relations;

                int load_config();
                int load_file(int v);
				
				int calc_sff2();
				
                int save_sf(double *funS, double *valK);				
				
				int calc_avgSk();
				

} ;



int main(int argc,char *argv[])
{
        int operation;
        if(argc > 1)
                operation = atoi(argv[1]);
        else
        {
			printf("HOW TO RUN\n");
			printf("PROGRAM type_of_operation [other parameters]\n\n");
			printf("TYPE OF OPERATIONS\n");
			printf("\t1 - calculating S(k) [sin(kr)/(kr)]: PROGRAM 1 xyz_input cfg_input\n");
			printf("\t2 - calculating average S(k): PROGRAM 2 avg_output\n");
        }
		box box1;
        switch(operation)
        {
			case 1:
				//HOW TO RUN
				//a.out 0 xyz_file cfg_file
				if(argc != 4)
				{
						printf("PROGRAM 1 xyz_input cfg_input k_input\n");
						exit(0);
				}
                strcpy(xyz_file,argv[2]);
                strcpy(cfg_file,argv[3]);				
				
				if (!box1.load_config())
				{
					printf("Error loading config file!\n");
					exit(0);
				}					
				if(box1.relations == 'x')
					box1.load_file(1); 
				else
					box1.load_file(0); 
				box1.calc_sff2();

			break;
			case 2:
				
				if(argc != 3)
				{
						printf("PROGRAM 2 S(k)_input\n");
						exit(0);
				}			
                strcpy(data_file,argv[2]);				
				box1.calc_avgSk();
			break;				
			default:
				//printf("Unknow operation!\n");
				exit(0);
        }
        return 1;
}

int box::load_config()
{
	printf("LOADING CONFIGURATION FILE ...\n");
	FILE *ff;
	ff = fopen(cfg_file,"r");
	if(ff)
	{
		fscanf(ff,"%d\n",&K);
		char tmp0;
		fscanf(ff,"%c\n",&(tmp0));
		relations = tmp0;
		double kmin, kmax;
		fscanf(ff,"%lf\n%lf\n",&kmin, &kmax);
		k_min = kmin;
		k_max = kmax;
		int tL;
		fscanf(ff,"%d\n",&tL);
		L = tL;
        float tScale;
        fscanf(ff, "%f\n", &tScale);
        scale = tScale;
		fclose(ff);
		
		if(k_min < 0)
			k_min = 2*M_PI/(double)L;
		if(k_max < 0)
		{	
			k_max = 2*M_PI/sqrt(2);
			k_max/=2.0;
		}
		printf("\tVectors k: kmax = %f; kmin = %f; size of k-array = %d\n", k_max, k_min, K);
		printf("\tSize of box: %d\n", L);
		
		return 1;
	}
	else
		return 0;
}
int box::load_file(int v)
{
	printf("LOADING XYZ FILE...\n");
	FILE *ff;
	ff = fopen(xyz_file,"r");
//	printf("%s\n", ff_name);
	if(ff)
	{
		printf("Reading file ...\n");
		int tmpN;
		fscanf(ff,"%d\n\n",&tmpN);
		allN = N = tmpN;
		printf("\tAll atoms in file: %d\n", N);

               
		int ii=0;
		tmpN = 0;
		while(!feof(ff))
		{
			if(tmpN >= N) {
				printf("Error! Too much records!\n");
				break;
			}
			char c;
			float x,y,z;
			fscanf(ff,"%c %f %f %f\n", &c, &x, &y, &z);
            
			if(!v)
			if (c != relations) continue;
            ii++;
		}
		printf("\tAtoms of type %c: %d\n",relations, ii);

		N = ii;
		atoms = new atom[N];
		fclose(ff);

		ff = fopen(xyz_file,"r");
		fscanf(ff,"%d\n\n",&tmpN);
		
		tmpN = 0;
		while(!feof(ff))
		{
			if(tmpN >= N) {
				break;
			}                
			char c;
			float x,y,z;
			fscanf(ff,"%c %f %f %f\n", &c, &x, &y, &z);
			
			if(!v)
			if (c != relations) continue;
			atoms[tmpN].type = c;
			atoms[tmpN].x = x * scale;
			atoms[tmpN].y = y * scale;
			atoms[tmpN].z = z * scale;

        	//printf("%c %f %f %f\n", atoms[tmpN].type, atoms[tmpN].x, atoms[tmpN].y, atoms[tmpN].z);
			tmpN++;
		}
		fclose(ff);

		if(tmpN != N) {
			printf("Error! Readed %d 'atoms'\n", tmpN);
			N = tmpN;
			exit(0);
		}
		return N;

	}
	else
	{
		printf("Error! Can not open file: %s!\n", xyz_file);
		return 0;
	}
}
int box::calc_sff2()
{
	printf("CALCULATING S(k)...\n");
	double *funS;
    double *valK;
	

	valK = new double[K];
	
	double k = k_min;
	double dk = (k_max - k_min)/(double)K;
	printf("\tSize of k-array : %d\n",K);
	for(int i=0; i<K; i++)
	{
		valK[i] = k;
		k += dk;	
		//printf("%f\n",valK[i]);
	}	
	
	funS = new double [K];
    for(int i=0; i<K; i++)
        funS[i] = 0;
	
	float dX, dY, dZ;         					//distance between two atoms
	float halfbox;
	halfbox = L/2.0;
	printf("half of box: %f\n",halfbox);
	double rmn;
	int irmn;
	
		
	for (int i = 0; i < N-1; ++i)
	for (int j = i+1; j < N; ++j) 
	{

		dX = fabs(atoms[j].x-atoms[i].x);
		if(dX>=halfbox) dX = L-dX;
		dY = fabs(atoms[j].y-atoms[i].y);
		if(dY>=halfbox) dY = L-dY;
		dZ = fabs(atoms[j].z-atoms[i].z);
		if(dZ>=halfbox) dZ = L-dZ;

		rmn = sqrt(dX*dX+dY*dY+dZ*dZ);
		//printf("%f\n",rmn);
		double tmp_x;
		if(rmn <= halfbox)
		for(int l = 0; l < K; ++l) 
		{
			//printf("%f\n",k);
			tmp_x = valK[l]*rmn;
			funS[l] += 2.0*sin(tmp_x)/tmp_x;
			//printf("%d %f\n", l, funS[l]);
		}

		//printf("%d\n",i);
	}    
	for (int i = 0; i < K; ++i) 
	{
		funS[i] = funS[i]/double(N) + 1;
		//printf("%lf %lf\n",valK[i],funS[i]);
	}
	
		
	//
	if(!save_sf(funS,valK))
		printf("Error with writing S(k) to file!\n");
	
	delete []valK;
	delete []funS;
	return 1;	
}
int box::save_sf(double *funS, double *valK)
{
	FILE *ff;
	char result_file[255];
	strcpy(result_file,xyz_file);

	result_file[strlen(result_file)-1] = 'f';
	result_file[strlen(result_file)-2] = 'f';
	result_file[strlen(result_file)-3] = 's';

	ff = fopen(result_file,"w");
	if(!ff) return 0;

	fprintf(ff,"#\n");

	for(int i=0; i<K; i++)
	{
		fprintf(ff,"%lf %lf\n",valK[i], funS[i]);
	}
	fclose(ff);
	return 1;  
}

int box::calc_avgSk()
{      
	printf("CALCULATING AVERAGE VALUE S(k)...\n");
	char line[255];
	double *sumSk;
	double *sumk;
	unsigned int count = 0;
	unsigned int cc;
	unsigned int n;
	double tmp1, tmp2;

	FILE *ff = fopen(data_file,"r");
	if(!ff)
	{
			printf("Error! Can not open file: %s!\n", data_file);
			exit(0);			
	}	
	fgets(line,255,ff);
	while(!feof(ff))
	{
		fgets(line,255,ff);
		if(line[0] != '#') count++;
		else break;
	}
	fclose(ff);
	
	sumSk = new double[count];
	sumk = new double[count]; 
	for(unsigned int i=0; i<count; i++)
		sumSk[i] = sumk[i] = 0;


	ff = fopen(data_file,"r");
	if(!ff)
	{
			printf("Error! Can not open file: %s!\n", data_file);
			exit(0);			
	}	
	cc = 0;
	n=0;
	while(!feof(ff))
	{
		if(cc == 0)
		{
			fgets(line,255,ff);
		}
		else
		{
			fscanf(ff,"%lf %lf\n",&tmp1, &tmp2);
			sumk[cc-1] += tmp1;
			sumSk[cc-1] += tmp2;
		}
		cc++;
		if(cc == count+1)
		{
			cc = 0;
			n++;
		}
	}
	fclose(ff);

	if(!n) n = 1;

	for(unsigned int i=0; i<count; i++)
	{
		sumSk[i] /= (double)n;
		sumk[i] /= (double)n;
	}
	
	char result_file[255];
	strcpy(result_file,data_file);

	result_file[strlen(result_file)-1] = 'g';
	result_file[strlen(result_file)-2] = 'v';
	result_file[strlen(result_file)-3] = 'a';
	
	ff = fopen(result_file,"w");
	if(!ff)
	{
			printf("Error! Can not create file: %s!\n", result_file);
			exit(0);			
	}
	
	printf("SAVING AVERAGE VALUE S(k)...\n");
	fprintf(ff,"#k S(k)\n");		
	for(unsigned int i=0; i<count; i++)
	{
		fprintf(ff,"%lf %lf\n",sumk[i],sumSk[i]);
	}
	fclose(ff);

	
	delete []sumSk;
	delete []sumk;
	printf("DONE\n");
	return 1;

}
