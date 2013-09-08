/*
PIOTR KNYCHA?A - UAM/POZNA? 2008

WYZNACZANIE CZYNNIKA STUKTURY DLA UK?AD?W ANIZOTROPOWYCH

WERSJA POPRAWNA - TEST 16 X 2009


WZ?R:
S(k) = 1/N [ (suma(sin k ri))^2 + (suma(cos k ri))^2 ]
	sin(k ri) oraz cos(k ri) - rozbijamy z wzoru na iloczyn skalarny 3D - sin(kx rx + ky ry + kz rz)
	dla tych wzor?w stosujemy wz?r na sin/cos sumy (stosujemy go dwa razy)
	
	sin(0 rx) = 0
	cos(0 rx) = 1
	
	wyznaczamy
	sin(k0 rx)
	cos(k0 rx)
	
	kolejne warto?ci sin/cos to wielokrotno?ci n * (k0 rx) - obliczamy je rekurencyjnie z wzoru na sin sumy	
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

//char bg_file[255] = "";	
//char k_file[255] = "";
char xyz_file[255] = "";	//buffer for file name
char cfg_file[255] = "";
char data_file[255] = "";
float scale_factor;
//char avg_file[255] = "";

//char array_file[255] = "";
inline double get_time(void) 
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + 1e-6 * (double)tv.tv_usec;
}

class point
{
	public:
		double x,y,z;
		point():x(0.0),y(0.0),z(0.0)
		{}
};

class atom
{
	public:
		char type;	//type of atom i
		point pos;	//position of atom i
		
		point *tSIN;//array of sin/cos from -kmax,...,k0,...,+kmax 	
		point *tCOS;//k0 is ath position Z
		atom()
		{
			type = 0;
			pos.x = pos.y = pos.z = 0;
		}
		~atom()
		{
			delete []tSIN;
			delete []tCOS;
		}
		int setXYZ(double x, double y, double z)
		{
			pos.x = x;
			pos.y = y;
			pos.z = z;
			
			return 1;
		}
		int setTYPE(char type)
		{
			this->type = type;
			return 1;
		}
		//calculating sin(nk*) and cos(nk*) for all vectors k
		int setSINCOS(int Z, int nZ, double k0);
		int getSINCOS(int Z, int nZ, double k0);
		int getSINCOS2(int kn, double k0, int Z);
		//calculating sin of the scalar product of vectors k and r
		double inline getSSIN(int kx, int ky, int kz);
		//calculating cos of the scalar product of vectors k and r					
		double inline getSCOS(int kx, int ky, int kz);
} ;
int atom::setSINCOS(int Z, int nZ, double k0)
{
	tSIN = new point[nZ];
	tCOS = new point[nZ];
	//sin(0)=0; cos(0)=1
	//sin[Z]=0.0; cos[Z]=1.0
	tSIN[Z].x = tSIN[Z].y = tSIN[Z].z = 0.0;
	tCOS[Z].x = tCOS[Z].y = tCOS[Z].z = 1.0;
	//sin(k*)=sin(k0*rx); cos(k*)=cos(k0*rx)
	//sin[Z-1]=-sin[Z+1]; cos[Z-1]=cos[Z+1]
	tSIN[Z+1].x = sin(k0*pos.x); tSIN[Z-1].x = -tSIN[Z+1].x;
	tSIN[Z+1].y = sin(k0*pos.y); tSIN[Z-1].y = -tSIN[Z+1].y;
	tSIN[Z+1].z = sin(k0*pos.z); tSIN[Z-1].z = -tSIN[Z+1].z;	
	tCOS[Z+1].x = cos(k0*pos.x); tCOS[Z-1].x = tCOS[Z+1].x;
	tCOS[Z+1].y = cos(k0*pos.y); tCOS[Z-1].y = tCOS[Z+1].y;
	tCOS[Z+1].z = cos(k0*pos.z); tCOS[Z-1].z = tCOS[Z+1].z;	
	
	
	double tmpX, tmpY, tmpZ;
	for(int i=2; i<Z+1; i++)
	{
		//here we use sin(x+y) and cos(x+y) to calculate value for k* > k0
		//sin(n*k0), where n>1
		//tSIN[Z+i]
		tmpX = tSIN[Z+i-1].x*tCOS[Z+1].x+tSIN[Z+1].x*tCOS[Z+i-1].x;
		tmpY = tSIN[Z+i-1].y*tCOS[Z+1].y+tSIN[Z+1].y*tCOS[Z+i-1].y;		
		tmpZ = tSIN[Z+i-1].z*tCOS[Z+1].z+tSIN[Z+1].z*tCOS[Z+i-1].z;
		tSIN[Z+i].x = tmpX; tSIN[Z-i].x = -tmpX;
		tSIN[Z+i].y = tmpY; tSIN[Z-i].y = -tmpY;
		tSIN[Z+i].z = tmpZ; tSIN[Z-i].z = -tmpZ;	
		//cos(n*k0), where n>1
		//tCOS[Z+i]
		tmpX = tCOS[Z+i-1].x*tCOS[Z+1].x-tSIN[Z+i-1].x*tSIN[Z+1].x;
		tmpY = tCOS[Z+i-1].y*tCOS[Z+1].y-tSIN[Z+i-1].y*tSIN[Z+1].y;
		tmpZ = tCOS[Z+i-1].z*tCOS[Z+1].z-tSIN[Z+i-1].z*tSIN[Z+1].z;
		tCOS[Z+i].x = tCOS[Z-i].x = tmpX;
		tCOS[Z+i].y = tCOS[Z-i].y = tmpY;
		tCOS[Z+i].z = tCOS[Z-i].z = tmpZ;
	}
	
	return 1;
}

int atom::getSINCOS(int Z, int nZ, double k0)
{
	
	printf("Atom position: %lf %lf %lf\n",pos.x,pos.y,pos.z);
	printf("sin xyz & cos xyz\n");
	for(int i=Z; i<nZ; i++)
		printf("%d (%d) %lf %lf %lf\t%lf %lf %lf\n",i-Z,i,tSIN[i].x,tSIN[i].y,tSIN[i].z,tCOS[i].x,tCOS[i].y,tCOS[i].z);
		
	return 1;
}
int atom::getSINCOS2(int kn, double k0,int Z)
{
	printf("(%lf %lf %lf)%d %lf: %lf %lf %lf, %lf %lf %lf\n",pos.x,pos.y,pos.z,kn-Z,(kn-Z)*k0,tSIN[kn].x,tSIN[kn].y,tSIN[kn].z,tCOS[kn].x,tCOS[kn].y,tCOS[kn].z);
	return 1;	
}

double atom::getSSIN(int kx, int ky, int kz)
{
	//sin of dot product vector k* and vector ri
	double tmp;
	tmp = tSIN[kx].x*tCOS[ky].y*tCOS[kz].z +
		  tCOS[kx].x*tSIN[ky].y*tCOS[kz].z +
		  tCOS[kx].x*tCOS[ky].y*tSIN[kz].z -
		  tSIN[kx].x*tSIN[ky].y*tSIN[kz].z;
	return tmp;
}
double atom::getSCOS(int kx, int ky, int kz)
{
	//cos of dot product vector k* and vector ri
	double tmp;
	tmp = tCOS[kx].x*tCOS[ky].y*tCOS[kz].z -
		  tSIN[kx].x*tSIN[ky].y*tCOS[kz].z -
		  tSIN[kx].x*tCOS[ky].y*tSIN[kz].z -
		  tCOS[kx].x*tSIN[ky].y*tSIN[kz].z;
	return tmp;	
}


class box
{
	public:
		atom *atoms;
		int allN;	//amount of all atoms in the box
		char rel;	//type relations
		int N;		//amount of the atoms of type rel in the box
		int L;		//size of box (box must be cube)
		int Z;		//range of k vector (Z=L/a) (a=sqrt(2) for FCC)
		int nZ;		//range of k vector nZ = 2*Z+1, Z=L/a
		double k0;	//value of k*=4PI/L
		
		box()
		{
			atoms = NULL;
		}
		~box()
		{
			if(atoms)
				delete []atoms;
		}


		int load_config();
		int load_file(int v);
		
		int calc_Sk();
		int static calc_avgSk();		

		
	
		

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
			printf("\t1 - calculating S(k): ROGRAM 1 xyz_input cfg_input\n");
			printf("\t2 - calculating avg S(k): PROGRAM 2 S(k)_input\n");			

        }

		box box1;			
        switch(operation)
        {
			case 1:
				//HOW TO RUN
				//a.out 1 xyz_file cfg_file
				if(argc != 4)
				{
						printf("PROGRAM 1 xyz_input cfg_input\n");
						exit(0);
				}
			
                strcpy(xyz_file,argv[2]);
                strcpy(cfg_file,argv[3]);				
				
				if (!box1.load_config())
				{
					printf("Error loading config file!\n");
					exit(0);
				}					
				box1.load_file(0); 
				box1.calc_Sk();

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
				printf("Error!\n");;
			break;
		}		
        return 1;
}

/*
configuration file:
type of relations (rel)
value of Z (if -1 then Z=L/a)
size of box (L)	
scale factor for jmol file
*/
int box::load_config()
{
	printf("LOADING CONFIGURATION FILE ...\n");
	FILE *ff;
	ff = fopen(cfg_file,"r");
	if(ff)
	{
		char tmp0;
		fscanf(ff,"%c\n",&(tmp0));
		rel = tmp0;
		int tmpZ;
		fscanf(ff,"%d\n",&tmpZ);
		Z = tmpZ;
		int tmpL;
		fscanf(ff,"%d\n",&tmpL);
		L = tmpL;
        fscanf(ff, "%f\n", &scale_factor);
		fclose(ff);
		
		tmpZ = int(double(L)/sqrt(2));
		//tmpZ = double(L);
		//***
		tmpZ /= 2;
		if(Z < 0)
			Z = tmpZ;
		if(Z>tmpZ)
			Z = tmpZ;
		nZ = 2*Z+1;
		k0 = 2*M_PI/double(L);
		//k0=M_PI/double(L);
		
		printf("\tType of relations: %c\n", rel);
		printf("\tValue of Z (and nZ): %d (%d)\n", Z, nZ);
		printf("\tSize of box: %d\n", L);
		printf("\tValue of k*: %lf\n",k0);
        printf("\tScale factor: %f\n", scale_factor);
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
//		printf("Reading file ...\n");
		int tmpN;
		fscanf(ff,"%d\n\n",&tmpN);
		allN = N = tmpN;
		printf("\tAll atoms in file: %d\n", N);

		int tmpX=0, tmpY=0, tmpZ=0;                
		int ii=0;
		tmpN = 0;
		while(!feof(ff))
		{
			if(tmpN >= N) {
				printf("Error! Too much records!\n");
				break;
			}
			char c;
			float fx,fy,fz;
            int x,y,z;
			fscanf(ff,"%c %f %f %f\n", &c, &fx, &fy, &fz);
            x = scale_factor*fx;
            y = scale_factor*fy;
            z = scale_factor*fz;

			if(x > tmpX) tmpX = x;
			if(y > tmpY) tmpY = y;
			if(z > tmpZ) tmpZ = z;
            
			if(!v)
			if (c != rel) continue;
            ii++;
		}
		fclose(ff);
		printf("\tAtoms of type %c: %d\n",rel, ii);

		N = ii;
		atoms = new atom[N];

		ff = fopen(xyz_file,"r");
		fscanf(ff,"%d\n\n",&tmpN);
		
		tmpN = 0;
		while(!feof(ff))
		{
			if(tmpN >= N) {
				break;
			}                
			char c;
			int fx,fy,fz;
            int x,y,z;
			fscanf(ff,"%c %d %d %d\n", &c, &fx, &fy, &fz);
		    
            x = fx*scale_factor;
            y = fy*scale_factor;
            z = fz*scale_factor;

			if(!v)
			if (c != rel) continue;
		
			atoms[tmpN].setTYPE(c);
			//atoms[tmpN].setXYZ(double(x),double(y),double(z));
			
			atoms[tmpN].setXYZ(double(x-L/2.0),double(y-L/2.0),double(z-L/2.0));
			//atoms[tmpN].setXYZ(double(x-tmpX/2.0),double(y-tmpY/2.0),double(z-tmpZ/2.0));
			//atoms[tmpN].setXYZ(double(x/double(tmpX)),double(y/double(tmpY)),double(z/double(tmpZ)));
			//atoms[tmpN].setXYZ(double(x/double(tmpX)-0.5),double(y/double(tmpY)-0.5),double(z/double(tmpZ)-0.5));
        	//printf("%c %lf %lf %lf\n", atoms[tmpN].type, atoms[tmpN].pos.x, atoms[tmpN].pos.y, atoms[tmpN].pos.z);
			tmpN++;
		}
		fclose(ff);
		
		if(tmpN != N) {
			printf("Error! Readed %d 'atoms'\n", tmpN);
			N = tmpN;
			exit(0);
		}
		int sizeX = tmpX + 1;
		int sizeY = tmpY + 1;
		int sizeZ = tmpZ + 1;
		
		int size = (sizeX + sizeY + sizeZ)/3;
		if(sizeX != size) 
		{
			printf("Error! The box must be cube!\n");
//			exit(0);
		}
		
		if(size != L)
		{
			printf("Error! Wrong size of box!");
//			exit(0);
		}		


		return N;

	}
	else
	{
		printf("Error! Can not open file: %s!\n", xyz_file);
		return 0;
	}
}

int box::calc_Sk()
{	
	printf("CALCULATING sin & cos FOR ALL ATOMS...\n");	
	for(int i=0; i<N; i++)
		atoms[i].setSINCOS(Z, nZ, k0);	
	
	printf("CALCULATING S(k)...\n");	
	double ***valS;
	valS = new double**[nZ];
	for(int i=0; i<nZ; i++)
		valS[i] = new double*[nZ];
	for(int i=0; i<nZ; i++)
	for(int j=0; j<nZ; j++)
		valS[i][j] = new double[nZ];
	
	for(int nx=0; nx<nZ; nx++)
	for(int ny=0; ny<nZ; ny++)
	for(int nz=0; nz<nZ; nz++)
	{
		double Ssin=0, Scos=0;
		for(int i=0; i<N; i++)
		{
				Ssin += atoms[i].getSSIN(nx,ny,nz);
				Scos += atoms[i].getSCOS(nx,ny,nz);
		}
		valS[nx][ny][nz] = (Ssin*Ssin + Scos*Scos)/double(N);
		//printf("%d %d %d\n",nx,ny,nz);
	}
	
	printf("SAVING S(k)...\n");	
	FILE *ff;
	char result_file[255];
/*	strcpy(result_file,xyz_file);

	result_file[strlen(result_file)-1] = 'd';
	result_file[strlen(result_file)-2] = '3';
	result_file[strlen(result_file)-3] = 'S';

	ff = fopen(result_file,"w");
	if(!ff) 
	{
		printf("Error! Can not create file: %s!\n", result_file);
		exit(0);
	}
	for(int nx=0; nx<nZ; nx++)
	for(int ny=0; ny<nZ; ny++)
	for(int nz=0; nz<nZ; nz++)
		fprintf(ff,"%d %d %d %lf\n",nx-Z,ny-Z,nz-Z,valS[nx][ny][nz]);
	fclose(ff);
	
	
	printf("CALCULATING AVERAGE S(k)...\n");	*/
	strcpy(result_file,xyz_file);

	result_file[strlen(result_file)-1] = 'f';
	result_file[strlen(result_file)-2] = 'f';
	result_file[strlen(result_file)-3] = 's';

	ff = fopen(result_file,"w");	
	//AVG BLOCK
	int avgN = 3*Z*Z+1;
	double *ss;
	int *cc;
	ss = new double[avgN];	
	cc = new int[avgN];
	for(int i=0; i<avgN; i++)
	{
		ss[i] = 0.0;
		cc[i] = 0;
	}	
	int dd;	
	for(int nx=0; nx<nZ; nx++)
	for(int ny=0; ny<nZ; ny++)
	for(int nz=0; nz<nZ; nz++)		
	{	
		dd = (nx-Z)*(nx-Z)+(ny-Z)*(ny-Z)+(nz-Z)*(nz-Z);
		cc[dd]++;
		ss[dd]+=valS[nx][ny][nz];
	}
	fprintf(ff,"#\n");
	for(int i=1; i<avgN; i++)
	{
		if(cc[i])
			fprintf(ff,"%lf %lf\n",sqrt(i)*k0,ss[i]/double(cc[i]));		
	}
	/*
	int dd;
	for(int i=0; i<avgN; i++)
	{
		double ss=0.0;
		int cc=0;
		for(int nx=0; nx<nZ; nx++)
		for(int ny=0; ny<nZ; ny++)
		for(int nz=0; nz<nZ; nz++)		
		{
			dd = (nx-Z)*(nx-Z)+(ny-Z)*(ny-Z)+(nz-Z)*(nz-Z);
			if(dd == i)
			{
				//printf("%d %d %d, %d - %lf\n",nx-Z,ny-Z,nz-Z,dd,valS[nx][ny][nz]);
				ss = ss + valS[nx][ny][nz];
				cc++;
			}
		}		
		if(cc)
		fprintf(ff,"%lf %lf\n",sqrt(i)*k0,ss/double(cc));
		//printf("\n\n");
	}
	*/
	fclose(ff);
	
	delete []ss;
	delete []cc;
	
	for(int i=0; i<nZ; i++)
	for(int j=0; j<nZ; j++)
		delete []valS[i][j];
	for(int i=0; i<nZ; i++)
		delete []valS[i];
	delete []valS;
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


	/*
	printf("AVG\n");
	double *avgSk;
	double *davgSk;
	int *ccSk;
	int avgN = 3*Z*Z+1;
	avgSk = new double [avgN];
	davgSk = new double [avgN];	
	ccSk= new int [avgN];
	for(int i=0; i<avgN; i++)
	{
		avgSk[i] = 0.0;
		davgSk[i] = 0.0;
		ccSk[i] = 0;
	}
	int dd;
	for(int nx=0; nx<nZ; nx++)
	for(int ny=0; ny<nZ; ny++)
	for(int nz=0; nz<nZ; nz++)
	{	
		//if((nx==ny)&&(ny==nz))
		{
			dd = (nx-Z)*(nx-Z)+(ny-Z)*(ny-Z)+(nz-Z)*(nz-Z);
			if(dd >= avgN)
				printf("\tWrong siez of table avgSk\n");
			if(valS[nx][ny][nz])
			ccSk[dd]++;
			avgSk[dd] += valS[nx][ny][nz];
		}
	}
	for(int i=0; i<avgN; i++)
	{
		if(ccSk[i])
			avgSk[i] /= double(ccSk[i]);
	}
	/*
	for(int nx=0; nx<nZ; nx++)
	for(int ny=0; ny<nZ; ny++)
	for(int nz=0; nz<nZ; nz++)
	{	
		//if((nx==ny)&&(ny==nz))
		{
			dd = (nx-Z)*(nx-Z)+(ny-Z)*(ny-Z)+(nz-Z)*(nz-Z);
			if(dd >= avgN)
				printf("\tWrong siez of table avgSk\n");
			double tmp = fabs(valS[nx][ny][nz]-avgSk[dd]);
			davgSk[dd] += tmp;
		}
	}	
	*/
	
	/*
	for(int i=0; i<avgN; i++)
	{
		if(avgSk[i])
		{
				printf("%lf %lf\n",sqrt(i), avgSk[i]);
		}
			//if(davgSk[i])
			//printf("%lf %lf\n",sqrt(i)*k0,davgSk[i]/double(ccSk[i]));		
	}
	/*
	printf("k*\n");	
	for(int nx=0; nx<nZ; nx++)
	{
	for(int ny=0; ny<nZ; ny++)
	//for(int nz=0; nz<nZ; nz++)
	{
		printf("%d %d %lf\n",nx-Z,ny-Z,valS[nx][ny][Z+1]); 
		
	}
	printf("\n");
	}
	printf("\n\nZk*\n");
	for(int nx=0; nx<nZ; nx++)
	{
	for(int ny=0; ny<nZ; ny++)
	//for(int nz=0; nz<nZ; nz++)
	{
		printf("%d %d %lf\n",nx-Z,ny-Z,valS[nx][ny][nZ-1]); 
		
	}
	printf("\n");
	}	
	printf("\n\nZ/2k*\n");
	for(int nx=0; nx<nZ; nx++)
	{
	for(int ny=0; ny<nZ; ny++)
	//for(int nz=0; nz<nZ; nz++)
	{
		printf("%d %d %lf\n",nx-Z,ny-Z,valS[nx][ny][Z+Z/2]); 
		
	}
	printf("\n");
	}	
	*/
	//delete []avgSk;	
	//delete []davgSk;	
	//delete []ccSk;
	
	
	//CALCULATING MOMENT OF INTERIA MATRIX
	/*
	double matrixI[3][3];
	for(int mi=0;mi<3;mi++)
	for(int mj=0;mj<3;mj++)		
		matrixI[mi][mj]=0;

	double tmp;
	for(int nx=0; nx<nZ; nx++)
	for(int ny=0; ny<nZ; ny++)
	for(int nz=0; nz<nZ; nz++)
	{	
		matrixI[0][0] += valS[nx][ny][nz]*((ny-Z)*(ny-Z)*k0*k0+(nz-Z)*(nz-Z)*k0*k0);
		matrixI[1][1] += valS[nx][ny][nz]*((nx-Z)*(nx-Z)*k0*k0+(nz-Z)*(nz-Z)*k0*k0);
		matrixI[2][2] += valS[nx][ny][nz]*((nx-Z)*(nx-Z)*k0*k0+(ny-Z)*(ny-Z)*k0*k0);		
		tmp = valS[nx][ny][nz]*(nx-Z)*(ny-Z)*k0*k0;
		matrixI[0][1] += tmp;
		matrixI[1][0] += tmp;
		tmp = valS[nx][ny][nz]*(nx-Z)*(nz-Z)*k0*k0;
		matrixI[0][2] += tmp;
		matrixI[2][0] += tmp;
		tmp = valS[nx][ny][nz]*(ny-Z)*(nz-Z)*k0*k0;
		matrixI[1][2] += tmp;
		matrixI[2][1] += tmp;				
	}
	matrixI[0][1] *= -1;	matrixI[1][0] *= -1;
	matrixI[0][2] *= -1;	matrixI[2][0] *= -1;
	matrixI[1][2] *= -1;	matrixI[2][1] *= -1;
	
	for(int mi=0;mi<3;mi++)
	{
		for(int mj=0;mj<3;mj++)		
			printf("%lf ",matrixI[mi][mj]);;	
		printf("\n");
	}
	*/


